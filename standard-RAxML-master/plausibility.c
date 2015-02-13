/*  RAxML-HPC, a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright March 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *
 *  When publishing work that is based on the results from RAxML please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#ifndef WIN32  
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include "axml.h"
#include "rmq.h" //include range minimum queries for fast plausibility checker
#include "bipartitionList.h"

#ifdef __SIM_SSE3

#include <xmmintrin.h>
#include <pmmintrin.h>

#endif

#ifdef _USE_PTHREADS
#include <pthread.h>
#endif

#ifdef _WAYNE_MPI
#include <mpi.h>
extern int processID;
extern int processes;
#endif

#define _NEW_MRE

extern FILE *INFILE;
extern char run_id[128];
extern char workdir[1024];
extern char bootStrapFile[1024];
extern char tree_file[1024];
extern char infoFileName[1024];
extern char resultFileName[1024];
extern char verboseSplitsFileName[1024];

extern double masterTime;

extern const unsigned int mask32[32];

extern volatile branchInfo      **branchInfos;
extern volatile int NumberOfThreads;
extern volatile int NumberOfJobs;


#define _USE_FAST_PLAUSIBILITY

#ifdef _USE_FAST_PLAUSIBILITY

/************************************  new fast code by David Dao *****************************************************************************/

/* Calculates the size of every subtree and extract its bipartition by seperating the subtree from the small tree 
 This Algorithm works on the induced bifurcating subtree only and needs therefore no multifurcation adaption 
 It then counts, how many bipartition is shared with the reference small tree*/
static int rec_findBipartitions(unsigned int ** bitvectors, int* seq, int arraysize, int* translate, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction)
{
  int 
    i, 
    j,
    o, 
    taxon,
    found = 0,
    before;
    
  unsigned int 
    k, 
    *toInsert = (unsigned int*)NULL,
    *bipartition = (unsigned int*)NULL;
    
  hashNumberType 
    position;
  
  int 
    numberOfSplits = 0, /* stop after extracting n-3 splits */ 
    firstTaxon = taxonToReduction[first - 1] + 1,   
    *V = (int *)rax_malloc((arraysize) * sizeof(int)),
    *bipartitions = (int *)rax_malloc((ntips - 3) * sizeof(int));

  for(i = arraysize - 1; i >= 0; i--) 
    {
      V[i] = 1;
    
      /* Extract Bipartiton from inner node! */
      if(!isTip(translate[seq[i]],numsp) && (numberOfSplits < (ntips - 3)))
	{
	  /* we can be sure of bifurcating trees, therefore j < deg(i) = 2 */
	  for(j = 0; j < 2; j++) 	    
	    V[i] = V[i] + V[i + V[i]];	   

	  /* find out the memory efficient index for this taxon which lies between 0.. (ntips - 1) instead of 1..mxtips  
	     We have to subtract 1 because array starts from 0.. */  
	  int 
	    index = taxonToReduction[translate[seq[i]] - 1];

	  /* save bipartition */
	  bipartitions[numberOfSplits] = index;

	  toInsert = bitvectors[index];

	  //printf("We are working with %i \n", index);
	  /* Calculate Bipartition */
	  for(j = 1; j < V[i]; j++) 
	    {   
	      /* look if Preorderlabel at index i + j is a tip! */
	      if(isTip(translate[seq[i + j]],numsp))
		{
		  /* translate taxon number to a number between 0 and smalltreesize */      
		  taxon = taxonToReduction[translate[seq[i + j]] - 1] + 1;

		  /* set bit to one */
		  toInsert[(taxon-1)  / MASK_LENGTH] |= mask32[(taxon-1) % MASK_LENGTH];
		}
	      else
		{      
		  before = taxonToReduction[translate[seq[i+j]] - 1];

		  bipartition = bitvectors[before];

		  for(k = 0; k < vLength; k++) 
		    toInsert[k] |= bipartition[k];

		  /* jump to the next relevant taxon */ 
		  j = j + V[i+j] - 1;
		}        
	    }
	  /* count number of splits and stop at n-3 */
	  numberOfSplits += 1;
	}
    }
  
  for(i=0;i < (ntips - 3); i++) 
    {
      toInsert = bitvectors[bipartitions[i]];

      /* if bitvector contains first taxon, use its complement */
      if(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]) 
	{                 
	  /* Padding the last bits! */
	  if(ntips % MASK_LENGTH != 0) 
	    {            
	      for(o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--)         
		toInsert[vLength - 1] |= mask32[o-1];       
	    }
	  
	  for(k=0;k < vLength; k++)         
	    toInsert[k] = ~toInsert[k];         
	}
      
      assert(!(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]));
      
      /* compute hash */
      position = oat_hash((unsigned char *)toInsert, sizeof(unsigned int) * vLength);
      
      position = position % hash->tableSize;
      
      /* re-hash to the new hash table that contains the bips of the large tree, pruned down 
	 to the taxa contained in the small tree
      */
      
      found = found + findHash(toInsert, hash, vLength, position);              
    }

  rax_free(V);
  rax_free(bipartitions);
  
  return found;
}

/* method adapted for multifurcating trees, changes are: 
we now need information about the degree of an inner node, because it is not 2 anymore 
we also can have less than n-3 splits and therefore need a new parameter telling us the actual split number */
static void rec_extractBipartitionsMulti(unsigned int** bitvectors, int* seq, int arraysize, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction, int* taxonHasDegree, int maxSplits)
{
  int 
    i,
    j,
    o,
    taxon,
    numberOfSplits = 0,
    firstTaxon = taxonToReduction[first - 1] + 1,
    *V = (int *)rax_malloc((arraysize) * sizeof(int)),
    *bipartitions = (int*)rax_malloc((ntips - 3) * sizeof(int));

  hashNumberType 
    position;
  

  unsigned int 
    k,
    *toInsert = (unsigned int*)NULL,
    *bipartition = (unsigned int*)NULL;


  for(i = arraysize - 1; i >= 0; i--) 
    {
      V[i] = 1;
      /* instead of n-3 we stop after maxSplits */
      if(!isTip(seq[i],numsp) && (numberOfSplits < maxSplits))
	{ 
	  /* instead of j < 2, for multifurcating trees, they can have an abitrarily degree. 
	     We save this information in an array called taxonHasDegree and look it up quickly*/
	  for(j = 0; j < taxonHasDegree[seq[i] - 1] ; j++)      
	    V[i] = V[i] + V[i + V[i]];              
	  
	  int 
	    index = taxonToReduction[seq[i] - 1];

	  bipartitions[numberOfSplits] = index;

	  toInsert = bitvectors[index];

	  /* Extract Bipartition */
	  for(j = 1; j < V[i]; j++) 
	    {
	      if(isTip(seq[i + j],numsp))
		{
		  taxon = taxonToReduction[seq[i + j] - 1] + 1;
		  
		  toInsert[(taxon-1) / MASK_LENGTH] |= mask32[(taxon-1) % MASK_LENGTH];
		}
	      else
		{		  
		  int 
		    before = taxonToReduction[seq[i+j] - 1];
		  
		  bipartition = bitvectors[before];
		  
		  for(k = 0; k < vLength; k++) 
		    toInsert[k] |= bipartition[k];
		  
		  /* jump to the next subtree */
		  j = j + V[i+j] - 1;
		}        
	    }
	  
	  numberOfSplits += 1;          
	}
    }
  /* we now iterate to maxSplits instead of n-3 */
  for(i=0;i < maxSplits; i++) 
    {
      toInsert = bitvectors[bipartitions[i]];
      
      if(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]) 
	{          
	  /* Padding the last bits! */
	  if(ntips % MASK_LENGTH != 0) 
	    {            
	      for(o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--)         
		toInsert[vLength - 1] |= mask32[o-1];       
	    }

	  for(k=0;k < vLength; k++)         
	    toInsert[k] = ~toInsert[k];         
	}
        
      assert(!(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]));

      assert(vLength > 0);

      for(k=0; k < vLength; k++)    
	/* compute hash */
	position = oat_hash((unsigned char *)toInsert, sizeof(unsigned int) * vLength);
    
      position = position % hash->tableSize;
    
      /* re-hash to the new hash table that contains the bips of the large tree, pruned down 
	 to the taxa contained in the small tree
      */
      insertHashPlausibility(toInsert, hash, vLength, position);     
    }

  rax_free(V);
  rax_free(bipartitions);  
}




/*Preordertraversal of the big tree using bitVectorInitrav as reference and taking start->back node, 
number of tips and start->number as parameter and delivers a TaxonToPreOrderLabel and LabelToTaxon Array*/
static void preOrderTraversal(nodeptr p, int numsp, int start, int* array, int* backarray, int* pos)
{
  if(isTip(p->number, numsp))
    {
      array[p->number - 1] = *pos; 
      
      backarray[*pos] = p->number;
      
      *pos = *pos + 1;
      
      return;
    }
  else
    {
      nodeptr q = p->next;
      
      array[p->number - 1] = *pos; 
      
      backarray[*pos] = p->number;
      
      *pos = *pos + 1;
      
      /* get start element */
      if(p->back->number == start) 	
	preOrderTraversal(p->back, numsp, start, array, backarray, pos);       

      do
        {
          preOrderTraversal(q->back, numsp, start, array, backarray, pos);

          q = q->next;
        }
      while(q != p); 
    }
}



/*extract all smalltree taxa and store a list of all Taxon*/
static void rec_extractTaxa(int* smallTreeTaxa, int* taxonToReduction, nodeptr p, int numsp, int* pos, int* pos2)
{ 
  if(isTip(p->number, numsp))
    {
      smallTreeTaxa[*pos] = p->number; 
      taxonToReduction[p->number - 1] = *pos;
      *pos = *pos + 1;
      
      return;
    }
  else
    {    
      nodeptr 
	q = p->next;
      
      taxonToReduction[p->number - 1] = *pos2;      
      *pos2 = *pos2 + 1;
      
      do
	{
	  rec_extractTaxa(smallTreeTaxa, taxonToReduction, q->back, numsp, pos, pos2);
	  q = q->next;
	}
      while(q != p);            
  }
}

/* traverses the reference small tree and additionally extracting for every node its degree. It is stored in int* deg */
static void rec_preOrderTraversalMulti(nodeptr p, int numsp, int start, int* backarray, int* deg, int* pos)
{
  int degree = 0;

  if(isTip(p->number, numsp))
    { 
      backarray[*pos] = p->number;
      
      *pos = *pos + 1;
      
      deg[p->number - 1] = degree;
      return;
    }
  else
    {
      nodeptr 
	q = p->next;

      backarray[*pos] = p->number;

      *pos = *pos + 1;

      if(p->back->number == start) 
	{
	  rec_preOrderTraversalMulti(p->back, numsp, start, backarray, deg, pos);
	  degree += 1;
	} 

      do
        {
          rec_preOrderTraversalMulti(q->back, numsp, start, backarray, deg, pos);

          q = q->next;

          degree += 1;
         }
      while(q != p); 

      deg[p->number - 1] = degree;
    }
}



/* special function inits bitvectors to store bipartitions for dynamic programming*/
static unsigned int **rec_initBitVector(tree *tr, unsigned int vectorLength)
{
  unsigned int 
    **bitVectors = (unsigned int **)rax_malloc(sizeof(unsigned int*) * 2 * tr->ntips - 1);
  
  int 
    i;
  
  for(i = 0; i < (2 * tr->ntips - 1); i++)    
    bitVectors[i] = (unsigned int *)rax_calloc(vectorLength, sizeof(unsigned int));
        
  return bitVectors;
}


static void rec_freeBitVector(tree *tr, unsigned int **bitVectors)
{  
  int 
    i;
  
  for(i = 0; i < (2 * tr->ntips - 1); i++)    
    rax_free(bitVectors[i]);
}

/*euler traversal for binary and rooted trees*/                                                                                                                                                              
static void eulerTour(nodeptr p, int numsp, int* array, int* reference, int* pos, int* taxonToEulerIndex)
{
  array[*pos] = reference[p->number - 1];

  if (isTip(p->number, numsp)) 
    {
      if (taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;
  
  if(!isTip(p->number, numsp))
    {
      eulerTour(p->next->back, numsp, array, reference, pos, taxonToEulerIndex);
      
      array[*pos] = reference[p->number - 1]; 
      
      *pos = *pos + 1;
      
      eulerTour(p->next->next->back, numsp, array, reference, pos, taxonToEulerIndex);
      
      array[*pos] = reference[p->number - 1];
      
      *pos = *pos + 1;
    }
}

/*For unrooted Trees there is a special case for the arbitrary root which has degree 3 */
static void unrootedEulerTour(nodeptr p, int numsp, int* array, int* reference, int* pos, int* taxonToEulerIndex)
{
  array[*pos] = reference[p->number - 1];
  
  if (isTip(p->number, numsp)) 
    {
      if(taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;
  
  if(!isTip(p->number, numsp))    
    eulerTour(p->back, numsp,array,reference, pos, taxonToEulerIndex);   

  array[*pos] = reference[p->number - 1];

  if(isTip(p->number, numsp)) 
    {
      if(taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;

  if(!isTip(p->number, numsp))    
    eulerTour(p->next->back, numsp,array,reference, pos, taxonToEulerIndex);    

  array[*pos] = reference[p->number - 1];

  if(isTip(p->number, numsp)) 
    {
      if(taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;

  if(!isTip(p->number, numsp))    
    eulerTour(p->next->next->back, numsp,array,reference, pos, taxonToEulerIndex);

  array[*pos] = reference[p->number - 1];

  if(isTip(p->number, numsp)) 
    {
      if(taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;
}

//function for built-in quicksort

static int sortIntegers(const void *a, const void *b)
{
  int 
    ia = *(int *)(a),
    ib = *(int *)(b);

  if(ia == ib)
    return 0;

  if(ib > ia)
    return -1;
  else
    return 1;
}

/************************************************************************************/
/************************************* RF-OPT functions *****************************/

/* method adapted for multifurcating trees, changes are: 
we now need information about the degree of an inner node, because it is not 2 anymore 
we also can have less than n-3 splits and therefore need a new parameter telling us the actual split number */
static unsigned int** RFOPT_extractBipartitionsMulti(unsigned int** bitvectors, int* seq, int arraysize, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction, int* taxonHasDegree, int maxSplits)
{
  int 
    i,
    j,
    o,
    taxon,
    numberOfSplits = 0,
    firstTaxon = taxonToReduction[first - 1] + 1,
    *V = (int *)rax_malloc((arraysize) * sizeof(int)),
    *bipartitions = (int*)rax_malloc((ntips - 3) * sizeof(int));

  hashNumberType 
    position;
  

  unsigned int 
    k,
    *toInsert = (unsigned int*)NULL,
    *bipartition = (unsigned int*)NULL;

  //Store the bipartition bitvector into a seperate array 
  unsigned int
    **returnInserts = (unsigned int**)rax_malloc(maxSplits*sizeof(unsigned int*));


  for(i = arraysize - 1; i >= 0; i--) 
    {
      V[i] = 1;
      /* instead of n-3 we stop after maxSplits */
      if(!isTip(seq[i],numsp) && (numberOfSplits < maxSplits))
  { 
    /* instead of j < 2, for multifurcating trees, they can have an abitrarily degree. 
       We save this information in an array called taxonHasDegree and look it up quickly*/
    for(j = 0; j < taxonHasDegree[seq[i] - 1] ; j++)      
      V[i] = V[i] + V[i + V[i]];              
    
    int 
      index = taxonToReduction[seq[i] - 1];

    bipartitions[numberOfSplits] = index;

    toInsert = bitvectors[index];

    /* Extract Bipartition */
    for(j = 1; j < V[i]; j++) 
      {
        if(isTip(seq[i + j],numsp))
    {
      taxon = taxonToReduction[seq[i + j] - 1] + 1;
      
      toInsert[(taxon-1) / MASK_LENGTH] |= mask32[(taxon-1) % MASK_LENGTH];
    }
        else
    {     
      int 
        before = taxonToReduction[seq[i+j] - 1];
      
      bipartition = bitvectors[before];
      
      for(k = 0; k < vLength; k++) 
        toInsert[k] |= bipartition[k];
      
      /* jump to the next subtree */
      j = j + V[i+j] - 1;
    }        
      }
    
    numberOfSplits += 1;          
  }
    }
  /* we now iterate to maxSplits instead of n-3 */
  for(i=0;i < maxSplits; i++) 
    {
      toInsert = bitvectors[bipartitions[i]];
      
      if(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]) 
  {          
    /* Padding the last bits! */
    if(ntips % MASK_LENGTH != 0) 
      {            
        for(o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--)         
    toInsert[vLength - 1] |= mask32[o-1];       
      }

    for(k=0;k < vLength; k++)         
      toInsert[k] = ~toInsert[k];         
  }
        
      assert(!(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]));

      assert(vLength > 0);

      for(k=0; k < vLength; k++)    
  /* compute hash */
  position = oat_hash((unsigned char *)toInsert, sizeof(unsigned int) * vLength);
    
      position = position % hash->tableSize;
    
      /* re-hash to the new hash table that contains the bips of the large tree, pruned down 
   to the taxa contained in the small tree
      */
      insertHashPlausibility(toInsert, hash, vLength, position);

      returnInserts[i] = toInsert;     
    }

  rax_free(V);
  rax_free(bipartitions);

  return returnInserts;  
}


/* Additionally to the method above, rec_findAddBipartitions add the bipartitions into a second hashtable ind_hash and returns the bitVector of the induced tree */
static unsigned int** RFOPT_findAddBipartitions(unsigned int ** bitvectors, int* seq, int arraysize, int* translate, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, hashtable* ind_hash, int* taxonToReduction)
{
  int 
    i, 
    j,
    o, 
    taxon,
    found = 0,
    before;
    
  unsigned int 
    k, 
    *toInsert = (unsigned int*)NULL,
    *bipartition = (unsigned int*)NULL;

  unsigned int
    **returnInserts = (unsigned int**)rax_malloc((ntips-3)*sizeof(unsigned int*));
    
  hashNumberType 
    position;
  
  int 
    numberOfSplits = 0, /* stop after extracting n-3 splits */ 
    firstTaxon = taxonToReduction[first - 1] + 1,   
    *V = (int *)rax_malloc((arraysize) * sizeof(int)),
    *bipartitions = (int *)rax_malloc((ntips - 3) * sizeof(int));

  for(i = arraysize - 1; i >= 0; i--) 
    {
      V[i] = 1;
    
      /* Extract Bipartiton from inner node! */
      if(!isTip(translate[seq[i]],numsp) && (numberOfSplits < (ntips - 3)))
  {
    /* we can be sure of bifurcating trees, therefore j < deg(i) = 2 */
    for(j = 0; j < 2; j++)      
      V[i] = V[i] + V[i + V[i]];     

    /* find out the memory efficient index for this taxon which lies between 0.. (ntips - 1) instead of 1..mxtips  
       We have to subtract 1 because array starts from 0.. */  
    int 
      index = taxonToReduction[translate[seq[i]] - 1];

    /* save bipartition */
    bipartitions[numberOfSplits] = index;

    toInsert = bitvectors[index];

    //printf("We are working with %i \n", index);
    /* Calculate Bipartition */
    for(j = 1; j < V[i]; j++) 
      {   
        /* look if Preorderlabel at index i + j is a tip! */
        if(isTip(translate[seq[i + j]],numsp))
    {
      /* translate taxon number to a number between 0 and smalltreesize */      
      taxon = taxonToReduction[translate[seq[i + j]] - 1] + 1;

      /* set bit to one */
      toInsert[(taxon-1)  / MASK_LENGTH] |= mask32[(taxon-1) % MASK_LENGTH];
    }
        else
    {      
      before = taxonToReduction[translate[seq[i+j]] - 1];

      bipartition = bitvectors[before];

      for(k = 0; k < vLength; k++) 
        toInsert[k] |= bipartition[k];

      /* jump to the next relevant taxon */ 
      j = j + V[i+j] - 1;
    }        
      }
    /* count number of splits and stop at n-3 */
    numberOfSplits += 1;
  }
    }
  
  for(i=0;i < (ntips - 3); i++) 
    {
      toInsert = bitvectors[bipartitions[i]];

      /* if bitvector contains first taxon, use its complement */
      if(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]) 
  {                 
    /* Padding the last bits! */
    if(ntips % MASK_LENGTH != 0) 
      {            
        for(o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--)         
    toInsert[vLength - 1] |= mask32[o-1];       
      }
    
    for(k=0;k < vLength; k++)         
      toInsert[k] = ~toInsert[k];         
  }
      
      assert(!(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]));
      
      /* compute hash */
      position = oat_hash((unsigned char *)toInsert, sizeof(unsigned int) * vLength);
      
      position = position % hash->tableSize;
      
      /* re-hash to the new hash table that contains the bips of the large tree, pruned down 
   to the taxa contained in the small tree
      */
      insertHashPlausibility(toInsert, ind_hash, vLength, position); 

      returnInserts[i] = toInsert;    

      found = found + findHash(toInsert, hash, vLength, position);              
    }

  rax_free(V);
  rax_free(bipartitions);
  
  return returnInserts;
}


/* 
Extract the set of the bitvector and stores the taxa into an array called set which it returns
TODO only for one bitvector! Stop Element is -1 
*/
static int* extractSet(int* bitvector, int* smallTreeTaxa){
		int numberOfOnes = __builtin_popcount(bitvector[0]);
		//plus one because of the terminal number -1 to determine the end of the array
		int* set = rax_malloc((numberOfOnes + 1) * sizeof(int));
		int i = 0;
		int numberOfZerosBefore = 0;
		int extract = bitvector[0];
		while(i < numberOfOnes) {

			//Extract the first bit from extract and identify the related taxa number in SmallTree
			numberOfZerosBefore = __builtin_ctz(extract) + numberOfZerosBefore;
			//printf("Number of Zeros: %i \n", __builtin_ctz(extract));

      set[i] = smallTreeTaxa[numberOfZerosBefore];

			//Now move extract to the next significant bit
      numberOfZerosBefore = numberOfZerosBefore + 1;
			extract = extract >> numberOfZerosBefore; //Never moved!
			i++;
		}
		//Add a stop element in the array
		set[numberOfOnes] = -1;

    //Sort to get a comparable sequence for steps to follow
    qsort(set, numberOfOnes, sizeof(int), sortIntegers);

		return set;
}

/* 
Can extract two bipartitons and merge their sets
Edit compared to extractSet: for bitvector 2 we restart the loop at position set[i + numberOfOnesBip1]
*/
static int* extractSets(int* bitvector, int* bitvector2, int* smallTreeTaxa){
    int numberOfOnesBip1 = __builtin_popcount(bitvector[0]);
    int numberOfOnesBip2 = __builtin_popcount(bitvector2[0]);
    int numberOfOnes = numberOfOnesBip1 + numberOfOnesBip2;

    //plus one because of the terminal number -1 to determine the end of the array
    int* set = rax_malloc((numberOfOnes + 1) * sizeof(int));
    int i = 0;
    int numberOfZerosBefore = 0;
    int extract = bitvector[0];
    while(i < numberOfOnesBip1) {

      //Extract the first bit from extract and identify the related taxa number in SmallTree
      numberOfZerosBefore = __builtin_ctz(extract) + numberOfZerosBefore;
      //printf("Number of Zeros: %i \n", __builtin_ctz(extract));
      set[i] = smallTreeTaxa[numberOfZerosBefore];

      numberOfZerosBefore = numberOfZerosBefore + 1;
      //Now move extract to the next significant bit
      extract = extract >> numberOfZerosBefore;
      i++;
    }


    //restart the whole process for bitvector 2
    i = 0;
    numberOfZerosBefore = 0;
    extract = bitvector2[0];

    while(i < numberOfOnesBip2) {

      //Extract the first bit from extract and identify the related taxa number in SmallTree
      numberOfZerosBefore = __builtin_ctz(extract) + numberOfZerosBefore;
      //printf("Number of Zeros: %i \n", __builtin_ctz(extract));
      set[i + numberOfOnesBip1] = smallTreeTaxa[numberOfZerosBefore];


      numberOfZerosBefore = numberOfZerosBefore + 1;
      //Now move extract to the next significant bit
      extract = extract >> numberOfZerosBefore;
      i++;
    }
    //Add a stop element in the array
    set[numberOfOnes] = -1;

    //Sort to get a comparable sequence for steps to follow
    qsort(set, numberOfOnes, sizeof(int), sortIntegers);
 
    return set;
}


/************************************* helper functions ***************************/

/* converts integer into binary representation */
static char *int2bin(int a, char *buffer, int buf_size) {
	    buffer += (buf_size - 1);

		    for (int i = 31; i >= 0; i--) {
				        *buffer-- = (a & 1) + '0';

						        a >>= 1;
								    }

			    return buffer;
}

/* function for built-in quicksort sorting after number of bits */
static int sortBipartitions(const void *a, const void *b) 
{
  int 
   ia = *(int *)(a),
   ib = *(int *)(b),
   bits_ia = __builtin_popcount(ia),
   bits_ib = __builtin_popcount(ib);
  
   if(bits_ia == bits_ib)
     return 0;
   if(bits_ib > bits_ia)
     return -1;
   else
     return 1;
}

/* sort multidimensional arrays with different size */
static int sortSets(const void *a, const void *b)
{
  int*
    ia = *(int**)(a);
  int*
    ib = *(int**)(b);
  int
    j = 0;

  //while it is equal, we look for the first comparison which is not equal
  while((ia[j] != -1) && (ib[j] != -1)){
  
  int
    ix = ia[j],
    iy = ib[j];

  //printf("iy: %i ix: %i \n",ix,iy);

  // if(ix == iy)
  //    return 0;

  if(ix > iy)
     return -1;
  if(ix < iy)
     return 1;

  //counter increment
  j++;
  
  }

  if((ia[j] == -1) || (ib[j] == -1)){
    if(ia[j] > ib[j])
     return -1;
    if(ia[j] < ib[j])
     return 1;
  }
  return 0;    
}

/* Checks if two arrays are identical and returns 1 and 0 */
static int isSameDropSet(int* a, int* b) {
  int j = 0;

  //while it is equal, we look for the first comparison which is not equal
  while((a[j] != -1) && (b[j] != -1)){
  int
    x = a[j],
    y = b[j];

  //if not equal, return with false
  if(x != y)
    return 0;

  j++;

  }

  if((a[j] == -1) || (b[j] == -1)){
    if(a[j] != b[j])
     return 0;
  }
  return 1; // made it through all tests, it is the same    
}

/* Checks if check already is inside sets between 0 ... numberOfSets 
    returns 0 if its not containing
    returns index+1 if it contains the element 
*/
static int contains(int* check, int** sets, int numberOfSets) {

  for (int i = 0; i < numberOfSets; i++) {
    
    int* dropset = sets[i];
    
    if(isSameDropSet(sets[i],check)){
      return (i+1);
    }
  }

  return 0;

}

//Get the index for which array arr is max
static int getMax(int* arr, int size) {

  int max = 0;

  for (int i = 1; i < size; i++) {

    if(arr[i] > arr[max]) {
      max = i;
    }

  }

  return max;
}

/******************************* Bit Manipulation functionality ***********************************/


static int setBit(int bitVector, int pos) {
	
	bitVector |= 1 << pos;

	return bitVector;
}

static int clearBit(int bitVector, int pos) {
	
	bitVector &= ~(1 << pos);

	return bitVector;
}

static int checkBit(int bitVector, int pos) {
	
	int bit = 0;
	
	bit = bitVector & (1 << pos);

	return bit;
}

//Method for setting bits in bitvectors
static int* setBitBV(int* bitVector, int pos, int vLength) {

  int vectorIndex = pos / MASK_LENGTH;
  int relativePos = pos % MASK_LENGTH;

  if(vectorIndex > vLength) {
    printf("setBITBV: Error! Position invalid \n");
    return 0;
  
  } else {
    bitVector[vectorIndex] = setBit(bitVector[vectorIndex], relativePos);
    return bitVector;
  } 
}

//Use to setup a mask to clear the offset bits
static int setOffSet(int mask, int offset) {
  
  for(int i = 0; i < offset; i++) {
    mask = setBit(mask, i);
  }

  return mask; 
}

//Use to setup a mask of existing bipartitions
static int getExistingBips(int mask, int offset, int bvecs_deletedBips) {

  mask = setOffSet(mask, offset);

  mask = mask & ~bvecs_deletedBips;

  return mask;
}

static void printBitVector(int bitVector) {
  char buffer[33];
  buffer[32] = '\0';
  int2bin(bitVector, buffer, 32);
  printf("\n==> BitVector = %s \n", buffer);
}

static void printSet(int* set) {
  int i = 0;
  printf("Set: ");
  while(set[i] != -1) {
    printf("%i ", set[i]);
    i++;
  }
  printf("\n");
}

//Takes as input a bitvector and returns a new bitvector OLD
static int getBipsOfDropSet(int bvec_bips, int dropsetNumber, int* numberOfBipsPerSet, int** bipsOfDropSet) {
    //Now iterate through bipsOfDropSet list
    for(int l = 0; l < numberOfBipsPerSet[dropsetNumber]; l++) {

      //Get the index list of bips for this Drop Set
      int* bips = bipsOfDropSet[dropsetNumber];
      int b_index = bips[l];

      //Generate a bitvector 
      bvec_bips = setBit(bvec_bips, b_index); 
    }

    return bvec_bips;
}

/**********************************************************************************/

static int* extractSetFromBitVector(unsigned int* bitvector, int* smallTreeTaxa, unsigned int vLength){

    int numberOfOnes = 0;

    for(int i = 0; i < vLength; i++) {
      numberOfOnes = numberOfOnes + __builtin_popcount(bitvector[i]);
    }

    //plus one because of the terminal number -1 to determine the end of the array
    int* set = rax_malloc((numberOfOnes + 1) * sizeof(int));

    int i = 0;

    int numberOfZerosBefore = 0;

    int startVector = 0;

    unsigned int extract = bitvector[startVector];

    while(i < numberOfOnes) {

      if(extract != 0) {

        //Extract the first bit from extract and identify the related taxa number in SmallTree
        numberOfZerosBefore = __builtin_ctz(extract) + numberOfZerosBefore;
        //printf("Number of Zeros: %i \n", __builtin_ctz(extract));

        //Consider all bitvectors before this vector. I.e. MASK_Length = 32 , 
        //this is the third index of the second bitvector, we have 3 + 2 * 32 as real taxon index
        int pastVectors = startVector * MASK_LENGTH;

        set[i] = smallTreeTaxa[numberOfZerosBefore + pastVectors];

        //Now move extract to the next significant bit
        numberOfZerosBefore = numberOfZerosBefore + 1;
        extract = extract >> numberOfZerosBefore;
        i++;

      } else {

        //if bitvector has no 1s anymore, go to the next one
        startVector++;
        extract = bitvector[startVector];

        assert(startVector <= vLength);

      }

    }

    //Add a stop element in the array
    set[numberOfOnes] = -1;

    //Sort to get a comparable sequence for steps to follow
    qsort(set, numberOfOnes, sizeof(int), sortIntegers);

    return set;
}

//Merge two dropsets into one set (ending element is -1)
static int* mergeSets(int* set1, int* set2) {

  int i = 0;
  int j = 0;
  while(set1[i] != -1) {
    i++;
  }
  while(set2[j] != -1) {
    j++;
  }

  int sizeOfSet = i + j;
  int* set = rax_malloc((sizeOfSet + 1) * sizeof(int));

  //Now copy both sets into set
  i = 0;
  j = 0;
  while(set1[i] != -1) {
    set[i] = set1[i];
    i++;
  }

  while(set2[j] != -1) {
    set[j + i] = set2[j];
    j++;
  }

  set[j+i] = -1;

  qsort(set, sizeOfSet, sizeof(int), sortIntegers);

  return set;

}

static int* extractSetsFromBitVector(unsigned int* bitvector, unsigned int* bitvector2, int* smallTreeTaxa, unsigned int vLength){

    int* set1; 
    int* set2; 
    int* set;

    set1 = extractSetFromBitVector(bitvector, smallTreeTaxa, vLength);
    set2 = extractSetFromBitVector(bitvector2, smallTreeTaxa, vLength);

    set = mergeSets(set1, set2);

    return set;
}

//Calculate the DropSet given two BitVectors of the same length
static int* getDropSetFromBitVectors(unsigned int* indBip, unsigned int* sBip, unsigned int vLength, int treenumber, int* taxaPerTree, int* smallTreeTaxa){

  //Calculate the dropsets by comparing two bipartitions
  //Determine the smallest of two sets (xANDx*,yANDy*) OR (xANDy*,yANDx*)
  unsigned int* x1_x2 = (unsigned int*)rax_malloc(sizeof(int) * vLength);
  unsigned int* y1_y2 = (unsigned int*)rax_malloc(sizeof(int) * vLength);

  unsigned int* x1_y2 = (unsigned int*)rax_malloc(sizeof(int) * vLength);
  unsigned int* y1_x2 = (unsigned int*)rax_malloc(sizeof(int) * vLength);

  unsigned int* selectedSet1 = (unsigned int*)NULL;
  unsigned int* selectedSet2 = (unsigned int*)NULL;
  unsigned int selectedCount1;
  unsigned int selectedCount2;

  unsigned int count1 = 0;
  unsigned int count2 = 0;
  unsigned int count3 = 0;
  unsigned int count4 = 0;

  //Calculate a logical operation on BitVectors
  for(int i = 0; i < vLength; i++) {

    //Calculating the two different Bipartition choices
    x1_x2[i] = indBip[i] & sBip[i];
    y1_y2[i] = ~(indBip[i]) & ~(sBip[i]);

    x1_y2[i] = indBip[i] & ~(sBip[i]);
    y1_x2[i] = ~(indBip[i]) & sBip[i];
  }

  //Toggle offset if required
  int ntips = taxaPerTree[treenumber];

  if(ntips % MASK_LENGTH != 0) 
      {            
        for(int o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--) {        
          x1_x2[vLength - 1] &= ~mask32[o-1];  
          y1_y2[vLength - 1] &= ~mask32[o-1];
          x1_y2[vLength - 1] &= ~mask32[o-1];
          y1_x2[vLength - 1] &= ~mask32[o-1];     
        }
      }


  for(int i = 0; i < vLength; i++) {
    //Calculate number of bits of the resulting set calculations  
    count1 = count1 + __builtin_popcount(x1_x2[i]);
    count2 = count2 + __builtin_popcount(y1_y2[i]);

    count3 = count3 + __builtin_popcount(x1_y2[i]); 
    count4 = count4 + __builtin_popcount(y1_x2[i]);
  }

  //Now choose which set we use for further processing
  if((count1 + count2)<(count3 + count4)) {

    selectedSet1 = x1_x2;
    selectedSet2 = y1_y2; 

    selectedCount1 = count1;
    selectedCount2 = count2;   

    free(x1_y2);
    free(y1_x2);

  } else {

    selectedSet1 = x1_y2;
    selectedSet2 = y1_x2;  

    selectedCount1 = count3;
    selectedCount2 = count4;  
    
    free(x1_x2);
    free(y1_y2);
  }

  //-1 is the stop element
  int* set;

  //Now we check how we extract the DropSets
  if(selectedCount1 == 0 & selectedCount2 == 0) {
    //This is a matching bipartition!
    int res[2] = {0,-1};
    set = res;

  } else {

    //We have three different cases
    if(selectedCount1 != 0 & selectedCount2 == 0) {
      //Only first dropset is non empty, therefore only look at the first selected Set

      set = extractSetFromBitVector(selectedSet1, smallTreeTaxa, vLength);
      
    } else if(selectedCount1 == 0 & selectedCount2 != 0) {
      //Only second dropset is non empty, therefore only look at the second selected Set

      set = extractSetFromBitVector(selectedSet2, smallTreeTaxa, vLength);

    } else {
      //Both dropsets are non empty

      set = extractSetsFromBitVector(selectedSet1, selectedSet2, smallTreeTaxa, vLength);

    }
  }

  return set;
}

//Filter for unique dropsets and return number of unique sets
static int getUniqueDropSets(int** sets, int** uniqSets, int* setsToUniqSets, int numberOfSets) {
  
  int numberOfUniqueSets = 0;

  int dropSetCount = 0;

  //Now generate Graph and calculate the Scores for each bipartitions
  for(int k = 0; k < numberOfSets; k++){
    
    //Get the dropset at position k
    int* dropset = sets[k];

    //Check if set is already in uniqSets
    int containIndex = contains(dropset,uniqSets,numberOfUniqueSets);
        
    if(!containIndex) {

      uniqSets[numberOfUniqueSets] = dropset;

      //Add the new index to the translation array
      setsToUniqSets[k] = numberOfUniqueSets; 

      numberOfUniqueSets++;

    } else {
      
      //If it is already inside uniqSet
      setsToUniqSets[k] = (containIndex - 1);
    }
  }

  return numberOfUniqueSets;
}

//Calculates all dropsets of two given bipartition lists
static void calculateDropSets(unsigned int*** indBipsPerTree, unsigned int*** sBipsPerTree, int** sets, int** smallTreeTaxaList, int* bipsPerTree, 
  int* taxaPerTree, unsigned int* vectorLengthPerTree, int numberOfTrees) {

  int countSets = 0;

  //First iterate through all trees
  for(int i = 0; i< numberOfTrees; i++) {

    //Get all induced Bips of this tree
    unsigned int **indBips = indBipsPerTree[i];

    //Get all small Bips of this tree
    unsigned int **sBips = sBipsPerTree[i];

    //Now go through all Bips of this tree
    for(int j = 0; j < bipsPerTree[i]; j++) {

      //Get the bitVector of this Bips
      unsigned int *indBip = indBips[j];

      //Now iterate through all small Bips in the tree
      for(int k = 0; k < bipsPerTree[i]; k++) {
        
        //get small Bip
        unsigned int *sBip = sBips[k];

        //extract DropSets from the comparision
        sets[countSets] = getDropSetFromBitVectors(indBip, sBip, vectorLengthPerTree[i], i, taxaPerTree, smallTreeTaxaList[i]);

        countSets++;

      }
    }
  }
}

static void detectInitialMatchings(int** sets, int* matchingVector, int* bipsPerTree, int numberOfTrees,  int vLength) { 

  //We use these variables to iterate through all sets and bips
  int countBips = 0;
  int dropSetCount = 0;
  //Now generate Graph and calculate the Scores for each bipartitions

  //First iterate through all trees 
  for(int i = 0; i < numberOfTrees; i++ ) {

    //Iterate through all bipartitions of each trees
    for(int j = 0; j < bipsPerTree[i]; j++) {

      //Compare each bip with all bips
      int dropSetsPerBip = bipsPerTree[i]; 
      
      //Will be set to one if dropset is matching, resulting in a score = 0
      int matching = 0; 
      
      //Check all dropsets of each bipartition
      for(int k = 0; k < dropSetsPerBip; k++){
        
        int* dropset = sets[dropSetCount + k];

        //Test if matching
        if(dropset[0] == 0){
          matching = 1;
        }
      }

      if (matching) {
        //don't do a thing 
        //bvec_scores[countBips] = 0; 
      } else {
        //Set BitVector on position of the bip
        matchingVector = setBitBV(matchingVector, countBips, vLength);
      }
      
      //Index of starting DropSets in Sets for the next Bip 
      dropSetCount = dropSetCount + dropSetsPerBip; 

      countBips++;
    }

  }

  printf("dropSetCount : %i \n", dropSetCount);
  printf("bipsCount : %i \n", countBips);
}


/**********************************************************************************/
/*************************** RF-OPT Algorihm **************************************/
/**********************************************************************************/

#define _USE_RF_OPT

#ifdef _USE_RF_OPT

//Use the plausibility checker overhead
void plausibilityChecker(tree *tr, analdef *adef)
{
  FILE
	*treeFile, 
    *treeFile2,
    *rfFile;
  
  tree 
    *smallTree = (tree *)rax_malloc(sizeof(tree));

  char 
    rfFileName[1024];

  int
    numberOfTreesAnalyzed = 0,
    i;

  double 
    avgRF = 0.0,
    sumEffectivetime = 0.0;

  /* set up an output file name */

  strcpy(rfFileName,         workdir);  
  strcat(rfFileName,         "RAxML_RF-Distances.");
  strcat(rfFileName,         run_id);

  rfFile = myfopen(rfFileName, "wb");  

  assert(adef->mode ==  PLAUSIBILITY_CHECKER);

  /* open the big reference tree file and parse it */

  treeFile = myfopen(tree_file, "r");

  printBothOpen("Parsing reference tree %s\n", tree_file);

  treeReadLen(treeFile, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);


  assert(tr->mxtips == tr->ntips);
  
  /*************************************************************************************/
  /* Preprocessing Step */

  double 
    preprocesstime = gettime();
  
  /* taxonToLabel[2*tr->mxtips - 2]; 
  Array storing all 2n-2 labels from the preordertraversal: (Taxonnumber - 1) -> (Preorderlabel) */
  int 
    *taxonToLabel  = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int)),

    /* taxonHasDeg[2*tr->mxtips - 2] 
    Array used to store the degree of every taxon, is needed to extract Bipartitions from multifurcating trees 
    (Taxonnumber - 1) -> (degree of node(Taxonnumber)) */

    *taxonHasDeg = (int *)rax_calloc((2*tr->mxtips - 2),sizeof(int)),

    /* taxonToReduction[2*tr->mxtips - 2]; 
  Array used for reducing bitvector and speeding up extraction: 

  (Taxonnumber - 1) -> Index in smallTreeTaxa (starting from 0)
  which is also:
  (Taxonnumber - 1) -> (0..1 (increment count of taxa appearing in small tree))
  (Taxonnumber - 1) -> (0..1 (increment count of inner nodes appearing in small tree)) */

    *taxonToReduction = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int));
    
  int 
    newcount = 0; //counter used for correct traversals

  /* labelToTaxon[2*tr->mxtips - 2];
  is used to translate between Perorderlabel and p->number: (Preorderlabel) -> (Taxonnumber) */
  int 
    *labelToTaxon = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int));
  
  /* Preorder-Traversal of the large tree */
  preOrderTraversal(tr->start->back,tr->mxtips, tr->start->number, taxonToLabel, labelToTaxon, &newcount);

  newcount = 0; //counter set to 0 to be now used for Eulertraversal

  /* eulerIndexToLabel[4*tr->mxtips - 5]; 
  Array storing all 4n-5 PreOrderlabels created during eulertour: (Eulerindex) -> (Preorderlabel) */
  int* 
    eulerIndexToLabel = (int *)rax_malloc((4*tr->mxtips - 5) * sizeof(int));

  /* taxonToEulerIndex[tr->mxtips]; 
  Stores all indices of the first appearance of a taxa in the eulerTour: (Taxonnumber - 1) -> (Index of the Eulertour where Taxonnumber first appears) 
  is used for efficient computation of the Lowest Common Ancestor during Reconstruction Step
  */
  int*
    taxonToEulerIndex  = (int *)rax_malloc((tr->mxtips) * sizeof(int));

  /* Init taxonToEulerIndex and taxonToReduction */
  int 
    ix;

  for(ix = 0; ix < tr->mxtips; ++ix)    
    taxonToEulerIndex[ix] = -1;    
  
  for(ix = 0; ix < (2*tr->mxtips - 2); ++ix)    
    taxonToReduction[ix] = -1;    


  /* Eulertraversal of the large tree*/
  unrootedEulerTour(tr->start->back,tr->mxtips, eulerIndexToLabel, taxonToLabel, &newcount, taxonToEulerIndex);

  /* Creating RMQ Datastructure for efficient retrieval of LCAs, using Johannes Fischers Library rewritten in C
  Following Files: rmq.h,rmqs.c,rmqs.h are included in Makefile.RMQ.gcc */
  RMQ_succinct(eulerIndexToLabel,4*tr->mxtips - 5);

  double 
    preprocessendtime = gettime() - preprocesstime;

  /* Proprocessing Step End */
  /*************************************************************************************/

  printBothOpen("The reference tree has %d tips\n", tr->ntips);

  fclose(treeFile);
  
  /***********************************************************************************/
  /* RF-OPT Preprocessing Step */
  /***********************************************************************************/

  /* now see how many small trees we have */
  treeFile = getNumberOfTrees(tr, bootStrapFile, adef);
  treeFile2 = getNumberOfTrees(tr, bootStrapFile, adef);

  checkTreeNumber(tr->numberOfTrees, bootStrapFile);

  /* allocate a data structure for parsing the potentially mult-furcating tree */

  allocateMultifurcations(tr, smallTree);

  /* Start Additional preprocessing step */

  int 
	  numberOfBips = 0,
	  numberOfSets = 0;

  //Stores the number of bips of each tree
  int *bipsPerTree = (int *)rax_malloc(tr->numberOfTrees * sizeof(int));

  //Stores the number of taxa for each tree
  int *taxaPerTree = (int *)rax_malloc(tr->numberOfTrees * sizeof(int));

  //To calculate all bipartitions, I created a new treeFile2 and a new getNumberOfTrees method!!
  for(i = 0; i < tr->numberOfTrees; i++) {

	  int this_treeBips = readMultifurcatingTree(treeFile2, smallTree, adef, TRUE);
	
	  numberOfBips = numberOfBips + this_treeBips;
	
	  numberOfSets = numberOfSets + this_treeBips * this_treeBips;

    bipsPerTree[i] = this_treeBips;
  }

  printf("numberOfBips: %i , numberOfSets: %i \n \n", numberOfBips, numberOfSets);	

  //stores induced bips (OLD?)
  unsigned int *ind_bips = (unsigned int *)rax_malloc(numberOfBips * sizeof(unsigned int));

  //stores smalltree bips (OLD?)
  unsigned int *s_bips = (unsigned int *)rax_malloc(numberOfBips * sizeof(unsigned int));

  //stores small bips per tree
  unsigned int ***sBipsPerTree = (unsigned int ***)rax_malloc(tr->numberOfTrees * sizeof(unsigned int**));

  //stores induced bips per tree
  unsigned int ***indBipsPerTree = (unsigned int ***)rax_malloc(tr->numberOfTrees * sizeof(unsigned int**));

  //stores vLength of each tree for processing bitVectors
  unsigned int *vectorLengthPerTree = (unsigned int *)rax_malloc(tr->numberOfTrees * sizeof(unsigned int*));

  //stores the corresponding tree number for each bip
  int *treenumberOfBip = (int *)rax_malloc(numberOfBips * sizeof(int));

  //Stores all dropsets of all trees 
  int **sets = (int **)rax_malloc(numberOfSets * sizeof(int*)); 
  
  //For each tree, stores a translation array from taxanumber smalltree->largetree
  int **smallTreeTaxaList = (int **)rax_malloc(tr->numberOfTrees * sizeof(int*)); 

  //For each tree, store a translation array from taxanumber largetree->smalltree
  int **taxonToReductionList = (int **)rax_malloc(tr->numberOfTrees * sizeof(int*));

  //I use these variables as global variables for all trees to determine the max number of possible sets to generate a static array
  int currentBips = 0;
  int currentSmallBips = 0;
  int currentSets = 0;

  //int currentTree = 0; already there in number of trees analyzed
	
  //Prefill sets with -1s
  for(int it = 0;it < (numberOfSets);it++){
	int fill[1] = {-1};
	sets[it] = fill; 
  }
  
  /***********************************************************************************/
  /* RF-OPT Preprocessing Step End */
  /***********************************************************************************/

  /* loop over all small trees */

  for(i = 0; i < tr->numberOfTrees;  i++)
    {      
      int
		numberOfSplits = readMultifurcatingTree(treeFile, smallTree, adef, TRUE);
      
      if(numberOfSplits > 0)
	{
	  int
	    firstTaxon;           

	  double
	    rec_rf,
	    maxRF;

	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Small tree %d has %d tips and %d bipartitions\n", i, smallTree->ntips, numberOfSplits);    
	  
	  /* compute the maximum RF distance for computing the relative RF distance later-on */
	  
	  /* note that here we need to pay attention, since the RF distance is not normalized 
	     by 2 * (n-3) but we need to account for the fact that the multifurcating small tree 
	     will potentially contain less bipartitions. 
	     Hence the normalization factor is obtained as n-3 + numberOfSplits, where n-3 is the number 
	     of bipartitions of the pruned down large reference tree for which we know that it is 
	     bifurcating/strictly binary */
	  
	  maxRF = (double)(2 * numberOfSplits);
	  
	  /* now get the index of the first taxon of the small tree.
	     we will use this to unambiguously store the bipartitions 
	  */
	  
	  firstTaxon = smallTree->start->number;

    //Saves the number of taxa in the tree (for RF-OPT)
    taxaPerTree[numberOfTreesAnalyzed] = smallTree->ntips; 
	  
	  /***********************************************************************************/
	  /* Reconstruction Step */
	  
	  double 
	    time_start = gettime();
	  
	  /* Init hashtable to store Bipartitions of the induced subtree T|t_i */
	  /* 
	     using smallTree->ntips instead of smallTree->mxtips yields faster code 
	     e.g. 120 versus 128 seconds for 20,000 small trees on my laptop 
	   */
	  hashtable
	    *s_hash = initHashTable(smallTree->ntips * 4);


    /* Init hashtable to store Bipartitions of the reference tree t_i*/
    hashtable
      *ind_hash = initHashTable(smallTree->ntips * 4);
	  
	  /* smallTreeTaxa[smallTree->ntips]; 
	     Stores all taxa numbers from smallTree into an array called smallTreeTaxa: (Index) -> (Taxonnumber)  */
	  int* 
	    smallTreeTaxa = (int *)rax_malloc((smallTree->ntips) * sizeof(int));
	  
	  /* counter is set to 0 for correctly extracting taxa of the small tree */
	  newcount = 0; 
	  
	  int 
	    newcount2 = 0;
	  
	  /* seq2[2*smallTree->ntips - 2]; 
	     stores PreorderSequence of the reference smalltree: (Preorderindex) -> (Taxonnumber) */
	  int* 
	    seq2 = (int *)rax_malloc((2*smallTree->ntips - 2) * sizeof(int));
	  
    /* used to store the vectorLength of the bitvector */
	  unsigned int 
	    vectorLength;
	  
	  /* extract all taxa of the smalltree and store it into an array, 
	     also store all counts of taxa and nontaxa in taxonToReduction */
	  rec_extractTaxa(smallTreeTaxa, taxonToReduction, smallTree->start, smallTree->mxtips, &newcount, &newcount2);
	  
	  rec_extractTaxa(smallTreeTaxa, taxonToReduction, smallTree->start->back, smallTree->mxtips, &newcount, &newcount2);
	  
	  /* counter is set to 0 to correctly preorder traverse the small tree */
	  newcount = 0;
	  
	  /* Preordertraversal of the small reference tree and save its sequence into seq2 for later extracting the bipartitions, it
	     also stores information about the degree of every node */
	  
	  rec_preOrderTraversalMulti(smallTree->start->back,smallTree->mxtips, smallTree->start->number, seq2, taxonHasDeg, &newcount);
	  
	  /* calculate the bitvector length */
	  if(smallTree->ntips % MASK_LENGTH == 0)
	    vectorLength = smallTree->ntips / MASK_LENGTH;
	  else
	    vectorLength = 1 + (smallTree->ntips / MASK_LENGTH); 


    /***********************************************************************************/
    /* RF-OPT Additional Preprocessing storing Bipartitions */
    /***********************************************************************************/    

    vectorLengthPerTree[numberOfTreesAnalyzed] = vectorLength;
	  
	  unsigned int 
	    **bitVectors = rec_initBitVector(smallTree, vectorLength);

    unsigned int
      **sBips;

	  /* store all non trivial bitvectors using an subtree approach for the reference subtree and 
	     store it into a hashtable, this method was changed for multifurcation */
	  sBips = RFOPT_extractBipartitionsMulti(bitVectors, seq2, newcount,tr->mxtips, vectorLength, smallTree->ntips, 
				       firstTaxon, s_hash, taxonToReduction, taxonHasDeg, numberOfSplits);

    sBipsPerTree[numberOfTreesAnalyzed] = sBips;

    /***********************************************************************************/
    /* End RF-OPT Additional Preprocessing storing Bipartitions */
    /***********************************************************************************/  
	  
	  /* counter is set to 0 to be used for correctly storing all EulerIndices */
	  newcount = 0; 
	  
	  /* smallTreeTaxonToEulerIndex[smallTree->ntips]; 
	     Saves all first Euler indices for all Taxons appearing in small Tree: 
	     (Index) -> (Index of the Eulertour where the taxonnumber of the small tree first appears) */
	  int* 
	    smallTreeTaxonToEulerIndex = (int *)rax_malloc((smallTree->ntips) * sizeof(int));
	  
	  /* seq[(smallTree->ntips*2) - 1] 
	     Stores the Preordersequence of the induced small tree */
	  int* 
	    seq = (int *)rax_malloc((2*smallTree->ntips - 1) * sizeof(int));
	  
	  
	  /* iterate through all small tree taxa */
	  for(ix = 0; ix < smallTree->ntips; ix++) 
	    {        
	      int 
	        taxanumber = smallTreeTaxa[ix];
	      
	      /* To create smallTreeTaxonToEulerIndex we filter taxonToEulerIndex for taxa in the small tree*/
	      smallTreeTaxonToEulerIndex[newcount] = taxonToEulerIndex[taxanumber-1]; 
	      
	      /* Saves all Preorderlabel of the smalltree taxa in seq*/
	      seq[newcount] = taxonToLabel[taxanumber-1];
	      
	      newcount++;
	    }
	  
	  /* sort the euler indices to correctly calculate LCA */
	  //quicksort(smallTreeTaxonToEulerIndex,0,newcount - 1);             
	  
	  qsort(smallTreeTaxonToEulerIndex, newcount, sizeof(int), sortIntegers);
	  
	  //printf("newcount2 %i \n", newcount2);      
	  /* Iterate through all small tree taxa */
	  for(ix = 1; ix < newcount; ix++)
	    {  
	      /* query LCAs using RMQ Datastructure */
	      seq[newcount - 1 + ix] =  eulerIndexToLabel[query(smallTreeTaxonToEulerIndex[ix - 1],smallTreeTaxonToEulerIndex[ix])]; 	 
	      
	      /* Used for dynamic programming. We save an index for every inner node:
		 For example the reference tree has 3 inner nodes which we saves them as 0,1,2.
		 Now we calculate for example 5 LCA to construct the induced subtree, which are also inner nodes. 
		 Therefore we mark them as 3,4,5,6,7  */
	      
	      taxonToReduction[labelToTaxon[seq[newcount - 1 + ix]] - 1] = newcount2;
	      
	      newcount2 += 1;
	    }
	  
	  /* sort to construct the Preordersequence of the induced subtree */
	  //quicksort(seq,0,(2*smallTree->ntips - 2));
	  
	  qsort(seq, (2 * smallTree->ntips - 2) + 1, sizeof(int), sortIntegers);
	  
	  /* calculates all bipartitions of the reference small tree and count how many bipartition it 
    shares with the induced small tree and stores those bipartitions in a additional hashtable called ind_hash */
	  int 
	    rec_bips = 0;

    unsigned int
      **indBips;

    indBips = RFOPT_findAddBipartitions(bitVectors, seq,(2*smallTree->ntips - 1), labelToTaxon, tr->mxtips, vectorLength, smallTree->ntips, firstTaxon, s_hash, ind_hash, taxonToReduction);
      
    indBipsPerTree[numberOfTreesAnalyzed] = indBips; 

    /* calculates all bipartitions of the reference small tree and put them into ind_hash*/
    // rec_extractBipartitionsMulti(bitVectors, seq2, (2*smallTree->ntips - 1),tr->mxtips, vectorLength, smallTree->ntips, 
             // firstTaxon, s_hash, taxonToReduction, taxonHasDeg, numberOfSplits);


	  /* Reconstruction Step End */
	  /***********************************************************************************/
	  
	  double 
	    effectivetime = gettime() - time_start;
	  
	  /*
	    if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Reconstruction time: %.10f secs\n\n", effectivetime);
	  */
	  
	  /* compute the relative RF */


    /***********************************************************************************/
	  /* RF-OPT Save Translation Vectors */
    /***********************************************************************************/
	    
    //copy array taxonToReduction because it is originally defined in preprocessing step
    int * taxonToReductionCopy = (int *)rax_malloc((tr->mxtips)*sizeof(int));

    memcpy(taxonToReductionCopy,taxonToReduction,(tr->mxtips)*sizeof(int));

	  //storing smallTree and taxonToReduction Arrays for further usage
	  smallTreeTaxaList[numberOfTreesAnalyzed] = smallTreeTaxa;

    taxonToReductionList[numberOfTreesAnalyzed] = taxonToReductionCopy;	  

    int this_currentSmallBips = 0; //Variable resets everytime for each tree analyzed
    
    
    /***********************************************************************************/
    /* End RF-OPT Save Translation Vectors */
    /***********************************************************************************/
  

	  rec_rf = (double)(2 * (numberOfSplits - rec_bips)) / maxRF;
	  
	  assert(numberOfSplits >= rec_bips);	  	 

	  avgRF += rec_rf;
	  sumEffectivetime += effectivetime;
	  
	  //if(numberOfTreesAnalyzed % 100 == 0)
	  printBothOpen("Relative RF tree %d: %f\n\n", i, rec_rf);
	  
	  fprintf(rfFile, "%d %f\n", i, rec_rf);
	  
	  //rax_free(smallTreeTaxa); //Need it for calculating the SmallTreeTaxaList after all iterations!
	  rax_free(seq);
	  rax_free(seq2);
	  rax_free(smallTreeTaxonToEulerIndex);

	  numberOfTreesAnalyzed++; //Counting the number of trees analyzed
	  }

  }// End of Small Tree Iterations
  

  /***********************************************************************************/
  /* RF-OPT DropSet Calculation using BitVectors */
  /***********************************************************************************/

  printf("===> BitVector Set Calculation \n");

  //Calculate dropsets of two given bips lists and extract all sets into array sets. Each set has following format
  //dropset = {taxa_1,taxa_2,...,taxa_n,-1};
  calculateDropSets(indBipsPerTree, sBipsPerTree, sets, smallTreeTaxaList, bipsPerTree, 
  taxaPerTree, vectorLengthPerTree, tr->numberOfTrees);

  /***********************************************************************************/
  /* RF-OPT Graph Construction */
  /***********************************************************************************/


  /*
    Filter for unique sets
  */
  printf("===> Filter for unique Sets (naive)...\n");
  /* unique sets array data structures */
  int** uniqSets = (int **) rax_malloc(sizeof(int*) * numberOfSets);
  int* setsToUniqSets = (int*) rax_malloc(sizeof(int) * numberOfSets);
  int numberOfUniqueSets = 0;
  int dropSetCount = 0;



  //stores the scores for each bips, we are using a bitvector approach (need to scale)
    
  //Legacy Code 
  int bvec_scores = 0;
  
  numberOfUniqueSets = getUniqueDropSets(sets, uniqSets, setsToUniqSets, numberOfSets);


  /*
    Detect initial matchings, we calculate them using bitvectors to represent our bipartitions
  */
  printf("===> Detect initial matchings...\n");
  int vLengthBip = 0;

  //determine the bitVector Length of our bitVector
  if(numberOfBips % MASK_LENGTH == 0)
    vLengthBip = numberOfBips / MASK_LENGTH; 
  else 
    vLengthBip = numberOfBips / MASK_LENGTH + 1;

  //Initialize a bvecScore vector with 0s
  int* bvecScores = (int*)rax_calloc(vLengthBip,sizeof(int));

  //Calculate Initial Matchings and save the result in bvecScores
  detectInitialMatchings(sets, bvecScores, bipsPerTree, numberOfTreesAnalyzed, vLengthBip);

  //Short summary until now:
  // - bipsPerTree consists of all bipartitions per tree
  // - bvecScores is the bitvector setting 1 to all bipartition indices which can score 
  // - taxaPerTree number of taxa per tree
  // - smallTreeTaxaList list of all smalltree->largetree translation arrays


  /*
    Generate useful data structures for calculating and updating scores
  */
  printf("===> Create data structures...\n");  
  //Stores the number of bips per Set and initialize it with 0s
  int* numberOfBipsPerSet = (int*)rax_calloc(numberOfUniqueSets,sizeof(int));

  //Stores all sets which includes this taxa
  int **setsOfTaxa = (int**)rax_malloc((tr->mxtips + 1) *sizeof(int*));
  
  //Now calculate number of bipartitions affected by each unique set
  for(int i = 0; i < numberOfSets; i++) {

    int setindex = setsToUniqSets[i];

    numberOfBipsPerSet[setindex]++;
  }

  //Now using the knowledge of how many bips there are per set, generate an array for each unique dropset containing all bips
  int** bipsOfDropSet = (int**)rax_malloc(sizeof(int*)*numberOfUniqueSets);
  
  //Allocate the space needed for storing all bips
  for(int i = 0; i < numberOfUniqueSets; i++) {

    bipsOfDropSet[i] = (int*)rax_malloc(sizeof(int)*numberOfBipsPerSet[i]); 
  }
  
  printf("==> Initialize the Bips Of Taxa \n");
  //Stores the number of bips each taxa is included (ABC|DE is stored by A,B,C,D and E)
  //It can be calculated by iterating through all trees and adding the taxa 
  int **bipsOfTaxa = (int**)rax_malloc((tr->mxtips + 1) * sizeof(int*));
  int *numberOfBipsPerTaxa = (int*)rax_calloc((tr->mxtips + 1), sizeof(int));
  int *taxaBipsCounter = (int*)rax_calloc((tr->mxtips + 1), sizeof(int));

  //Now add up all
  for (int tree = 0; tree < tr->numberOfTrees; tree++) {

    int* list = smallTreeTaxaList[tree];

    for (int j = 0; j < taxaPerTree[tree]; j++) {

      int taxa = list[j];

      numberOfBipsPerTaxa[taxa] = numberOfBipsPerTaxa[taxa] + bipsPerTree[tree];
    } 
  }

  //Now create dummy arrays inside bipsOfTaxa
  for(int i = 1; i < tr->mxtips+1; i++) {
    bipsOfTaxa[i] = (int*)rax_malloc(sizeof(int)*numberOfBipsPerTaxa[i]);
  }

  printf("==> Storing all bip indices of a certain dropset into an array \n");
  //For checking if all dropsets are iterated
  dropSetCount = 0;
  //Arrays of counter to keep in track
  int* counterOfSet = (int*)rax_malloc(sizeof(int)*numberOfUniqueSets);
  for(int i = 0; i < numberOfUniqueSets; i++) {
    counterOfSet[i] = 0;
  }

  currentBips = 0; //Need to keep in track of the number of bips
  //First iterate through all trees 
  for(int i = 0; i < numberOfTreesAnalyzed; i++ ) {

    //get the correct smallTreeTaxa List
    int* list = smallTreeTaxaList[i];

    //For each bipartition in the tree
    for(int j = 0; j < bipsPerTree[i]; j++) {

      //Look at all bips it is compared too
      int dropSetsPerBip = bipsPerTree[i];

      for(int k = 0; k < dropSetsPerBip; k++){

        int indexOfUniqDropSet = setsToUniqSets[dropSetCount + k];

        int* bips_array = bipsOfDropSet[indexOfUniqDropSet]; 

        //add bipartition j into the bips array of its dropset
        bips_array[counterOfSet[indexOfUniqDropSet]] = currentBips; 

        //increment the internal array index 
        counterOfSet[indexOfUniqDropSet]++;
      }
    //Jump to the next correct dropSetCount!
    dropSetCount = dropSetCount + dropSetsPerBip;

    //now insert the bip into bipsOfTaxa Array
    for(int ix = 0; ix < taxaPerTree[i]; ix++) {

      //get the taxa number
      int stree_Taxa = list[ix];

      //get the bips list of this taxa number
      int* bipsList = bipsOfTaxa[stree_Taxa];

      //now get the position of the biplist and put in our bip index
      bipsList[taxaBipsCounter[stree_Taxa]] = currentBips;

      //increment the counter 
      taxaBipsCounter[stree_Taxa]++;

    }

    //increment currentBips
    currentBips++; 
    }

  }

  /***********************************************************************************/
  /* End RF-OPT Graph Construction */
  /***********************************************************************************/

  /* Short summary :
    sets - array of all dropsets
    uniqSets - array of all unique dropsets
    bipsPerTree - bips per tree
    setsToUniqSets - translates the index of sets to the index of its unique dropset index
    bipsOfDropSets - all bips which disappear when dropset i is pruned
    scores - has all scores between 0 and 1 for the bips (however 0s can be found out by looking at all dropsets with link to dropset 0 (because we sort and it will always be the lowest))  
  */


  /***********************************************************************************/
  /* RF-OPT Initial Score Calculation */
  /***********************************************************************************/

  unsigned int bipsVectorLength;

  /* calculate the bitvector length for bips bitvector */
  if(numberOfBips % MASK_LENGTH == 0)
    bipsVectorLength = numberOfBips / MASK_LENGTH;
  else
    bipsVectorLength = 1 + (numberOfBips / MASK_LENGTH); 

  //Starting from index 1 (because 0 stands for all who already matches)
  //We need a score array saving the scores for each uniqset
  int* rf_score = (int*)rax_calloc(numberOfUniqueSets,sizeof(int));

  printf("==> Calculating the score for the first iteration \n \n");

  //Store all bvecs of all merged and destroyed bipartitions per DropSet 
  int* bvecs_bips = (int*)rax_malloc(sizeof(int)*numberOfUniqueSets);
  int* bvecs_destroyed = (int*)rax_malloc(sizeof(int)*numberOfUniqueSets);


  //Iterate through all sets
  for(int i = 0; i < numberOfUniqueSets; i++) {

  	//Bitvectors of merged and destroyed
  	int bvec_destroyed = 0;

  	int* set = uniqSets[i]; //Get the dropset, first dropset is 0 (if something is matching)

  	//printf(" ==> Analyze Unique DropSet %i \n", i);

  	//We use this data structure to keep track of the to toggled bits
  	int* toggleBits = (int*)rax_calloc(numberOfBips, sizeof(int));

  	//Now iterate through the set
  	int j = 0;

    //Stores the affected bips into a bitvector
    int bvec_bips = 0;

  	while(set[j] != -1) {

  		int taxa = set[j]; //Get the taxa
  		//printf("  Taxa number is %i \n",taxa);

      //Check if set[j] is itself already a set
      int test[2] = {taxa,-1}; 

      //0 if it is not a set, index + 1 otherwise
      int test_index = contains(test, uniqSets, numberOfUniqueSets);

      if(test_index){
        //printf("  It also is in uniqSet %i \n", test_index - 1);
        bvec_bips = getBipsOfDropSet(bvec_bips, (test_index - 1), numberOfBipsPerSet, bipsOfDropSet);

      }

      //Get all bips of this taxa to detect which one will be destroyed
  		int* listOfBips = bipsOfTaxa[taxa]; 

  		//Go through all bipartitions containing this taxa
  		for(int k = 0; k < numberOfBipsPerTaxa[taxa]; k++){

  			int bipindex = listOfBips[k]; //Get the index of the Bipartition

  			int bip = ind_bips[bipindex];

			  //Now analyze this Bipartition

			  //Which tree does this bipartition belongs too?
			  int treenumber = treenumberOfBip[bipindex];

			  //Get the taxonToSmallTree Array of this tree
			  int* stTaxa = taxonToReductionList[treenumber];

			  //Translate the global taxon number it into the local index used by our bips
			  int translated_index = stTaxa[taxa - 1]; //We use taxa - 1 because we start counting at taxa 1 = 0 !

			  //Save the to toggle index into toggleBits vector
			  toggleBits[bipindex] |= 1 << translated_index;

			  //Sort for bits set on one side of the bip and on the other side
			  int leftBits = __builtin_popcount(toggleBits[bipindex] & bip);
			  int rightBits = __builtin_popcount(toggleBits[bipindex]) - leftBits;

			  //Check for the number of bits set in the original bip 
			  int leftBip = __builtin_popcount(bip);
			  int rightBip = taxaPerTree[treenumber] - leftBip;

			  //Subtract the total number of bits set on one side of the bip with the bits we have to toggle
			  int leftBip_after = leftBip - leftBits;
			  int rightBip_after = rightBip - rightBits;

			  //Check if bipartition gets trivial/destroyed due to pruning the taxa and set the bit (representing the bip) which is destroyed
			  if((leftBip_after <= 1) | (rightBip_after <=1)) {

        //Add bips to the bits which represent destroyed bipartitions
        bvec_destroyed = setBit(bvec_destroyed,bipindex);

			  }
			
  		}	

  		j++;

  	}//End iterate through the set

    int penality = 0;
    int score = 0;

    int bvec_mask = 0;
    bvec_mask = setOffSet(bvec_mask, numberOfBips);

    //Bitvector of already matching bips
    int bvec_tmp = 0;
    bvec_tmp = ~bvec_scores & bvec_mask;

    //Penality score are all bitvectors who were matching but is destroyed 
    penality = __builtin_popcount(bvec_destroyed & bvec_tmp);

  	//Now iterate through bipsOfDropSet list and extract all bips which will merge into a bitVector
    bvec_bips = getBipsOfDropSet(bvec_bips, i, numberOfBipsPerSet, bipsOfDropSet);

    //Calculate the bitvectors which remains
    bvec_tmp = ~bvec_destroyed & bvec_mask;

    bvec_tmp = bvec_bips & bvec_tmp;

    score = __builtin_popcount(bvec_scores & bvec_tmp);

    rf_score[i] = score - penality;

    //Save our results for convenience into an array
    bvecs_bips[i] = bvec_bips;
    bvecs_destroyed[i] = bvec_destroyed;


  }//End Score Calculation

  printf("======> Scores:\n");
  for(int i = 0; i < numberOfUniqueSets; i++) {
  	printf("RF Score for %i : %i \n", i, rf_score[i]);
    //printBitVector(bvecs_bips[i]);
    //printBitVector(bvecs_destroyed[i]);
  }

  int maxDropSet = getMax(rf_score, numberOfUniqueSets);
  printf("Max Element is %i \n", maxDropSet);


  /***********************************************************************************/
  /* RF-OPT Create Update Data Structures */
  /***********************************************************************************/


  printf("====> Delete DropSet from all bips and update numbers \n");

  //Create a bitVector to store all deleted taxa
  int bvec_deletedTaxa = 0;

  //Create a bitVector to store all still existing bips
  int bvec_existingBips = 0;

  //Create a bitvector to store deleted dropsets
  int bvec_deletedDropSets = 0;

  //Get the dropset
  int* deleteDropSet = uniqSets[maxDropSet];

  //Store it into a BitVector
  bvec_deletedDropSets = setBit(bvec_deletedDropSets,maxDropSet);

  //Select all bips destroyed by removing this dropset
  int bvec_destroyedBips = bvecs_destroyed[maxDropSet];

  //Select all bips that are now matching
  int bvec_matchingBips = bvecs_bips[maxDropSet];

  //Filter for existent bipartitions
  bvec_existingBips = getExistingBips(bvec_existingBips, numberOfBips, bvec_destroyedBips);

  //Iterate through its taxa
  int iterSet = 0;
  while(deleteDropSet[iterSet] != -1) {

    //Get taxon
    int taxon = deleteDropSet[iterSet];

    //Store the taxon into deletedTaxa BitVector
    bvec_deletedTaxa = setBit(bvec_deletedTaxa, taxon - 1);

    //Check if taxon is inside
    int test[2] = {taxon, -1};

    int index = contains(test, uniqSets, numberOfUniqueSets);

    iterSet++;
  }

  //printBitVector(bvec_existingBips);
  //printBitVector(bvec_deletedTaxa);

  //Update the scores with now matching bips
  bvec_scores = bvec_scores & (~bvec_matchingBips);

  //printBitVector(bvec_scores);

  /* Short summary :
    bvec_existingBips - bitVector of all still existing bips
    bvec_deletedTaxa - bitVector to keep track of destroyed taxa
  */

  /***********************************************************************************/
  /* TODO RF-OPT Update function */
  /***********************************************************************************/

  
  /***********************************************************************************/
  /* End RF-OPT Update function */
  /***********************************************************************************/


  printf("==> Unique Sets: ");
  for(int i = 0; i < numberOfUniqueSets; i++) {
    int j = 0;
    int* set = uniqSets[i];
    while(set[j] > -1) {
      printf("%i ",set[j]);
      j++;
    }
    printf("; ");
  }
  printf("\n");

  printf("Number of total Bips %i \n", numberOfBips);



  //printf("Small Bipartitions?: ");

  for(int i = 0; i < tr->numberOfTrees; i++) {
    unsigned int **sBips = sBipsPerTree[i];

    for(int j = 0; j < bipsPerTree[i]; j++) {
      unsigned int *bip = sBips[j];
    }
  }

  //printf("Ind Bipartitions?: ");


  // printf("Induced Bipartitions: ");

  // printBitVector(ind_bips[0]);
  // printBitVector(ind_bips[1]);
  // printBitVector(ind_bips[2]);
  // printBitVector(ind_bips[3]);
  // printBitVector(ind_bips[4]);
  // printBitVector(ind_bips[5]);
  // printBitVector(ind_bips[6]);


  /***********************************************************************************/
  /* Console Logs for debugging */
  /***********************************************************************************/

  //Printing if
  if(0){


    for(int i = 0; i < numberOfUniqueSets; i++) {
      printf("Bips of Set %i: ", i);
        for(int j = 0; j < numberOfBipsPerSet[i]; j++) {
          int* bips = bipsOfDropSet[i];
          printf("%i ", bips[j]);
        }
      printf("\n");
    }


    printf("Induced Bips! \n");
    // Now checking which dropset would destroy which bipartition 
    for(int i = 0 ; i < numberOfBips; i++) {
      printf("Bip %i is %i \n",i,ind_bips[i]);
    }


    printf("Taxa Names : \n");
    for(int i = 0; i < tr->mxtips + 1; i++) {
      printf("%s ",tr->nameList[i]);
    }
    printf("\n");

    printf("Small Tree Taxa Names 0 : \n");
    for(int i = 0; i < taxaPerTree[0]; i++) {
      int* list = smallTreeTaxaList[0];
      int taxa = list[i];	
      printf("%s ",tr->nameList[taxa]);
    }
    printf("\n");

    printf("Small Tree Taxa Names 1 : \n");
    for(int i = 0; i < taxaPerTree[1]; i++) {
      int* list = smallTreeTaxaList[1];
      int taxa = list[i];	
      printf("%s ",tr->nameList[taxa]);
    }
    printf("\n");

    printf("Small Tree Taxa Names 2 : \n");
    for(int i = 0; i < taxaPerTree[2]; i++) {
      int* list = smallTreeTaxaList[2];
      int taxa = list[i];	
      printf("%s ",tr->nameList[taxa]);
    }
    printf("\n");

    printf("Number of DropSets extracted%i \n",dropSetCount);
    printf("Number of Bips extracted %i \n",currentBips);

    //Testing ...
    printf("Number of Sets is %i \n",numberOfSets);
    printf("Number of Unique Sets is %i \n",numberOfUniqueSets);

    printf("==> Testing bips of unique sets \n");
    for(int i = 0; i < numberOfUniqueSets; i++) {
      printf("Bips of Set %i: ", i);
        for(int j = 0; j < numberOfBipsPerSet[i]; j++) {
          int* bips = bipsOfDropSet[i];
          printf("%i ", bips[j]);
        }
      printf("\n");
    }

    printf("==> Testing bips of taxa \n");
    for(int i = 1; i < tr->mxtips + 1; i++) {
      printf("Bips of Taxa %i: ", i);
        for(int j = 0; j < numberOfBipsPerTaxa[i]; j++) {
        int* bips = bipsOfTaxa[i];
        printf("%i ", bips[j]);
        }
      printf("\n");
    }



  printf("==> Unique Sets: ");
  for(int i = 0; i < numberOfUniqueSets; i++) {
    int j = 0;
    int* set = uniqSets[i];
    while(set[j] > -1) {
      printf("%i ",set[j]);
      j++;
    }
    printf("; ");
  }
  printf("\n");

  printf("==> setsToUniqSets: ");
  for(int i = 0; i < numberOfSets; i++) {
    printf("%i ",setsToUniqSets[i]);
  }
  printf("\n");

  //=== TREE GRAPH CONSTRUCTION ENDS ===
  printf("Scores: ");
  printBitVector(bvec_scores);
  
  printf("BipsPerTree: ");
  for(int foo = 0; foo < tr->numberOfTrees; foo++) {

    printf("%i ",bipsPerTree[foo]);

  } 

  printf("\nInduced Bips: ");
  for(int foo = 0;foo < numberOfBips; foo++) {
	  
    printf("%u ",ind_bips[foo]);
  
  }

  printf("\nSmall Tree Bips: ");
  for(int foo = 0;foo < numberOfBips; foo++) {
  
    printf("%u ",s_bips[foo]);

  }

	printf("\n == Sets == \n");
  for(int fooo = 0; fooo < numberOfSets; fooo++){
    printf("Set %i: ", fooo);
    int i = 0;
    while(sets[fooo][i] > -1) {
	   printf("%i ",sets[fooo][i]);
     i++;
    }
    printf("\n");
  }
  printf("\n");

}//End Printif


  printBothOpen("Number of small trees skipped: %d\n\n", tr->numberOfTrees - numberOfTreesAnalyzed);
  
  printBothOpen("Average RF distance %f\n\n", avgRF / (double)numberOfTreesAnalyzed);
  
  printBothOpen("Large Tree: %i, Number of SmallTrees analyzed: %i \n\n", tr->mxtips, numberOfTreesAnalyzed); 
  
  printBothOpen("Total execution time: %f secs\n\n", gettime() - masterTime);
   
  printBothOpen("File containing all %d pair-wise RF distances written to file %s\n\n", numberOfTreesAnalyzed, rfFileName);

  printBothOpen("execution stats:\n\n");
  printBothOpen("Accumulated time Effective algorithm: %.5f sec \n", sumEffectivetime);
  printBothOpen("Average time for effective: %.10f sec \n",sumEffectivetime / (double)numberOfTreesAnalyzed);
  printBothOpen("Preprocessingtime: %0.5f sec \n\n", preprocessendtime);
 

  fclose(treeFile);
  fclose(rfFile);    
  
  /* free the data structure used for parsing the potentially multi-furcating tree */

  freeMultifurcations(smallTree);
  rax_free(smallTree);

  rax_free(taxonToLabel);
  rax_free(taxonToEulerIndex);
  rax_free(labelToTaxon);
  rax_free(eulerIndexToLabel);
  rax_free(taxonToReduction);
  rax_free(taxonHasDeg);
}


/**********************************************************************************/
/**********************************************************************************/

#else

/************************************* efficient multifurcating algorithm **************************************/

/* Edits compared to the bifurcating algorithm: 
Array storing the degree of every node, called taxonHasDeg,
rec_preOrderTraversalMulti now extracts information about degree of inner nodes,
rec_extractBipartitionsMulti needs two additional parameters: array storing the degree of the nodes and max number of splits */

void plausibilityChecker(tree *tr, analdef *adef)
{
  FILE 
    *treeFile,
    *rfFile;
  
  tree 
    *smallTree = (tree *)rax_malloc(sizeof(tree));

  char 
    rfFileName[1024];

  int
    numberOfTreesAnalyzed = 0,
    i;

  double 
    avgRF = 0.0,
    sumEffectivetime = 0.0;

  /* set up an output file name */

  strcpy(rfFileName,         workdir);  
  strcat(rfFileName,         "RAxML_RF-Distances.");
  strcat(rfFileName,         run_id);

  rfFile = myfopen(rfFileName, "wb");  

  assert(adef->mode ==  PLAUSIBILITY_CHECKER);

  /* open the big reference tree file and parse it */

  treeFile = myfopen(tree_file, "r");

  printBothOpen("Parsing reference tree %s\n", tree_file);

  treeReadLen(treeFile, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);

  assert(tr->mxtips == tr->ntips);
  
  /*************************************************************************************/
  /* Preprocessing Step */

  double 
    preprocesstime = gettime();
  
  /* taxonToLabel[2*tr->mxtips - 2]; 
  Array storing all 2n-2 labels from the preordertraversal: (Taxonnumber - 1) -> (Preorderlabel) */
  int 
    *taxonToLabel  = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int)),

    /* taxonHasDeg[2*tr->mxtips - 2] 
    Array used to store the degree of every taxon, is needed to extract Bipartitions from multifurcating trees 
    (Taxonnumber - 1) -> (degree of node(Taxonnumber)) */

    *taxonHasDeg = (int *)rax_calloc((2*tr->mxtips - 2),sizeof(int)),

    /* taxonToReduction[2*tr->mxtips - 2]; 
  Array used for reducing bitvector and speeding up extraction: 
  (Taxonnumber - 1) -> (0..1 (increment count of taxa appearing in small tree))
  (Taxonnumber - 1) -> (0..1 (increment count of inner nodes appearing in small tree)) */

    *taxonToReduction = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int));
    
  int 
    newcount = 0; //counter used for correct traversals

  /* labelToTaxon[2*tr->mxtips - 2];
  is used to translate between Perorderlabel and p->number: (Preorderlabel) -> (Taxonnumber) */
  int 
    *labelToTaxon = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int));
  
  /* Preorder-Traversal of the large tree */
  preOrderTraversal(tr->start->back,tr->mxtips, tr->start->number, taxonToLabel, labelToTaxon, &newcount);

  newcount = 0; //counter set to 0 to be now used for Eulertraversal

  /* eulerIndexToLabel[4*tr->mxtips - 5]; 
  Array storing all 4n-5 PreOrderlabels created during eulertour: (Eulerindex) -> (Preorderlabel) */
  int* 
    eulerIndexToLabel = (int *)rax_malloc((4*tr->mxtips - 5) * sizeof(int));

  /* taxonToEulerIndex[tr->mxtips]; 
  Stores all indices of the first appearance of a taxa in the eulerTour: (Taxonnumber - 1) -> (Index of the Eulertour where Taxonnumber first appears) 
  is used for efficient computation of the Lowest Common Ancestor during Reconstruction Step
  */
  int*
    taxonToEulerIndex  = (int *)rax_malloc((tr->mxtips) * sizeof(int));

  /* Init taxonToEulerIndex and taxonToReduction */
  int 
    ix;

  for(ix = 0; ix < tr->mxtips; ++ix)    
    taxonToEulerIndex[ix] = -1;    
  
  for(ix = 0; ix < (2*tr->mxtips - 2); ++ix)    
    taxonToReduction[ix] = -1;    


  /* Eulertraversal of the large tree*/
  unrootedEulerTour(tr->start->back,tr->mxtips, eulerIndexToLabel, taxonToLabel, &newcount, taxonToEulerIndex);

  /* Creating RMQ Datastructure for efficient retrieval of LCAs, using Johannes Fischers Library rewritten in C
  Following Files: rmq.h,rmqs.c,rmqs.h are included in Makefile.RMQ.gcc */
  RMQ_succinct(eulerIndexToLabel,4*tr->mxtips - 5);

  double 
    preprocessendtime = gettime() - preprocesstime;

  /* Proprocessing Step End */
  /*************************************************************************************/

  printBothOpen("The reference tree has %d tips\n", tr->ntips);

  fclose(treeFile);
  
  /* now see how many small trees we have */

  treeFile = getNumberOfTrees(tr, bootStrapFile, adef);

  checkTreeNumber(tr->numberOfTrees, bootStrapFile);

  /* allocate a data structure for parsing the potentially mult-furcating tree */

  allocateMultifurcations(tr, smallTree);

  /* loop over all small trees */

  for(i = 0; i < tr->numberOfTrees;  i++)
    {      
      int
	numberOfSplits = readMultifurcatingTree(treeFile, smallTree, adef, TRUE);
      
      if(numberOfSplits > 0)
	{
	  int
	    firstTaxon;           

	  double
	    rec_rf,
	    maxRF;

	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Small tree %d has %d tips and %d bipartitions\n", i, smallTree->ntips, numberOfSplits);    
	  
	  /* compute the maximum RF distance for computing the relative RF distance later-on */
	  
	  /* note that here we need to pay attention, since the RF distance is not normalized 
	     by 2 * (n-3) but we need to account for the fact that the multifurcating small tree 
	     will potentially contain less bipartitions. 
	     Hence the normalization factor is obtained as n-3 + numberOfSplits, where n-3 is the number 
	     of bipartitions of the pruned down large reference tree for which we know that it is 
	     bifurcating/strictly binary */
	  
	  maxRF = (double)(2 * numberOfSplits);
	  
	  /* now get the index of the first taxon of the small tree.
	     we will use this to unambiguously store the bipartitions 
	  */
	  
	  firstTaxon = smallTree->start->number;
	  
	  /***********************************************************************************/
	  /* Reconstruction Step */
	  
	  double 
	    time_start = gettime();
	  
	  /* Init hashtable to store Bipartitions of the induced subtree */
	  /* 
	     using smallTree->ntips instead of smallTree->mxtips yields faster code 
	     e.g. 120 versus 128 seconds for 20,000 small trees on my laptop 
	   */
	  hashtable
	    *s_hash = initHashTable(smallTree->ntips * 4);
	  
	  /* smallTreeTaxa[smallTree->ntips]; 
	     Stores all taxa numbers from smallTree into an array called smallTreeTaxa: (Index) -> (Taxonnumber)  */
	  int* 
	    smallTreeTaxa = (int *)rax_malloc((smallTree->ntips) * sizeof(int));
	  
	  /* counter is set to 0 for correctly extracting taxa of the small tree */
	  newcount = 0; 
	  
	  int 
	    newcount2 = 0;
	  
	  /* seq2[2*smallTree->ntips - 2]; 
	     stores PreorderSequence of the reference smalltree: (Preorderindex) -> (Taxonnumber) */
	  int* 
	    seq2 = (int *)rax_malloc((2*smallTree->ntips - 2) * sizeof(int));
	  /* used to store the vectorLength of the bitvector */
	  unsigned int 
	    vectorLength;
	  
	  /* extract all taxa of the smalltree and store it into an array, 
	     also store all counts of taxa and nontaxa in taxonToReduction */
	  rec_extractTaxa(smallTreeTaxa, taxonToReduction, smallTree->start, smallTree->mxtips, &newcount, &newcount2);
	  
	  rec_extractTaxa(smallTreeTaxa, taxonToReduction, smallTree->start->back, smallTree->mxtips, &newcount, &newcount2);
	  
	  /* counter is set to 0 to correctly preorder traverse the small tree */
	  newcount = 0;
	  
	  /* Preordertraversal of the small tree and save its sequence into seq2 for later extracting the bipartitions, it
	     also stores information about the degree of every node */
	  
	  rec_preOrderTraversalMulti(smallTree->start->back,smallTree->mxtips, smallTree->start->number, seq2, taxonHasDeg, &newcount);
	  
	  /* calculate the bitvector length */
	  if(smallTree->ntips % MASK_LENGTH == 0)
	    vectorLength = smallTree->ntips / MASK_LENGTH;
	  else
	    vectorLength = 1 + (smallTree->ntips / MASK_LENGTH); 
	  
	  unsigned int 
	    **bitVectors = rec_initBitVector(smallTree, vectorLength);
	  
	  /* store all non trivial bitvectors using an subtree approach for the induced subtree and 
	     store it into a hashtable, this method was changed for multifurcation */
	  rec_extractBipartitionsMulti(bitVectors, seq2, newcount,tr->mxtips, vectorLength, smallTree->ntips, 
				       firstTaxon, s_hash, taxonToReduction, taxonHasDeg, numberOfSplits);
	  
	  /* counter is set to 0 to be used for correctly storing all EulerIndices */
	  newcount = 0; 
	  
	  /* smallTreeTaxonToEulerIndex[smallTree->ntips]; 
	     Saves all first Euler indices for all Taxons appearing in small Tree: 
	     (Index) -> (Index of the Eulertour where the taxonnumber of the small tree first appears) */
	  int* 
	    smallTreeTaxonToEulerIndex = (int *)rax_malloc((smallTree->ntips) * sizeof(int));
	  
	  /* seq[(smallTree->ntips*2) - 1] 
	     Stores the Preordersequence of the induced small tree */
	  int* 
	    seq = (int *)rax_malloc((2*smallTree->ntips - 1) * sizeof(int));
	  
	  
	  /* iterate through all small tree taxa */
	  for(ix = 0; ix < smallTree->ntips; ix++) 
	    {        
	      int 
		taxanumber = smallTreeTaxa[ix];
	      
	      /* To create smallTreeTaxonToEulerIndex we filter taxonToEulerIndex for taxa in the small tree*/
	      smallTreeTaxonToEulerIndex[newcount] = taxonToEulerIndex[taxanumber-1]; 
	      
	      /* Saves all Preorderlabel of the smalltree taxa in seq*/
	      seq[newcount] = taxonToLabel[taxanumber-1];
	      
	      newcount++;
	    }
	  
	  /* sort the euler indices to correctly calculate LCA */
	  //quicksort(smallTreeTaxonToEulerIndex,0,newcount - 1);             
	  
	  qsort(smallTreeTaxonToEulerIndex, newcount, sizeof(int), sortIntegers);
	  
	  //printf("newcount2 %i \n", newcount2);      
	  /* Iterate through all small tree taxa */
	  for(ix = 1; ix < newcount; ix++)
	    {  
	      /* query LCAs using RMQ Datastructure */
	      seq[newcount - 1 + ix] =  eulerIndexToLabel[query(smallTreeTaxonToEulerIndex[ix - 1],smallTreeTaxonToEulerIndex[ix])]; 	 
	      
	      /* Used for dynamic programming. We save an index for every inner node:
		 For example the reference tree has 3 inner nodes which we saves them as 0,1,2.
		 Now we calculate for example 5 LCA to construct the induced subtree, which are also inner nodes. 
		 Therefore we mark them as 3,4,5,6,7  */
	      
	      taxonToReduction[labelToTaxon[seq[newcount - 1 + ix]] - 1] = newcount2;
	      
	      newcount2 += 1;
	    }
	  
	  /* sort to construct the Preordersequence of the induced subtree */
	  //quicksort(seq,0,(2*smallTree->ntips - 2));
	  
	  qsort(seq, (2 * smallTree->ntips - 2) + 1, sizeof(int), sortIntegers);
	  
	  /* calculates all bipartitions of the reference small tree and count how many bipartition it shares with the induced small tree */
	  int 
	    rec_bips = rec_findBipartitions(bitVectors, seq,(2*smallTree->ntips - 1), labelToTaxon, tr->mxtips, vectorLength, smallTree->ntips, firstTaxon, s_hash, taxonToReduction);
	  
	  /* Reconstruction Step End */
	  /***********************************************************************************/
	  
	  double 
	    effectivetime = gettime() - time_start;
	  
	  /*
	    if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Reconstruction time: %.10f secs\n\n", effectivetime);
	  */
	  
	  /* compute the relative RF */
	  
	  rec_rf = (double)(2 * (numberOfSplits - rec_bips)) / maxRF;
	  
	  assert(numberOfSplits >= rec_bips);	  	 

	  avgRF += rec_rf;
	  sumEffectivetime += effectivetime;
	  
	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Relative RF tree %d: %f\n\n", i, rec_rf);
	  
	  fprintf(rfFile, "%d %f\n", i, rec_rf);
	  
	  /* free masks and hast table for this iteration */
	  rec_freeBitVector(smallTree, bitVectors);
	  rax_free(bitVectors);
	  
	  freeHashTable(s_hash);
	  rax_free(s_hash);
	  
	  rax_free(smallTreeTaxa);
	  rax_free(seq);
	  rax_free(seq2);
	  rax_free(smallTreeTaxonToEulerIndex);

	  numberOfTreesAnalyzed++;
	}
    }
  
  printBothOpen("Number of small trees skipped: %d\n\n", tr->numberOfTrees - numberOfTreesAnalyzed);
  
  printBothOpen("Average RF distance %f\n\n", avgRF / (double)numberOfTreesAnalyzed);
  
  printBothOpen("Large Tree: %i, Number of SmallTrees analyzed: %i \n\n", tr->mxtips, numberOfTreesAnalyzed); 
  
  printBothOpen("Total execution time: %f secs\n\n", gettime() - masterTime);
   
  printBothOpen("File containing all %d pair-wise RF distances written to file %s\n\n", numberOfTreesAnalyzed, rfFileName);

  printBothOpen("execution stats:\n\n");
  printBothOpen("Accumulated time Effective algorithm: %.5f sec \n", sumEffectivetime);
  printBothOpen("Average time for effective: %.10f sec \n",sumEffectivetime / (double)numberOfTreesAnalyzed);
  printBothOpen("Preprocessingtime: %0.5f sec \n\n", preprocessendtime);
 

  fclose(treeFile);
  fclose(rfFile);    
  
  /* free the data structure used for parsing the potentially multi-furcating tree */

  freeMultifurcations(smallTree);
  rax_free(smallTree);

  rax_free(taxonToLabel);
  rax_free(taxonToEulerIndex);
  rax_free(labelToTaxon);
  rax_free(eulerIndexToLabel);
  rax_free(taxonToReduction);
  rax_free(taxonHasDeg);
}

#endif

#else

/************************************* old slow code below ********************************************************************************/

/* function to extract the bit mask for the taxa that are present in the small tree */

static void setupMask(unsigned int *smallTreeMask, nodeptr p, int numsp)
{
  if(isTip(p->number, numsp))
    smallTreeMask[(p->number - 1) / MASK_LENGTH] |= mask32[(p->number - 1) % MASK_LENGTH];
  else
    {    
      nodeptr 
	q = p->next;

      /* I had to change this function to account for mult-furcating trees.
	 In this case an inner node can have more than 3 cyclically linked 
	 elements, because there might be more than 3 outgoing branches 
	 from an inner node */

      while(q != p)
	{
	  setupMask(smallTreeMask, q->back, numsp);
	  q = q->next;
	}
      
      //old code below 
      //setupMask(smallTreeMask, p->next->back, numsp);	  
      //setupMask(smallTreeMask, p->next->next->back, numsp);      
    }
}


static void newviewBipartitions(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength)
{
  if(isTip(p->number, numsp))
    return;
  {
    nodeptr 
      q = p->next->back, 
      r = p->next->next->back;
    unsigned int       
      *vector = bitVectors[p->number],
      *left  = bitVectors[q->number],
      *right = bitVectors[r->number];
    unsigned 
      int i;           

    while(!p->x)
      {	
	if(!p->x)
	  getxnode(p);
      }

    p->hash = q->hash ^ r->hash;

    if(isTip(q->number, numsp) && isTip(r->number, numsp))
      {		
	for(i = 0; i < vectorLength; i++)
	  vector[i] = left[i] | right[i];	  	
      }
    else
      {	
	if(isTip(q->number, numsp) || isTip(r->number, numsp))
	  {
	    if(isTip(r->number, numsp))
	      {	
		nodeptr tmp = r;
		r = q;
		q = tmp;
	      }	   
	    	    
	    while(!r->x)
	      {
		if(!r->x)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	    	 
	  }
	else
	  {	    
	    while((!r->x) || (!q->x))
	      {
		if(!q->x)
		  newviewBipartitions(bitVectors, q, numsp, vectorLength);
		if(!r->x)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   	    	    	    	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	 
	  }

      }     
  }     
}


/* this function actually traverses the small tree, generates the bit vectors for all 
   non-trivial bipartitions and simultaneously counts how many bipartitions (already stored in the has table) are shared with the big tree
*/

static int bitVectorTraversePlausibility(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h,
					 int *countBranches, int firstTaxon, tree *tr, boolean multifurcating)
{

  /* trivial bipartition */

  if(isTip(p->number, numsp))
    return 0;
  else
    {
      int 
	found = 0;

      nodeptr 
	q = p->next;          

      /* recursively descend into the tree and get the bips of all subtrees first */

      do 
	{
	  found = found + bitVectorTraversePlausibility(bitVectors, q->back, numsp, vectorLength, h, countBranches, firstTaxon, tr, multifurcating);
	  q = q->next;
	}
      while(q != p);
           
      /* compute the bipartition induced by the current branch p, p->back,
	 here we invoke two different functions, depending on whether we are dealing with 
	 a multi-furcating or bifurcating tree.
      */

      if(multifurcating)
	newviewBipartitionsMultifurcating(bitVectors, p, numsp, vectorLength);
      else
	newviewBipartitions(bitVectors, p, numsp, vectorLength);
      
      assert(p->x);      

      /* if p->back does not lead to a tip this is an inner branch that induces a non-trivial bipartition.
	 in this case we need to lookup if the induced bipartition is already contained in the hash table 
      */

      if(!(isTip(p->back->number, numsp)))
	{	
	  /* this is the bit vector to insert into the hash table */
	  unsigned int 
	    *toInsert = bitVectors[p->number];
	  
	  /* compute the hash number on that bit vector */
	  hashNumberType 
	    position = oat_hash((unsigned char *)toInsert, sizeof(unsigned int) * vectorLength) % h->tableSize;	 	 

	  /* each bipartition can be stored in two forms (the two bit-wise complements
	     we always canonically store that version of the bit-vector that does not contain the 
	     first taxon of the small tree, we use an assertion to make sure that all is correct */

	  assert(!(toInsert[(firstTaxon - 1) / MASK_LENGTH] & mask32[(firstTaxon - 1) % MASK_LENGTH]));	 
	  	      
	  /* increment the branch counter to assure that all inner branches are traversed */
	  
	  *countBranches =  *countBranches + 1;	
	  	 	
	  /* now look up this bipartition in the hash table, If it is present the number of 
	     shared bipartitions between the small and the big tree is incremented by 1 */
	   
	  found = found + findHash(toInsert, h, vectorLength, position);	
	}
      return found;
    }
}



////multifurcating trees 

void plausibilityChecker(tree *tr, analdef *adef)
{
  FILE 
    *treeFile,
    *rfFile;
  
  tree 
    *smallTree = (tree *)rax_malloc(sizeof(tree));

  char 
    rfFileName[1024];
 
  /* init hash table for big reference tree */
  
  hashtable
    *h      = initHashTable(tr->mxtips * 2 * 2);
  
  /* init the bit vectors we need for computing and storing bipartitions during 
     the tree traversal */
  unsigned int 
    vLength, 
    **bitVectors = initBitVector(tr, &vLength);
   
  int
    numberOfTreesAnalyzed = 0,
    branchCounter = 0,
    i;

  double 
    avgRF = 0.0;

  /* set up an output file name */

  strcpy(rfFileName,         workdir);  
  strcat(rfFileName,         "RAxML_RF-Distances.");
  strcat(rfFileName,         run_id);

  rfFile = myfopen(rfFileName, "wb");  

  assert(adef->mode ==  PLAUSIBILITY_CHECKER);

  /* open the big reference tree file and parse it */

  treeFile = myfopen(tree_file, "r");

  printBothOpen("Parsing reference tree %s\n", tree_file);

  treeReadLen(treeFile, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);

  assert(tr->mxtips == tr->ntips);

  printBothOpen("The reference tree has %d tips\n", tr->ntips);

  fclose(treeFile);
  
  /* extract all induced bipartitions from the big tree and store them in the hastable */
  
  bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, BIPARTITIONS_RF, (branchInfo *)NULL,
			  &branchCounter, 1, FALSE, FALSE);
     
  assert(branchCounter == tr->mxtips - 3);   
  
  /* now see how many small trees we have */

  treeFile = getNumberOfTrees(tr, bootStrapFile, adef);

  checkTreeNumber(tr->numberOfTrees, bootStrapFile);

  /* allocate a data structure for parsing the potentially mult-furcating tree */

  allocateMultifurcations(tr, smallTree);

  /* loop over all small trees */

  for(i = 0; i < tr->numberOfTrees;  i++)
    {          
      int           
	numberOfSplits = readMultifurcatingTree(treeFile, smallTree, adef, TRUE);

      if(numberOfSplits > 0)
	{
	  unsigned int
	    entryCount = 0,
	    k,
	    j,	
	    *masked    = (unsigned int *)rax_calloc(vLength, sizeof(unsigned int)),
	    *smallTreeMask = (unsigned int *)rax_calloc(vLength, sizeof(unsigned int));

	  hashtable
	    *rehash = initHashTable(tr->mxtips * 2 * 2);

	  double
	    rf,
	    maxRF;

	  int 
	    bCounter = 0,  
	    bips,
	    firstTaxon,
	    taxa = 0;

	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Small tree %d has %d tips and %d bipartitions\n", i, smallTree->ntips, numberOfSplits);    

	  /* compute the maximum RF distance for computing the relative RF distance later-on */
	  
	  /* note that here we need to pay attention, since the RF distance is not normalized 
	     by 2 * (n-3) but we need to account for the fact that the multifurcating small tree 
	     will potentially contain less bipartitions. 
	     Hence the normalization factor is obtained as 2 * numberOfSplits, where numberOfSplits is the number of bipartitions
	     in the small tree.
	  */
	  
	  maxRF = (double)(2 * numberOfSplits);
	  
	  /* now set up a bit mask where only the bits are set to one for those 
	     taxa that are actually present in the small tree we just read */
	  
	  /* note that I had to apply some small changes to this function to make it work for 
	     multi-furcating trees ! */

	  setupMask(smallTreeMask, smallTree->start,       smallTree->mxtips);
	  setupMask(smallTreeMask, smallTree->start->back, smallTree->mxtips);

	  /* now get the index of the first taxon of the small tree.
	     we will use this to unambiguously store the bipartitions 
	  */
	  
	  firstTaxon = smallTree->start->number;
	  
	  /* make sure that this bit vector is set up correctly, i.e., that 
	     it contains as many non-zero bits as there are taxa in this small tree 
	  */
	  
	  for(j = 0; j < vLength; j++)
	    taxa += BIT_COUNT(smallTreeMask[j]);
	  assert(taxa == smallTree->ntips);
	  
	  /* now re-hash the big tree by applying the above bit mask */
	  
	  
	  /* loop over hash table */
	  
	  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
	    {    
	      if(h->table[k] != NULL)
		{
		  entry *e = h->table[k];
		  
		  /* we resolve collisions by chaining, hence the loop here */
		  
		  do
		    {
		      unsigned int 
			*bitVector = e->bitVector; 
		      
		      hashNumberType 
			position;
		      
		      int 
			count = 0;
		      
		      /* double check that our tree mask contains the first taxon of the small tree */
		      
		      assert(smallTreeMask[(firstTaxon - 1) / MASK_LENGTH] & mask32[(firstTaxon - 1) % MASK_LENGTH]);
		      
		      /* if the first taxon is set then we will re-hash the bit-wise complement of the 
			 bit vector.
			 The count variable is used for a small optimization */
		      
		      if(bitVector[(firstTaxon - 1) / MASK_LENGTH] & mask32[(firstTaxon - 1) % MASK_LENGTH])		    
			{
			  //hash complement
			  
			  for(j = 0; j < vLength; j++)
			    {
			      masked[j] = (~bitVector[j]) & smallTreeMask[j];			     
			      count += BIT_COUNT(masked[j]);
			    }
			}
		      else
			{
			  //hash this vector 
			  
			  for(j = 0; j < vLength; j++)
			    {
			      masked[j] = bitVector[j] & smallTreeMask[j];  
			      count += BIT_COUNT(masked[j]);      
			    }
			}
		      
		      /* note that padding the last bits is not required because they are set to 0 automatically by smallTreeMask */	
		      
		      /* make sure that we will re-hash  the canonic representation of the bipartition 
			 where the bit for firstTaxon is set to 0!
		      */
		      
		      assert(!(masked[(firstTaxon - 1) / MASK_LENGTH] & mask32[(firstTaxon - 1) % MASK_LENGTH]));
		      
		      /* only if the masked bipartition of the large tree is a non-trivial bipartition (two or more bits set to 1 
			 will we re-hash it */
		      
		      if(count > 1)
			{
			  /* compute hash */
			  position = oat_hash((unsigned char *)masked, sizeof(unsigned int) * vLength);
			  position = position % rehash->tableSize;
			  
			  /* re-hash to the new hash table that contains the bips of the large tree, pruned down 
			     to the taxa contained in the small tree
			  */
			  insertHashPlausibility(masked, rehash, vLength, position);
			}		
		      
		      entryCount++;
		      
		      e = e->next;
		    }
		  while(e != NULL);
		}
	    }
	  
	  /* make sure that we tried to re-hash all bipartitions of the original tree */
	  
	  assert(entryCount == (unsigned int)(tr->mxtips - 3));
	  
	  /* now traverse the small tree and count how many bipartitions it shares 
	     with the corresponding induced tree from the large tree */
	  
	  /* the following function also had to be modified to account for multi-furcating trees ! */
	  
	  bips = bitVectorTraversePlausibility(bitVectors, smallTree->start->back, smallTree->mxtips, vLength, rehash, &bCounter, firstTaxon, smallTree, TRUE);
	  
	  /* compute the relative RF */
	  
	  rf = (double)(2 * (numberOfSplits - bips)) / maxRF;           
	  
	  assert(numberOfSplits >= bips);

	  assert(rf <= 1.0);
	  
	  avgRF += rf;
	  
	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Relative RF tree %d: %f\n\n", i, rf);

	  fprintf(rfFile, "%d %f\n", i, rf);
	  
	  /* I also modified this assertion, we nee to make sure here that we checked all non-trivial splits/bipartitions 
	     in the multi-furcating tree whech can be less than n - 3 ! */
	  
	  assert(bCounter == numberOfSplits);         
	  
	  /* free masks and hast table for this iteration */
	  
	  rax_free(smallTreeMask);
	  rax_free(masked);
	  freeHashTable(rehash);
	  rax_free(rehash);
	  numberOfTreesAnalyzed++;
	}
    }

  printBothOpen("Number of small trees skipped: %d\n\n", tr->numberOfTrees - numberOfTreesAnalyzed);
  
  printBothOpen("Average RF distance %f\n\n", avgRF / (double)numberOfTreesAnalyzed);

  printBothOpen("Total execution time: %f secs\n\n", gettime() - masterTime);

  printBothOpen("\nFile containing all %d pair-wise RF distances written to file %s\n\n", numberOfTreesAnalyzed, rfFileName);

  

  fclose(treeFile);
  fclose(rfFile);    
  
  /* free the data structure used for parsing the potentially multi-furcating tree */

  freeMultifurcations(smallTree);
  rax_free(smallTree);

  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  
  freeHashTable(h);
  rax_free(h);
}

#endif

