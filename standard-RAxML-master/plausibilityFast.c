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


#include "bipartitionList.h" //legacy code


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


/************************************  new fast code by David Dao *****************************************************************************/

/* Calculates the size of every subtree and extract its bipartition by seperating the subtree from the small tree 
 This Algorithm works on the induced bifurcating subtree only and needs therefore no multifurcation adaption 
 It then counts, how many bipartition is shared with the reference small tree*/
int rec_findBipartitions(FILE * file, unsigned int ** bitvectors, int* seq, int arraysize, int* translate, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction)
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

  fprintf(file, "ind_bip_start\n");
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

      // Save the small bitvectors
      for(k=0; k < vLength; k++)
        fprintf(file,"%u ",toInsert[k]);  
      fprintf(file,"\n");  
      
      found = found + findHash(toInsert, hash, vLength, position);              
    }

    fprintf(file, "ind_bip_end\n");

  rax_free(V);
  rax_free(bipartitions);
  
  return found;
}

/* method adapted for multifurcating trees, changes are: 
we now need information about the degree of an inner node, because it is not 2 anymore 
we also can have less than n-3 splits and therefore need a new parameter telling us the actual split number */
void rec_extractBipartitionsMulti(FILE *file, unsigned int** bitvectors, int* seq, int arraysize, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction, int* taxonHasDegree, int maxSplits)
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
  fprintf(file, "s_bip_start\n");
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

      //Save the bitvector for the bip from left to right
      
      for(k=0; k < vLength; k++)
        fprintf(file,"%u ",toInsert[k]);  
      fprintf(file,"\n");  
    }
    fprintf(file, "s_bip_end\n");

  rax_free(V);
  rax_free(bipartitions);  
}




/*Preordertraversal of the big tree using bitVectorInitrav as reference and taking start->back node, 
number of tips and start->number as parameter and delivers a TaxonToPreOrderLabel and LabelToTaxon Array*/
void preOrderTraversal(nodeptr p, int numsp, int start, int* array, int* backarray, int* pos)
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
void rec_extractTaxa(int* smallTreeTaxa, int* taxonToReduction, nodeptr p, int numsp, int* pos, int* pos2)
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
void rec_preOrderTraversalMulti(nodeptr p, int numsp, int start, int* backarray, int* deg, int* pos)
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
unsigned int **rec_initBitVector(tree *tr, unsigned int vectorLength)
{
  unsigned int 
    **bitVectors = (unsigned int **)rax_malloc(sizeof(unsigned int*) * 2 * tr->ntips - 1);
  
  int 
    i;
  
  for(i = 0; i < (2 * tr->ntips - 1); i++)    
    bitVectors[i] = (unsigned int *)rax_calloc(vectorLength, sizeof(unsigned int));
        
  return bitVectors;
}


void rec_freeBitVector(tree *tr, unsigned int **bitVectors)
{  
  int 
    i;
  
  for(i = 0; i < (2 * tr->ntips - 1); i++)    
    rax_free(bitVectors[i]);
}

/*euler traversal for binary and rooted trees*/                                                                                                                                                              
void eulerTour(nodeptr p, int numsp, int* array, int* reference, int* pos, int* taxonToEulerIndex)
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
void unrootedEulerTour(nodeptr p, int numsp, int* array, int* reference, int* pos, int* taxonToEulerIndex)
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
int sortIntegers(const void *a, const void *b)
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

