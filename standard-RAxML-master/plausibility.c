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


#define _USE_FAST_PLAUSIBILITY

#ifdef _USE_FAST_PLAUSIBILITY

/************************************  new fast code by David Dao *****************************************************************************/

#include "plausibilityFast.h"

/************************************* rf-opt functions *****************************/


/**********************************************************************************/
/*************************** RF-OPT Algorihm **************************************/
/**********************************************************************************/

//#define _USE_RF_OPT

#ifdef _USE_RF_OPT

#include "rfoptAlgo.h"

//Seperating plausibility analysis from the original plausibility checker
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
    *rfFile,
    *nameFile,
    *i_bipsFile;
  
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


  /* this file contains mask length, then tree i, taxa number, splits*/
  i_bipsFile = myfopen("ind_bips.txt","wb");
  /* this file contains the name vector of the large tree */

  fprintf(i_bipsFile,"MASK_LENGTH %i\n",MASK_LENGTH);

  printBothOpen("Parsing reference tree %s\n", tree_file);

  treeReadLen(treeFile, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);

  assert(tr->mxtips == tr->ntips);

  nameFile = myfopen("names.txt","wb");

  int nameCount = 0;

  //saving the tree taxon names into file names.txt
  for (nameCount = 0; nameCount < tr->mxtips + 1; nameCount++) {
  	fprintf(nameFile, "%s ", tr->nameList[nameCount+1]);
  }
  
  fclose(nameFile);

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

	  /* writes information to ind_bips.txt */
	  fprintf(i_bipsFile,"Tree\n%d %d %d\n", i, smallTree->ntips, numberOfSplits);    
	  
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
	  
	  /* store all non trivial bitvectors using an subtree approach for the small subtree and 
	     store it into a hashtable, this method was changed for multifurcation */
	  rec_extractBipartitionsMulti(i_bipsFile, bitVectors, seq2, newcount,tr->mxtips, vectorLength, smallTree->ntips, 
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
		/* write small tree taxa array into ind_bip.txt*/
		fprintf(i_bipsFile,"%i ",smallTreeTaxa[ix]);
	      
	      /* To create smallTreeTaxonToEulerIndex we filter taxonToEulerIndex for taxa in the small tree*/
	      smallTreeTaxonToEulerIndex[newcount] = taxonToEulerIndex[taxanumber-1]; 
	      
	      /* Saves all Preorderlabel of the smalltree taxa in seq*/
	      seq[newcount] = taxonToLabel[taxanumber-1];
	      
	      newcount++;
	    }
	    fprintf(i_bipsFile,"\n");

	  // for(ix = 0; ix < smallTree->ntips; ix++)
	  //   {
	  //     int
	  //   _taxanumber = smallTree[ix];
	  //   fprintf(i_bipsFile,"%i",taxonToReduction[ix])
	  //   }
	  
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
	    rec_bips = rec_findBipartitions(i_bipsFile, bitVectors, seq,(2*smallTree->ntips - 1), labelToTaxon, tr->mxtips, vectorLength, smallTree->ntips, firstTaxon, s_hash, taxonToReduction);
	  
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
 
  fclose(i_bipsFile);
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

