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

//Helper methods
#include "plausibilityFast.h"
#include "rfoptMethods.h"

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
  //int **sets = NULL;

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
  detectInitialMatchings(sets, bvecScores, bipsPerTree, numberOfTreesAnalyzed, vLengthBip); //TODO:Buggy

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

      
    //#define _PRINT_
      
    #ifdef _PRINT_

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

  #endif

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
