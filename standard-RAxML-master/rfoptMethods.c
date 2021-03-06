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

#include "dbg.h" //debug stuff
#include "hashmap.h" //include darray and hashmaps
#include "bipartitionList.h" //legacy code
#include "rfoptMethods.h" //include structures

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

//used for including sortInteger in qsort
#include "plausibilityFast.h"



/********************************** Init functions ******************************/

Bipartition* Bipartition_create(unsigned int* bitvector, int matching, int treenumber) {

    Bipartition* bip = rax_malloc(sizeof(Bipartition));
    check_mem(bip);

    bip->bitvector = bitvector;
    bip->matching = matching;
    bip->treeNumber = treenumber;
    bip->vLength = 0;
    bip->leftSize = 0;
    bip->rightSize = 0;
    bip->destroyed = 0;
    bip->predictDestroyed = 0;
    bip->firstTaxon = 0;

    return bip;

error:
    if(bip) {
        rax_free(bip);
    }

    return NULL;
};


Dropset* Dropset_create(int* dropset) {

    Dropset* drop = rax_malloc(sizeof(Dropset));
    check_mem(drop);

    drop->set = dropset;
    //Create pointer list pointing to bipartitions in the tree
    drop->bipartitions = NULL;
    drop->score = 0; 
    drop->destroyed = 0;

    return drop;

error:
    if(drop) {
        rax_free(drop);
    }

    return NULL;
};


RTaxon* RTaxon_create(int taxonNumber) {

    RTaxon* taxon = rax_malloc(sizeof(RTaxon));
    check_mem(taxon);

    taxon->taxonNumber = taxonNumber;
    //Create pointer list pointing to bipartitions in the tree and the dropsets containing
    //this taxon
    taxon->dropsets = NULL;
    taxon->trees = NULL;
    taxon->deleted = 0; 

    return taxon;

error:
    if(taxon) {
        rax_free(taxon);
    }

    return NULL;

};



/************************************  rf-opt functions *****************************************************************************/


/* method adapted for multifurcating trees, changes are: 
we now need information about the degree of an inner node, because it is not 2 anymore 
we also can have less than n-3 splits and therefore need a new parameter telling us the actual split number */
unsigned int** RFOPT_extractBipartitionsMulti(unsigned int** bitvectors, int* seq, int arraysize, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction, int* taxonHasDegree, int maxSplits)
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
 unsigned int** RFOPT_findAddBipartitions(unsigned int ** bitvectors, int* seq, int arraysize, int* translate, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, hashtable* ind_hash, int* taxonToReduction)
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

/************************************* helper functions ***************************/

/* converts integer into binary representation */
 char *int2bin(int a, char *buffer, int buf_size) {
      buffer += (buf_size - 1);
        int i;
        for (i = 31; i >= 0; i--) {
                *buffer-- = (a & 1) + '0';

                    a >>= 1;
                    }

          return buffer;
}

/* function for built-in quicksort sorting after number of bits */
 int sortBipartitions(const void *a, const void *b) 
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
 int sortSets(const void *a, const void *b)
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
int isSameDropSet(int* a, int* b) {
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
int contains(int* check, int** sets, int numberOfSets) {
  int i;
  for (i = 0; i < numberOfSets; i++) {
    
    int* dropset = sets[i];
    
    if(isSameDropSet(sets[i],check)){
      return (i+1);
    }
  }

  return 0;

}

//Get the index for which array arr is max
int getMax(int* arr, int size) {

  int max = 0;
  int i;
  for (i = 1; i < size; i++) {

    if(arr[i] > arr[max]) {
      max = i;
    }

  }

  return max;
}

/******************************* Bit Manipulation functionality ***********************************/


int setBit(int bitVector, int pos) {
  
  bitVector |= 1 << pos;

  return bitVector;
}

int clearBit(int bitVector, int pos) {
  
  bitVector &= ~(1 << pos);

  return bitVector;
}

int checkBit(int bitVector, int pos) {
  
  int bit = 0;
  
  bit = bitVector & (1 << pos);

  return bit;
}

//Method for setting bits in bitvectors
unsigned int* setBitBV(unsigned int* bitVector, int pos) {

  int vectorIndex = pos / MASK_LENGTH;
  int relativePos = pos % MASK_LENGTH;
  
  bitVector[vectorIndex] = setBit(bitVector[vectorIndex], relativePos);
  
  return bitVector;
   
}

//Use to setup a mask to clear the offset bits
int setOffSet(int mask, int offset) {
  int i;
  for(i = 0; i < offset; i++) {
    mask = setBit(mask, i);
  }

  return mask; 
}

//Use to setup a mask of existing bipartitions
int getExistingBips(int mask, int offset, int bvecs_deletedBips) {

  mask = setOffSet(mask, offset);

  mask = mask & ~bvecs_deletedBips;

  return mask;
}

void printBitVector(int bitVector) {
  char buffer[33];
  buffer[32] = '\0';
  int2bin(bitVector, buffer, 32);
  printf("\n==> BitVector = %s \n", buffer);
}

void printSet(int* set) {
  int i = 0;
  printf("Set: ");
  while(set[i] != -1) {
    printf("%i ", set[i]);
    i++;
  }
  printf("\n");
}

//Takes as input a bitvector and returns a new bitvector OLD
 int getBipsOfDropSet(int bvec_bips, int dropsetNumber, int* numberOfBipsPerSet, int** bipsOfDropSet) {
    //Now iterate through bipsOfDropSet list
    int l;
    for(l = 0; l < numberOfBipsPerSet[dropsetNumber]; l++) {

      //Get the index list of bips for this Drop Set
      int* bips = bipsOfDropSet[dropsetNumber];
      int b_index = bips[l];

      //Generate a bitvector 
      bvec_bips = setBit(bvec_bips, b_index); 
    }

    return bvec_bips;
}

/**********************************************************************************/

int* extractSetFromBitVector(unsigned int* bitvector, int* smallTreeTaxa, unsigned int vLength){


    int numberOfOnes = 0;
    int i;
    for(i = 0; i < vLength; i++) {
      numberOfOnes = numberOfOnes + __builtin_popcount(bitvector[i]);
      //printf("%i - %i ones \n", i, __builtin_popcount(bitvector[i]));
    }

    //plus one because of the terminal number -1 to determine the end of the array
    int* set = rax_malloc((numberOfOnes + 1) * sizeof(int));

    i = 0;

    int numberOfZerosBefore = 0;

    int startVector = 0;

    unsigned int extract = bitvector[startVector];

    int countBits = 0; //to count seen taxa

    while(i < numberOfOnes) {

      if(extract != 0) {

        if(checkBit(extract,0)) {

          //Consider all bitvectors before this vector. I.e. MASK_Length = 32 , 
          //this is the third index of the second bitvector, we have 3 + 2 * 32 as real taxon index
          int pastBits = startVector * MASK_LENGTH;
          set[i] = smallTreeTaxa[pastBits + countBits];
          i++;

        }

        extract = extract >> 1;

        countBits++;

      } else {

        countBits = 0; //reset count bit

        //if bitvector has no 1s anymore, go to the next one
        startVector++;
        extract = bitvector[startVector];

        //printf("numberOfOnes: %i, vL: %i sV: %i i: %i \n", numberOfOnes, vLength, startVector, i);
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
int* mergeSets(int* set1, int* set2) {

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

int* extractSetsFromBitVector(unsigned int* bitvector, unsigned int* bitvector2, int* smallTreeTaxa, unsigned int vLength){

    int* set1; 
    int* set2; 
    int* set;

    set1 = extractSetFromBitVector(bitvector, smallTreeTaxa, vLength);
    set2 = extractSetFromBitVector(bitvector2, smallTreeTaxa, vLength);

    set = mergeSets(set1, set2);

    return set;
}

//Calculate the DropSet given two BitVectors of the same length
int* getDropSetFromBitVectors(unsigned int* indBip, unsigned int* sBip, unsigned int vLength, int treenumber, int* taxaPerTree, int* smallTreeTaxa){

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

  int i;
  //Calculate a logical operation on BitVectors
  for(i = 0; i < vLength; i++) {

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
        int o;       
        for(o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--) {        
          x1_x2[vLength - 1] &= ~mask32[o-1];  
          y1_y2[vLength - 1] &= ~mask32[o-1];
          x1_y2[vLength - 1] &= ~mask32[o-1];
          y1_x2[vLength - 1] &= ~mask32[o-1];     
        }
      }


  for(i = 0; i < vLength; i++) {
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

  //for matching bipartitons we store them like this, we need to allocate memory so the variable exists after function scope
  int* matching = (int*)rax_malloc(2*sizeof(int));
  matching[0] = 0;
  matching[1] = -1;

  //Now we check how we extract the DropSets
  if(selectedCount1 == 0 & selectedCount2 == 0) {
    //This is a matching bipartition!
    set = matching;

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

RTaxon** createRTaxonList(int numberOfTaxa) {

  //initial maxSize
  int maxSize = 600;

  RTaxon** map = (RTaxon**)rax_malloc(sizeof(RTaxon*)*numberOfTaxa+1);
  check_mem(map);
  //TaxonList starts with 1
  int i = 0;
  DArray* trees = NULL;
  DArray* dropsets = NULL;

  for(i = 0; i < numberOfTaxa+1; i++) {

    //Create an RTaxon struct
    map[i] = RTaxon_create(i);

    //Init Tree DArray in RTaxon
    trees = DArray_create(sizeof(int),maxSize);
    map[i]->trees = trees;

    //Init Dropset DArray in RTaxon
    dropsets = DArray_create(sizeof(Dropset),maxSize);
    map[i]->dropsets = dropsets;

  }

  return map;

error:
  
  if(map) {
    rax_free(map);
  }

  return NULL;

}

//Set all initial RTaxon->trees
void initRTaxonList(RTaxon** map, int** smallTreeTaxaList, int numberOfTrees, int* taxaPerTree) {

  RTaxon* taxon = NULL;
  int i = 0;
  int j = 0;
  //iterate through all trees
  for(i = 0; i < numberOfTrees; i++) {
    int* translationArray = smallTreeTaxaList[i];
    //printf("Tree: %i \n", i);

    for(j = 0; j < taxaPerTree[i]; j++) {
      int globalTaxonNumber = translationArray[j];
      //printf("Taxon: %i \n", globalTaxonNumber);

      taxon = map[globalTaxonNumber];

      int* val1 = DArray_new(taxon->trees);

      *val1 = i;

      DArray_push(taxon->trees, val1);

    }
  }
}


//uses the first taxon to create unambiguous bitvectors. First taxon is set to 0.
//As input we put the index of the bit, which should serve as unique identifier (aka firstIndex)
static unsigned int* createUniqueKey(unsigned int* bitVector, int firstIndex, int vLength, int ntips) { 
  
  int o = 0;
  int k = 0;

  //for checking purposes 
  int firstTaxon = firstIndex + 1;   

  /* if bitvector contains first taxon, use its complement */
  if(bitVector[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]) {                 
    /* Padding the last bits! */
    if(ntips % MASK_LENGTH != 0) {            
      for(o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--) {        
        bitVector[vLength - 1] |= mask32[o-1];       
      }
    }
    
    for(k=0;k < vLength; k++) {         
      bitVector[k] = ~bitVector[k];
    }         
  }

  return bitVector;
}

static int vHashLength = 0;


/***************************** hash stuff **************************************/ 

//Set the private variable HashLength for use in the compareBitVectors function
void setHashLength(int vLength) {
  vHashLength = vLength;

}

//Iterate through all elements of the array. Because we need to know the length of the array we set a private variable
static int compareBitVectors(unsigned int* a, unsigned int* b) {
  int i = 0;

  while(i < vHashLength) {

    if(a[i] != b[i]) {
      return 1;
    }

    i++;
  }

  //Fall through
  return 0;
}

//returns 0 if same and 1 otherwise
static int compareHash(void* a, void* b) {

  return compareBitVectors((unsigned int*)a, (unsigned int*)b);

}

/************************ hash stuff ends ***************************************/

//Calculates all dropsets of two given bipartition lists, we don't have to create a unique representation :)
void calculateDropSets(RTaxon** taxonList, Hashmap** mapArray, Hashmap* map, unsigned int*** indBipsPerTree, unsigned int*** sBipsPerTree, int** sets, int** smallTreeTaxaList, int* bipsPerTree, 
  int* taxaPerTree, unsigned int* vectorLengthPerTree, int numberOfTrees) {

  int countSets = 0;
  int countBips = 0;
  int uniqueDropSetCount = 0;
  int i;
  //First iterate through all trees
  for(i = 0; i< numberOfTrees; i++) {

    //Get all induced Bips of this tree
    unsigned int **indBips = indBipsPerTree[i];

    //Get all small Bips of this tree
    unsigned int **sBips = sBipsPerTree[i];

    //Create a Hashmap for each tree storing the smallBips using bipartition hash functions
    Hashmap* treeMap = Hashmap_create(compareHash, NULL);

    mapArray[i] = treeMap;
    int j;
    //Now go through all Bips of this tree
    for(j = 0; j < bipsPerTree[i]; j++) {

      //Get the bitVector of this Bips
      unsigned int* indBip = indBips[j];
      int k;
      //Now iterate through all small Bips in the tree
      for(k = 0; k < bipsPerTree[i]; k++) {
        
        //get small Bip
        unsigned int *sBip = sBips[k];

        //extract DropSets from the comparision
        sets[countSets] = getDropSetFromBitVectors(indBip, sBip, vectorLengthPerTree[i], i, taxaPerTree, smallTreeTaxaList[i]);

        //================== Dropset generation ============================//

        //Our key for the hashtable is the sorted dropset taxa for this bip
        int* key = sets[countSets];

        //If the dropset is 0, the bipartition is initally matching
        int matching = 0; //matching initally 0
        if(key[0] == 0){
          matching = 1;
        }  

        //Get the existing dropset
        Dropset* drop = Hashmap_get(map, key);

        //Declare bips
        DArray* bips = NULL;

        //If dropset doesn't exist we create a dropset struct 
        //Else we just push another bip into the dropset bipartitions list
        if(drop == NULL) {

          drop = Dropset_create(key);

          //Create an DArray which is able to store bips structures
          bips = DArray_create(sizeof(Bipartition),taxaPerTree[i]);

          //Dropset now saves a darray of bips
          drop->bipartitions = bips;

          assert(0 == Hashmap_set(map, drop->set, drop));

          uniqueDropSetCount++;

        } else {

          //Get the existing dropset
          drop = Hashmap_get(map, drop->set);

          //Get the existing bipartitions list
          bips = drop->bipartitions;

        }

        //Create a pointer which will be given to the dropset
        Bipartition* bip = NULL; 

        //Check if we already have it in the hashmap
        Bipartition* res = Hashmap_get(treeMap,sBip);

        int sBipLength = vectorLengthPerTree[i];
        
        //Set vHashLength for HashComparision
        setHashLength(sBipLength);
        
        //================== Tree Hashtable generation =======================//

        //This is done only once per bitvector
        if(res == NULL) { 
          //Create bipartition stuct only when the bip is not inside the hashtable
          bip = Bipartition_create(sBip,matching,i);

          int l = 0;
          int leftCount = 0;
          int rightCount = 0;
          for(l = 0; l < sBipLength; l++) {
            leftCount = leftCount + __builtin_popcount(sBip[l]);
          }

          //store popcount into the left side
          bip->leftSize = leftCount;

          //leftover bits are from right side
          rightCount = taxaPerTree[i] - leftCount;
          bip->rightSize = rightCount;
          //saving vector length
          bip->vLength = sBipLength;

          //Set bipartition if its not in hashtable
          Hashmap_set(treeMap,bip->bitvector,bip);
          //printf("new bip set for tree %i \n", i);
          countBips++;

        } else {
          bip = res;
          //check if matching was set to 1 and set it
          if(matching) { 
            bip->matching = 1;
          }
        }

        //Push bipartition into array (not sure if we don't need a control var for double)
        assert(0 == DArray_push(bips,bip));

        //================== RTaxon Dropset ==================================//
        //Saves a list of dropsets 
        int l = 0;

        while(key[l] != -1) { 
          int globalTaxonNumber = key[l];
          RTaxon* taxon = taxonList[globalTaxonNumber];
          DArray_push(taxon->dropsets,drop); 
          l++;
        }

        countSets++;
        if(countSets % 5000 == 0) {
          printf("looked at %i sets, current unique dropsets: %i \n", countSets, uniqueDropSetCount);
        }
        //printf("%i \n",countSets);

      }
    }
  }

  printf("\n=> Unique Dropsets found: %i \n", uniqueDropSetCount);
}

//Search for the first unset bit
static int findIndexOfFirstZero(unsigned int* bitVector, int vLength) {
  int i = 0;
  int j = 0;
  int index = 0;
  //Iterate through all bits
  for(i = 0; i < vLength; i++) {
    for(j = 0; j < MASK_LENGTH; j++) {
      //check if its not setted
      if(!checkBit(bitVector[i],j)) {
        index = (MASK_LENGTH * i + j);
        return index;
      }
    }
  }

  return -1; //sth went wrong, couldnt find a zero
}

//Create BitVectors
unsigned int** createBitVectors(int numberOfTrees, unsigned int* vectorLengthPerTree) {
  
  unsigned int** bitVectors = (unsigned int**)rax_malloc(sizeof(unsigned int*) * numberOfTrees);
  
  int i = 0;
  
  for(i = 0; i < numberOfTrees; i++)    
    bitVectors[i] = (unsigned int*)rax_calloc(vectorLengthPerTree[i], sizeof(unsigned int));
        
  return bitVectors;
}

//calculates which bips are going to be destroyed using rehashing and size of the bips
//special care for deleting the first taxon - we have to choose a new taxon for unambigious representation
static void removeTaxonFromTree(unsigned int** deletedTaxa, int treeNumber, Hashmap** mapArray, int vLength, int ntips) {

  int i = 0;
  int j = 0;
  int k = 0;
  //Create new hashtable to check for merging bips
  Hashmap* testMap = NULL;

  //Hashtable for comparing bitvectors
  testMap = Hashmap_create(compareHash,NULL);

  //Get the tree hashmap
  Hashmap* treeMap = mapArray[treeNumber];

  //get the deleted taxons
  unsigned int* deletedTaxaList = deletedTaxa[treeNumber];

  //check if we deleted the firstTaxon (index 0)
  int deletedFirstTaxon = checkBit(deletedTaxaList[0],0);

  int nextFirstTaxonIndex = 0;
  //if we delete the first Taxon, we search a new one
  if(deletedFirstTaxon) {
    nextFirstTaxonIndex = findIndexOfFirstZero(deletedTaxaList,vLength);

  }

  //iterate through bips in tree
  for(i = 0; i < DArray_count(treeMap->buckets); i++) {
        DArray* bucket = DArray_get(treeMap->buckets, i);
        if(bucket) {
            for(j = 0; j < DArray_count(bucket); j++) {
                HashmapNode* node = DArray_get(bucket, j);

                Bipartition* bip = node->data;
                
                //if bip is not destroyed
                if(!(bip->destroyed)){

                  bip->predictDestroyed = 0; 
                  //We assume it won't be destroyed -- if this bip survives all checks it falls through

                  //total number of deleted taxa
                  int numberOfDeletedTaxa = 0;
                  int newLeftSize = 0;
                  int vLength = bip->vLength;
                  //get the deleted taxa
                  unsigned int* res = (unsigned int*)rax_malloc(vLength * sizeof(unsigned int));

                  //used to keep track of original bitVectors left and right size
                  unsigned int* resCount = (unsigned int*)rax_malloc(vLength * sizeof(unsigned int));

                  unsigned int* bitvector = bip->bitvector;

                  memcpy(res,bitvector,(vLength)*sizeof(unsigned int));

                  //Check if we destroyed the firstTaxon
                  if(deletedFirstTaxon) {
                    //We create a new unambigious representation using the next found taxon in res
                    res = createUniqueKey(res,nextFirstTaxonIndex,bip->vLength,ntips);
                      
                  }

                  //calculate how much taxa are deleted on the left side
                  for(k = 0; k < bip->vLength; k++) {

                    numberOfDeletedTaxa = numberOfDeletedTaxa + __builtin_popcount(deletedTaxaList[k]);
                    
                    //Delete taxa in the new unambigious representation
                    res[k] = res[k] & ~(deletedTaxaList[k]);

                    //Use this to keep track of the sizes
                    resCount[k] = bitvector[k] & ~(deletedTaxaList[k]);
                    newLeftSize = newLeftSize + __builtin_popcount(resCount[k]);

                  }
                  
                  assert(bip->leftSize >= newLeftSize);

                  int leftDeletedSize = bip->leftSize - newLeftSize;

                  int rightDeletedSize = 0;

                  //Simple math to figure our how many taxa are subtracted from the right side
                  rightDeletedSize = numberOfDeletedTaxa - leftDeletedSize;

                  //Get the new size 
                  int newRightSize = bip->rightSize - rightDeletedSize;

                  //Check if destroyed
                  if(newRightSize < 2 || newLeftSize < 2) {
                    bip->predictDestroyed = 1;

                  } else {  
                    //Now rehash and check if its already there
                    setHashLength(bip->vLength);

                    //we look for the new reduced bitvector
                    //due to each tree and umabigious representation, if we delete taxa they still stay unambigious
                    unsigned int* val = Hashmap_get(testMap,res);
                    
                    //not yet in hashtable testMap, so we rehash
                    if(val == NULL) {
                      Hashmap_set(testMap,res,res);
                    } else {
                      //already inside! So we set this bipartition to destroyed
                      bip->predictDestroyed = 1;
                    }

                  }

                } 


            }
        }
    }

  //TODO: free testMap;
  Hashmap_destroy(testMap);
}

//translate from global to local index
static int getLocalIndex(int** taxonToReductionList, int treeNumber, int globalTaxonNumber) {
  //Get the converter of the tree i  
  int* localTaxa = taxonToReductionList[treeNumber];

  //Translate the global taxon number it into the local index used by our bips
  //We use taxa - 1 because we start counting for localTaxa at taxa 1 = 0 !
  int localIndex = localTaxa[globalTaxonNumber - 1];

  return localIndex; 
}

//Copy bitvectors
static unsigned int** copyBitVectors(int numberOfTrees, unsigned int* vectorLengthPerTree, unsigned int** originalBitVectors) {

  unsigned int** bitVectors = (unsigned int**)rax_malloc(sizeof(unsigned int*) * numberOfTrees);
  
  int i = 0;
  
  //Copy all bits
  for(i = 0; i < numberOfTrees; i++) {  

    bitVectors[i] = (unsigned int*)rax_malloc(vectorLengthPerTree[i] * sizeof(unsigned int));
    memcpy(bitVectors[i],originalBitVectors[i],vectorLengthPerTree[i] * sizeof(unsigned int));

  }
        
  return bitVectors;  
}

//Function for calculating the predicted penality
static int calculateLoss(int numberOfTrees, Hashmap** mapArray) {

  int loss = 0;
  int i = 0;

  for(i = 0; i < numberOfTrees; i++) {
    //printf("tree %i \n",i);
    Hashmap* treeHash = mapArray[i];
    int k = 0;
    int j = 0;
    for(k = 0; k < DArray_count(treeHash->buckets); k++) {
      DArray* bucket = DArray_get(treeHash->buckets,k);
      if(bucket) {
        
        for(j = 0; j < DArray_count(bucket); j++) {
          HashmapNode* node = DArray_get(bucket, j);

          Bipartition* bip = node->data;
                
          //if it is not destroyed yet
          if(bip->destroyed == 0) {
            //if it is getting destroyed ..
            if(bip->predictDestroyed == 1) {
              //and it was matching before ..
              if(bip->matching == 1) { 
                //we increment loss!
                loss++;
              }
            }
          }

        }
      }
    }
  }

  return loss;
}

//Set deletedTaxaCopy as array containing the bitvectors where deleted taxa are set to one and encoded with local indices
unsigned int** setDeletedBitVectors(unsigned int** deletedTaxaCopy, int* taxaList, 
  RTaxon** RTaxonList, int** taxonToReductionList) {

  RTaxon* rtaxon = NULL;
  int i = 0;
  int j = 0;
  int globalIndex = 0;
  int treeNumber = 0;
  int localIndex = 0;
  unsigned int* deletedTaxaTree = NULL; 

  while(taxaList[i] != -1) {
    globalIndex = taxaList[i];
    rtaxon = RTaxonList[globalIndex];

    DArray* trees = rtaxon->trees;

    for(j = 0; j < DArray_count(trees); j++) {

      int* tree = DArray_get(trees,j);
      treeNumber = *tree; 

      //translate global index into local index
      localIndex = getLocalIndex(taxonToReductionList, treeNumber, globalIndex);

      //Get the copied deleted taxa list for each tree
      deletedTaxaTree = deletedTaxaCopy[treeNumber];

      //Set bit at position localIndex to 1 on the copied version
      setBitBV(deletedTaxaTree, localIndex);
    }

    i++;
  }

  return deletedTaxaCopy;
}

//Predict the score
int Dropset_score(Dropset* drop, RTaxon** RTaxonList, unsigned int** deletedTaxa, Hashmap** mapArray, 
  int** taxonToReductionList, int numberOfTrees, unsigned int* vectorLengthPerTree, int* taxaPerTree) {

  int i = 0;
  int j = 0;
  int scoreGain = 0;
  int scorePenalty = 0;
  int treeNumber = 0;
  int localIndex = 0;

  int globalIndex = 0;
  RTaxon* rtaxon = NULL;

  Hashmap* predict = NULL;
  unsigned int* deletedTaxaTree = NULL;
  unsigned int** deletedTaxaCopy = NULL;

  //Get taxa
  int* taxaList = drop->set;

  //Performance: Maybe only copy trees we need instead of all?
  deletedTaxaCopy = copyBitVectors(numberOfTrees, vectorLengthPerTree, deletedTaxa);
  //-1 is ending symbol for this array
  deletedTaxaCopy = setDeletedBitVectors(deletedTaxaCopy, taxaList, RTaxonList, taxonToReductionList);

  //Now iterate through all trees and remove the taxons
  for(i = 0; i < numberOfTrees; i++) {

    removeTaxonFromTree(deletedTaxaCopy, i, mapArray, vectorLengthPerTree[i], taxaPerTree[i]); 
  }

  //Traverse through all possible RF gains 
  for(i = 0; i < DArray_count(drop->bipartitions); i++) {
    Bipartition* bip = DArray_get(drop->bipartitions, i);
    
    //If it is not destroyed before
    if(bip->destroyed == 0) {
      
      if(bip->matching == 0) {
        //if it was not matching before..
        if(bip->predictDestroyed == 0) {
          scoreGain++; //..and it its not destroyed, we have a score gain
        }
      }
    }
  }

  //Now traverse through all bips and see if we would destroy any matching
  scorePenalty = calculateLoss(numberOfTrees,mapArray);

  int totalScore = scoreGain - scorePenalty;
  //TODO: We need to free copied bitvectors!

  return totalScore;
}












