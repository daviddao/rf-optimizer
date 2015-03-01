/************************************ data structs ************************************/

typedef struct Dropset {
	int* set; //dropset
	int score; //score for this dropset
	DArray* bipartitions; //Affected bipartitions

} Dropset;

//Each tree has a own hashmap 
//structure to store bipartitions
//key: deleted set DArray 
//if bipartition - deleted 
typedef struct Bipartition {
	unsigned int* bitvector; //stores bitvector of smalltree
	int matching; //matching 1, not matching 0
	int treenumber; //stores number of tree where bip is from

} Bipartition;

//Struct storing information of dropsets and trees where this taxon occurs.
//These datastructures have to be modified when deleting this taxon.
typedef struct RTaxon {
	int taxonNumber;
	DArray* dropsets;
	DArray* trees;
	int deleted;

} RTaxon;


/********************************** Init functions ******************************/

Bipartition* Bipartition_create(unsigned int* bitvector, int matching, int treenumber);

Dropset* Dropset_create(int* dropset);

RTaxon* RTaxon_create(int taxonNumber);

/************************************  rf-opt functions *****************************************************************************/

/* method adapted for multifurcating trees, changes are: 
we now need information about the degree of an inner node, because it is not 2 anymore 
we also can have less than n-3 splits and therefore need a new parameter telling us the actual split number */
unsigned int** RFOPT_extractBipartitionsMulti(unsigned int** bitvectors, int* seq, int arraysize, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction, int* taxonHasDegree, int maxSplits);


/* Additionally to the method above, rec_findAddBipartitions add the bipartitions into a second hashtable ind_hash and returns the bitVector of the induced tree */
unsigned int** RFOPT_findAddBipartitions(unsigned int** bitvectors, int* seq, int arraysize, int* translate, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, hashtable* ind_hash, int* taxonToReduction);


/* 
Extract the set of the bitvector and stores the taxa into an array called set which it returns
TODO only for one bitvector! Stop Element is -1 
*/
 int* extractSet(int* bitvector, int* smallTreeTaxa);

/* 
Can extract two bipartitons and merge their sets
Edit compared to extractSet: for bitvector 2 we restart the loop at position set[i + numberOfOnesBip1]
*/
int* extractSets(int* bitvector, int* bitvector2, int* smallTreeTaxa);

/************************************* helper functions ***************************/

/* converts integer into binary representation */
char *int2bin(int a, char *buffer, int buf_size);

/* function for built-in quicksort sorting after number of bits */
int sortBipartitions(const void *a, const void *b); 

/* sort multidimensional arrays with different size */
int sortSets(const void *a, const void *b);

/* Checks if two arrays are identical and returns 1 and 0 */
int isSameDropSet(int* a, int* b);

/* Checks if check already is inside sets between 0 ... numberOfSets 
    returns 0 if its not containing
    returns index+1 if it contains the element 
*/
int contains(int* check, int** sets, int numberOfSets);

//Get the index for which array arr is max
int getMax(int* arr, int size);

/******************************* Bit Manipulation functionality ***********************************/

int setBit(int bitVector, int pos);

int clearBit(int bitVector, int pos);

int checkBit(int bitVector, int pos);

//Method for setting bits in bitvectors
int* setBitBV(int* bitVector, int pos, int vLength);

//Use to setup a mask to clear the offset bits
int setOffSet(int mask, int offset);

//Use to setup a mask of existing bipartitions
int getExistingBips(int mask, int offset, int bvecs_deletedBips);

void printBitVector(int bitVector);

void printSet(int* set);

//Takes as input a bitvector and returns a new bitvector OLD
int getBipsOfDropSet(int bvec_bips, int dropsetNumber, int* numberOfBipsPerSet, int** bipsOfDropSet);

/**********************************************************************************/

int* extractSetFromBitVector(unsigned int* bitvector, int* smallTreeTaxa, unsigned int vLength);

//Merge two dropsets into one set (ending element is -1)
int* mergeSets(int* set1, int* set2);

int* extractSetsFromBitVector(unsigned int* bitvector, unsigned int* bitvector2, int* smallTreeTaxa, unsigned int vLength);

//Calculate the DropSet given two BitVectors of the same length
int* getDropSetFromBitVectors(unsigned int* indBip, unsigned int* sBip, unsigned int vLength, int treenumber, int* taxaPerTree, int* smallTreeTaxa);

//Filter for unique dropsets and return number of unique sets
int getUniqueDropSets(int** sets, int** uniqSets, int* setsToUniqSets, int numberOfSets);

//Calculates all dropsets of two given bipartition lists
void calculateDropSets(Hashmap** mapArray, Hashmap* map, unsigned int*** indBipsPerTree, unsigned int*** sBipsPerTree, int** sets, int** smallTreeTaxaList, int* bipsPerTree, 
  int* taxaPerTree, unsigned int* vectorLengthPerTree, int numberOfTrees);

RTaxon** createRTaxonList(int numberOfTaxa);

void initRTaxonList(RTaxon** map, int** smallTreeTaxaList, int numberOfTrees, int* taxaPerTree);

void detectInitialMatchings(int** sets, int* matchingVector, int* bipsPerTree, int numberOfTrees,  int vLength);
