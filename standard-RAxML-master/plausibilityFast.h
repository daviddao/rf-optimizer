/************************************  new fast code by David Dao *****************************************************************************/

/* Calculates the size of every subtree and extract its bipartition by seperating the subtree from the small tree 
 This Algorithm works on the induced bifurcating subtree only and needs therefore no multifurcation adaption 
 It then counts, how many bipartition is shared with the reference small tree*/
int rec_findBipartitions(unsigned int ** bitvectors, int* seq, int arraysize, int* translate, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction);

/* method adapted for multifurcating trees, changes are: 
we now need information about the degree of an inner node, because it is not 2 anymore 
we also can have less than n-3 splits and therefore need a new parameter telling us the actual split number */
void rec_extractBipartitionsMulti(unsigned int** bitvectors, int* seq, int arraysize, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction, int* taxonHasDegree, int maxSplits);

/*Preordertraversal of the big tree using bitVectorInitrav as reference and taking start->back node, 
number of tips and start->number as parameter and delivers a TaxonToPreOrderLabel and LabelToTaxon Array*/
void preOrderTraversal(nodeptr p, int numsp, int start, int* array, int* backarray, int* pos);

/*extract all smalltree taxa and store a list of all Taxon*/
void rec_extractTaxa(int* smallTreeTaxa, int* taxonToReduction, nodeptr p, int numsp, int* pos, int* pos2);

/* traverses the reference small tree and additionally extracting for every node its degree. It is stored in int* deg */
void rec_preOrderTraversalMulti(nodeptr p, int numsp, int start, int* backarray, int* deg, int* pos);

/* special function inits bitvectors to store bipartitions for dynamic programming*/
unsigned int **rec_initBitVector(tree *tr, unsigned int vectorLength);

/* free bitvector */
void rec_freeBitVector(tree *tr, unsigned int **bitVectors);

/*euler traversal for binary and rooted trees*/                                                                                                                                                              
void eulerTour(nodeptr p, int numsp, int* array, int* reference, int* pos, int* taxonToEulerIndex);

/*For unrooted Trees there is a special case for the arbitrary root which has degree 3 */
void unrootedEulerTour(nodeptr p, int numsp, int* array, int* reference, int* pos, int* taxonToEulerIndex);

//function for built-in quicksort
int sortIntegers(const void *a, const void *b);
