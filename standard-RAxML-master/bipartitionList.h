void checkTreeNumber(int numberOfTrees, char *fileName);
int findHash(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, hashNumberType position);
void insertHashPlausibility(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, hashNumberType position);
hashNumberType oat_hash(unsigned char *p, int len);
void newviewBipartitionsMultifurcating(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength);
void freeBitVectors(unsigned int **v, int n);
