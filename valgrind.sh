#!/bin/bash
rm *.T1
cd standard-RAxML-master/
make -f Makefile.gcc
rm *.o
cd ..
valgrind --track-origins=yes ./standard-RAxML-master/raxmlHPC -f R -m GTRCAT -t ./test/largetree -z ./test/reference -n T1
