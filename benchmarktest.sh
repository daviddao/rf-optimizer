#!/bin/bash
rm *.T1
cd standard-RAxML-master/
make -f Makefile.gcc
rm *.o
cd ..
./standard-RAxML-master/raxmlHPC -f R -m GTRCAT -t ./test/largetree -z ./test/30_reference -n T1
