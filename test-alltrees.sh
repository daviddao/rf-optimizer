#!/bin/bash
rm *.T1
cd standard-RAxML-master/
make -f Makefile.gcc
cd ..
./standard-RAxML-master/raxmlHPC -f R -m GTRCAT -t ./largetest/largetree.txt -z ./largetest/175830-Trees-from-STbase.txt -n T1
