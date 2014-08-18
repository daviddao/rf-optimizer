#!/bin/bash
rm standard-RAxML-master/*.o
rm *.T1
./standard-RAxML-master/make -f Makefile.gcc
./standard-RAxML-master/raxmlHPC -f R -m GTRCAT -t ./test/largetree -z ./test/reference -n T1
