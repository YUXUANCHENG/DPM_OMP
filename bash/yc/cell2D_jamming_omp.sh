#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
# maindir=$cellsdir/main

# compile into binary using packing.h
binf=$(pwd)/jamming.o
# mainf=$maindir/jamming/cellJamming.cpp

# run compiler
rm -f $binf
g++ --std=c++11 -fopenmp -I $srcdir $srcdir/*.cpp -o $binf 

# $binf



