#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
# maindir=$cellsdir/main

# directory for all output for cell simulations
outputdir=~/project/cells

# directory for simulations specific to jamming
simtypedir=$outputdir/dpjam

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin

# compile into binary using packing.h
binf=bin/jamming.o
# mainf=$maindir/jamming/cellJamming.cpp
# echo Running $numSeeds jamming sims of $NCELLS cells with $NV verts, with dispersity $sizeDisp , calA0 = $calA0 , and bending energy kb = $kb

# run compiler
rm -f $binf
g++ --std=c++11 -I $srcdir $mainf $srcdir/*.cpp -o $binf 
echo compiling with : g++ --std=c++11 -I $srcdir $srcdir/*.cpp -o $binf 

./$binf



