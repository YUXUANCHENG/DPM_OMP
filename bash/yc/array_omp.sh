#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
# maindir=$cellsdir/main

# compile into binary using packing.h
workdir=$(pwd)
binf=$(pwd)/jamming.o
jobnumber=10;
# mainf=$maindir/jamming/cellJamming.cpp

# run compiler
rm -f $binf
g++ --std=c++11 -fopenmp -I $srcdir $srcdir/*.cpp -o $binf 

taskf=$workdir/task.txt

let range=$jobnumber-1
for index in `seq 0 $range`; do
#    cd $workdir
   current=$workdir/$index/
 #   mkdir -p $current   
 #   cd $current
 #   touch hello.world
    runString="mkdir -p $current;cd $current;$binf $index;"
    echo "$runString" >> $taskf
done

slurmf=$workdir/slurm.sh
partition=pi_ohern
job_name=DPM

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=10 >> $slurmf
echo \#SBATCH --array=1-$jobnumber >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

sbatch -t 5-00:00:00 $slurmf



