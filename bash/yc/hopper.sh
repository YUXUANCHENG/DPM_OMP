#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
# maindir=$cellsdir/main

# compile into binary using packing.h
workdir=$(pwd)
binf=$(pwd)/jamming.o
jobnumber=20
factor1=5
factor2=2
# mainf=$maindir/jamming/cellJamming.cpp

# run compiler
rm -f $binf
g++ -O3 --std=c++11 -fopenmp -I $srcdir $srcdir/*.cpp -o $binf 

taskf=$workdir/task.txt
rm -f $taskf

let range2=$jobnumber*$factor2-1
let range1=$jobnumber*$factor1-1
for index_i in `seq 0 $range1`; do
    for index_j in `seq 0 $range2`; do
        current=$workdir/"$index_i"_"$index_j"/
        runString="mkdir -p $current;cd $current;$binf $index_i $index_j;"
        echo "$runString" >> $taskf
    done
done

slurmf=$workdir/slurm.sh
rm -f $slurmf

partition=pi_ohern
job_name=DPM
let total_job=$jobnumber*$jobnumber*$factor1*$factor2

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --mem-per-cpu=512 >> $slurmf
echo \#SBATCH --array=1-$total_job >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#sbatch --mail-user=yuxuan.cheng@yale.edu >> $slurmf
echo \#sbatch --mail-type=ALL >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

sbatch -t 14-00:00:00 $slurmf