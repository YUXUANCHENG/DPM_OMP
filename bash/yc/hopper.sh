#!/bin/bash

# directories with code
cellsdir=~/DPM_OMP_Hopper/DPM_OMP
srcdir=$cellsdir/src
# maindir=$cellsdir/main

# compile into binary using packing.h
workdir=$(pwd)
binf=$(pwd)/jamming.o
jobnumber=1
# factor1=10
# factor2=40


factor1=19
factor1s=1
factor2=20

# factor1=19
# factor1s=18
# factor2=100

# mainf=$maindir/jamming/cellJamming.cpp

# run compiler
rm -f $binf
ml GCC/10
g++ -O3 -std=c++17 -fopenmp -I $srcdir $srcdir/*.cpp -o $binf 

taskf=$workdir/task.txt
rm -f $taskf

let range2=$jobnumber*$factor2-1
let range1=$jobnumber*$factor1-1
let range1s=$factor1s-1
# let range1=$factor1
for index_i in `seq $range1s $range1`; do
    for index_j in `seq 0 $range2`; do
        current=$workdir/"$index_i"_"$index_j"/
        runString="mkdir -p $current;cd $current;$binf $index_i $index_j;"
        echo "$runString" >> $taskf
    done
done

slurmf=$workdir/slurm.sh
rm -f $slurmf

# partition=scavenge
# partition=pi_ohern
# partition=general
partition=day
job_name=hopper
let range1t=$factor1-$factor1s+1
let total_job=$jobnumber*$jobnumber*$range1t*$factor2

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
# echo \#SBATCH --cpus-per-task=4 >> $slurmf
# echo \#SBATCH --cpus-per-task=6 >> $slurmf
echo \#SBATCH --mem-per-cpu=512 >> $slurmf
echo \#SBATCH --array=1-$total_job >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#sbatch --mail-user=yuxuan.cheng@yale.edu >> $slurmf
echo \#sbatch --mail-type=ALL >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# sbatch -t 07:00:00 $slurmf
sbatch -t 1-00:00:00 $slurmf
# sbatch -t 3-00:00:00 $slurmf