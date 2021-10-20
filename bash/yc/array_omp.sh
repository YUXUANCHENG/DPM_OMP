#!/bin/bash

# directories with code
cellsdir=~/DPM_OMP
srcdir=$cellsdir/src
# maindir=$cellsdir/main

# compile into binary using packing.h
workdir=$(pwd)
binf=$(pwd)/jamming.o
jobnumber=10;
# mainf=$maindir/jamming/cellJamming.cpp

# run compiler
rm -f $binf
g++ -O3 --std=c++11 -fopenmp -I $srcdir $srcdir/*.cpp -o $binf 

taskf=$workdir/task.txt
rm -f $taskf

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
rm -f $slurmf
partition=pi_ohern
#partition=general
job_name=DPM

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=10 >> $slurmf
echo \#SBATCH --mem-per-cpu=5G >> $slurmf
echo \#SBATCH --array=1-$jobnumber >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

sbatch -t 14-00:00:00 $slurmf



