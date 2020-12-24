#!/bin/bash

### PBS - Pro Job Scheduler

## Job Name
#PBS -N distributed_sp

## PBS output file, errorfile
#PBS -o pbs.log
#PBS -e errorfile.err

## Job Time Limit (including being in the queue) HH:MM:SS
#PBS -l walltime=24:00:00

## Job TIme Limit (CPU) HH:MM:SS
#PBS -l cput=00:10:00

## Number of nodes and processor per node to request
## Requesting 6 nodes and 1 processor per node
#PBS -l select=20:ncpus=1

## Get the job id
job_id=`echo $PBS_JOBID | cut -f 1 -d .`

## Get number of processors
np=$(wc -l < $PBS_NODEFILE)

## Run the job
I_MPI_HYDRA_TOPOLIB=ipl mpiexec.hydra -np $np -genv OMP_NUM_THREADS=1 -genv I_MPI_PIN=1 -genv I_MPI_FABRICS=shm:ofi -hostfile $PBS_NODEFILE $PBS_O_WORKDIR/main > $PBS_O_WORKDIR/output.txt

#echo "$np" > out.txt

#cd "$PBS_O_WORKDIR"

#echo "$PBS_O_WORKDIR" > out.txt

#tpdir=`echo $PBS_JOBID | cut -f 1 -d .`

#tempdir=$HOME/scratch/job$tpdir

#mkdir -p $tempdir
#cd $tempdir
#cp -R $PBS_O_WORKDIR/* .
#mpicc inputfile.c
#mpiexec.hydra -np 6 -genv OMP_NUM_THREADS=1 -genv I_MPI_PIN=1 -genv I_MPI_FABRICS=shm:ofi -hostfile $PBS_NODEFILE ./a.out > output.txt
#rm a.out
#mv ../job$tpdir $PBS_O_WORKDIR/.