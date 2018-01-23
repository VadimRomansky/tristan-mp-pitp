#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH -p tornado
#SBATCH -J Tristan
#SBATCH -o pic-%j.out
#SBATCH -e pic-%j.err
#SBATCH -t 240:00:00
mpirun -np 28 /home/ipntsr/romansky/Tristan/tristan-mp2d -i input.shock
