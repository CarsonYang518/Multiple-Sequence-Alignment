#!/bin/bash
#SBATCH --partition=physical
#SBATCH --constraint=physg5 
#SBATCH --time=2:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=kaixunyang_kseqalign

module load gcc/10.1.0
g++ -Wall -fopenmp -o kaixunyang_kseqalign kaixunyang_kseqalign.cpp -O3
cat mseq-big13-example.dat | OMP_PLACES=sockets ./kaixunyang_kseqalign