#!/bin/bash

rm matrix.o matrix.e sorted_matrix.txt


module load gcc openmpi/4.1.4_ft3

#mpicc main.c -o matrix -D DEBUG
mpicc main.c -o matrix


sbatch -J matrix -o matrix.o -e matrix.e -N 1 -n $1 --mem=10GB --time=00:10:00 task.sh matrix $1