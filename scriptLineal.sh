#!/bin/bash
rm matrix.o matrix.e sorted_matrix.txt


module load gcc openmpi/4.1.4_ft3

mpicc main.c -o matrix

sbatch -J matrix -o matrix.o -e matrix.e -N 1 -n 2 --mem=32GB --time=01:59:00 taskVariasEjecuciones.sh matrix 2
