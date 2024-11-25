#!/bin/bash

module load gcc openmpi/4.1.4_ft3
for p in $(seq 2 2 64)
do
    nombreTrabajo="matrix$p"
    sbatch -J $nombreTrabajo -o matrix.o -e matrix.e -N 2 -n $p --mem=64GB --time=01:59:00 taskVariasEjecuciones.sh matrix $p 
done