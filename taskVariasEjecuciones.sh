#!/bin/bash
for i in $(seq 5051 50 6551)
do
    FILENAME="matrices/matrix${i}.txt"
    srun $1 $2 $FILENAME
done
