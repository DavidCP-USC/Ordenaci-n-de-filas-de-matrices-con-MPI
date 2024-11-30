#!/bin/bash

for i in $(seq 5051 50 10001)
do
    stdbuf -oL echo "Creando matriz $i x $i"
    ./createMatrix $i
done