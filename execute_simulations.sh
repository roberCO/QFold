#!/bin/bash

# check the input parameters ($1 path of proteins | $2 initial number of bits | $3 final number of bits)
if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] || [ "$4" == "" ] ; then
    echo "Any parameter is missing"
    echo "Parameter 1 (proteins name): " $1
    echo "Parameter 2 (aminoacids): " $2
    echo "Parameter 3 (initial number bits): " $3
    echo "Parameter 4 (final number bits): " $4

    exit 1
fi

for index in $(seq $3 $4); do
    
    python main.py $1 $2 $index minifold simulation
    python main.py $1 $2 $index random simulation

done