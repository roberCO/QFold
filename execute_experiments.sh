#!/bin/bash

# check the input parameters ($1 path of proteins | $2 initial number of bits | $3 final number of bits)
if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] || [ "$4" == "" ]; then
    echo "Any parameter is missing"
    echo "Parameter 1 (protein name): " $1
    echo "Parameter 2 (amminoacids): " $2 
    echo "Parameter 3 (initial number bits): " $3
    echo "Parameter 4 (final number bits): " $4

    exit 1
fi

for index in $(seq $3 $4)
do
    python main.py $1 $2 $index
    printf "Precalculated energies file for protein $1 with $index bits created!\n\n"
done