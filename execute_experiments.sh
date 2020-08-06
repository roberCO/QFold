#!/bin/bash

# check the input parameters ($1 path of proteins | $2 initial number of bits | $3 final number of bits)
if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ]; then
    echo "Any parameter is missing"
    echo "Parameter 1 (proteins path): " $1
    echo "Parameter 2 (initial number bits): " $2
    echo "Parameter 3 (final number bits): " $3

    exit 1
fi

NUMBER_CORES=$(grep -c ^processor /proc/cpuinfo)

while read p; do
    protein="$(cut -d' ' -f1 <<<"$p")"
    aa="$(cut -d' ' -f2 <<<"$p")"
    for index in $(seq $2 $3)
    do
    
        $(python main.py $protein $aa $index)
        echo "Precalculated energies file for protein $protein with $index bits created!\n"
    done
done <$1