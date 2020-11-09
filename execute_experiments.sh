#!/bin/bash

# check the input parameters ($1 path of proteins | $2 initial number of bits | $3 final number of bits)
if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ]; then
    echo "Any parameter is missing"
    echo "Parameter 1 (proteins path): " $1
    echo "Parameter 2 (initial number bits): " $2
    echo "Parameter 3 (final number bits): " $3

    exit 1
fi

# read all lines of the file
lines=()
while read line; do
    lines+=("$line")
done <$1

for index in $(seq $2 $3); do

    for beta in 100 1000 10000; do

        mv config/config.json config/config_temp.json
        jq -r '.beta |= $beta' config/config_temp.json > config/config.json
        rm config/config_temp.json
    
        for line in "${lines[@]}"; do

            proteinline="$(tr -s ' ' <<< "$line")"
            protein="$(cut -d' ' -f1 <<<"$proteinline")"
            aa="$(cut -d' ' -f2 <<<"$proteinline")"
            id="$(cut -d' ' -f3 <<<"$proteinline")"
            id=${id/"("/""}
            id=${id/")"/""}

            if [[ $id == *"#"* ]]; then

                for init in minifold random; do
                    python main.py $protein $aa $index $init
                done
            else
                for init in minifold random; do
                    python main.py $protein $aa $index $init -i $id
                done
            fi
        done
    done
done
