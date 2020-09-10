#!/bin/bash

python main.py glycylglycine GG 3 minifold
python main.py glycylglycine GG 3 random
python main.py glycylglycine GG 4 minifold
python main.py glycylglycine GG 4 random


python main.py alanylcysteine AC 3 minifold
python main.py alanylcysteine AC 3 random
python main.py alanylcysteine AC 4 minifold
python main.py alanylcysteine AC 4 random

mv config/config.json config/config_temp.json
jq -r '.final_step |= 10' config/config_temp.json > config/config.json

python main.py serylalanine SA 3 minifold
python main.py serylalanine SA 3 random
python main.py serylalanine SA 4 minifold
python main.py serylalanine SA 4 random

python main.py valylglycine VG 3 minifold
python main.py valylglycine VG 3 random
python main.py valylglycine VG 4 minifold
python main.py valylglycine VG 4 random

python main.py glycylcysteine GC 3 minifold
python main.py glycylcysteine GC 3 random
python main.py glycylcysteine GC 4 minifold
python main.py glycylcysteine GC 4 random

python main.py glycylglycylglycine GGG 3 minifold
python main.py glycylglycylglycine GGG 3 random
python main.py glycylglycylglycine GGG 4 minifold
python main.py glycylglycylglycine GGG 4 random
