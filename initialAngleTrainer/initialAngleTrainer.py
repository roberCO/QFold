import sys
import minifoldTrainer

if len(sys.argv) == 2 and sys.argv[1] == 'list_trainers':

    print('Available trainer options:')
    print('0 -> Minifold (requires a CASP file) => Comand: python initialAngleTrainer.py trainer 0')

elif len(sys.argv) == 3 and sys.argv[1] == 'trainer':

    option = int(sys.argv[2])

    if option == 0:

        pathInputFile = input("Path input file: ")

        minifold = minifoldTrainer.MinifoldTrainer(pathInputFile) 
        minifold.train()

else: 

    print("<*> ERROR: Wrong number of parameters - Usage options:")
    print(" ■ Request available trainers: python initialAngleTrainer.py list_trainers") 
    print(" ■ Execute a trainer: python initialAngleTrainer.py trainer 0") 

    sys.exit(0)

