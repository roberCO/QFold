import sys
import minifoldTrainer
import sys 
sys.path.insert(1, '/home/roberto/Escritorio/qfold/QFold/')
import utils

if len(sys.argv) == 2 and sys.argv[1] == 'list_trainers':

    print('Available trainer options:')
    print('0 -> Minifold (requires a CASP file) => Comand: python initialAngleTrainer.py trainer 0')

elif len(sys.argv) == 3 and sys.argv[1] == 'trainer':

    option = int(sys.argv[2])

    #Read config file with the QFold configuration variables
    config_path = './config/config.json'

    tools = utils.Utils(config_path)
    config_variables = tools.get_config_variables()

    if option == 0:

        print('############################################################')
        print('##                    MINIFOLD TRAINER                    ##')
        print('##                                                        ##')
        print('## It generates a knowledge model from a CASP input file  ##')
        print('############################################################\n')

        pathInputFile = input("Path input file: ")

        minifold = minifoldTrainer.MinifoldTrainer(
            pathInputFile,
            config_variables['model_path'],
            config_variables['maximum_aminoacid_length'],
            config_variables['window_size'],
            config_variables['epochs'],
            config_variables['batch_size'],
            )
        minifold.train()

else: 

    print("<*> ERROR: Wrong number of parameters - Usage options:")
    print(" ■ Request available trainers: python initialAngleTrainer.py list_trainers") 
    print(" ■ Execute a trainer: python initialAngleTrainer.py trainer 0") 

    sys.exit(0)

