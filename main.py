import sys
import initializer
import angleCalculator
import psiFour
import utils

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)

#HARDCODED
aminoacids = 'GG'

proteinName = sys.argv[1].lower()
numberBitsRotation = int(sys.argv[2])
rotationSteps = pow(2, int(numberBitsRotation))

#Read config file with the QFold configuration variables
config_path = './config/config.json'

tools = utils.Utils(config_path)
config_variables = tools.get_config_variables()
angleInitializer = initializer.Initializer(
    config_variables['psi4_path'],
    config_variables['input_filename_energy_psi4'], 
    config_variables['output_filename_energy_psi4'],
    config_variables['energy_method'],
    config_variables['precalculated_energies_path'], 
    config_variables['model_path'], 
    config_variables['window_size'], 
    config_variables['maximum_aminoacid_length'],
    config_variables['initialization_option']
    )

angleCalculator = angleCalculator.AngleCalculator(numberBitsRotation, config_variables['ancilla_bits'], config_variables['scaling_factor'], config_variables['number_iterations'])
psi = psiFour.PsiFour(config_variables['psi4_path'], config_variables['input_filename_energy_psi4'], config_variables['output_filename_energy_psi4'], config_variables['precalculated_energies_path'], config_variables['energy_method'])

#Check if it existes a precalculated energy file with the same parameters, if not call initializer to calculate it
#The format should be energies[proteinName][numberBitsForRotation] ex: energiesGlycylglycine2.json
try:
    f = open(config_variables['precalculated_energies_path']+'energies_'+proteinName+'_'+str(numberBitsRotation)+'.json')
    f.close()
except IOError:
    print('<!> Info: No precalculated energies file found => Calculating energies')
    angleInitializer.calculateEnergies(proteinName, numberBitsRotation, aminoacids)

#Create an empty list of enery list
#HARDCODED for proteins with only two aminoacids
#TODO modify to any number of aminoacids (it should a list of list, each position of the list contains a list of phi and psi values of this list position)
[energyList, phi_angle_psi4, psi_angle_psi4] = psi.readEnergyJson(proteinName, numberBitsRotation)

quantum_p_t = angleCalculator.calculate3DStructure(energyList, config_variables['steps'], config_variables['beta_max'], 0)[0][0]
classical_matrix = angleCalculator.calculate3DStructure(energyList, config_variables['steps'], config_variables['beta_max'], 1)
classical_p_t = classical_matrix[0][0]

quantum_TTS = tools.calculateTTS(config_variables['precision_solution'], config_variables['steps'], quantum_p_t)
classical_TTS = tools.calculateTTS(config_variables['precision_solution'], config_variables['steps'], classical_p_t)

print('Quantum Metropolis => TTS:', quantum_TTS)
print('Classical Metropolis => TTS:', classical_TTS)
