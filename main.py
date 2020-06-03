import sys
import initializer
import angleCalculator
import psiFour
import utils

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)

print('\n###################################################################')
print('##                             QFOLD                             ##')
print('##                                                               ##')
print('## Tool that combines AI and QC to solve protein folding problem ##')
print('###################################################################\n')

#HARDCODED
aminoacids = 'GG'

proteinName = sys.argv[1].lower()
numberBitsRotation = int(sys.argv[2])
rotationSteps = 2**(int(numberBitsRotation))

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
    config_variables['initialization_option'],
    config_variables['basis']
    )

psi = psiFour.PsiFour(config_variables['psi4_path'], config_variables['input_filename_energy_psi4'], config_variables['output_filename_energy_psi4'], config_variables['precalculated_energies_path'], config_variables['energy_method'])

#Check if it existes a precalculated energy file with the same parameters, if not call initializer to calculate it
#The format should be energies[proteinName][numberBitsForRotation] ex: energiesGlycylglycine2.json
try:
    f = open(config_variables['precalculated_energies_path']+'energies_'+proteinName+'_'+str(numberBitsRotation)+'.json')
    f.close()
except IOError:
    print('<!> Info: No precalculated energies file found => Calculating energies\n')
    angleInitializer.calculate_delta_energies(proteinName, numberBitsRotation, aminoacids)

#Create an empty list of enery list
#HARDCODED for proteins with only two aminoacids
#TODO modify to any number of aminoacids (it should a list of list, each position of the list contains a list of phi and psi values of this list position)
[deltas_dict, psi4_min_energy, initial_min_energy, phi_position_min_energy, psi_position_min_energy] = psi.readEnergyJson(proteinName, numberBitsRotation)

print('## 3D STRUCTURE CALCULATOR ##\n')

angleCalculator = angleCalculator.AngleCalculator(numberBitsRotation, config_variables['ancilla_bits'], config_variables['scaling_factor'], config_variables['number_iterations'])

q_accumulated_tts = []
c_accumulated_tts = []
x_axis = []

min_q_tts = {'step': 0, 'value': -1}
min_c_tts = {'step': 0, 'value': -1}

for step in range(config_variables['initial_step'], config_variables['final_step']):

    quantum_matrix = angleCalculator.calculate3DStructure(deltas_dict, step, config_variables['beta_max'], 0)
    classical_matrix = angleCalculator.calculate3DStructure(deltas_dict, step, config_variables['beta_max'], 1)

    quantum_p_t = quantum_matrix[phi_position_min_energy][psi_position_min_energy]
    classical_p_t = classical_matrix[phi_position_min_energy][psi_position_min_energy]

    quantum_TTS = tools.calculateTTS(config_variables['precision_solution'], step, quantum_p_t)
    classical_TTS = tools.calculateTTS(config_variables['precision_solution'], step, classical_p_t)

    if quantum_TTS < min_q_tts['value'] or min_q_tts['value'] == -1:
        
        min_q_tts['value'] = quantum_TTS
        min_q_tts['step'] = step

    if classical_TTS < min_c_tts['value'] or min_c_tts['value'] == -1:
        
        min_c_tts['value'] = classical_TTS
        min_c_tts['step'] = step

    q_accumulated_tts.append(quantum_TTS)
    c_accumulated_tts.append(classical_TTS)
    x_axis.append(step)

    tools.plot_tts(x_axis, q_accumulated_tts, c_accumulated_tts, proteinName, numberBitsRotation)

# Difference between the minimum energy of initializer minus the minimum energy of psi4
min_energy_difference = (1 - (initial_min_energy - psi4_min_energy)) *100
delta_mean = tools.calculate_delta_mean(deltas_dict)
std_dev_deltas = tools.calculate_std_dev_deltas(deltas_dict)

# Compare the difference between the minimum energy of initializer minus the minimum energy of psi4 with the mean of energy deltas
precision_vs_delta_mean = tools.calculate_diff_vs_mean_diffs(min_energy_difference, delta_mean)


print('\n\n********************************************************')
print('**                       RESULTS                      **')
print('********************************************************')
print('**                                                    **')
print('** Quantum Metropolis   => Min TTS:', '{:.10f}'.format(min_q_tts['value']), 'at step:', min_q_tts['step'], ' **')
print('** Classical Metropolis => Min TTS:', '{:.10f}'.format(min_c_tts['value']), 'at step:', min_c_tts['step'], ' **')
print('**                                                    **')
print('** -------------------------------------------------- **')
print('**                                                    **')
print('** Precision QFold     =>', min_energy_difference,'%        **')
print('** Precision vs Δ mean =>', precision_vs_delta_mean ,'     **')
print('** Mean Δ              =>', delta_mean, '                  **')
print('** Standard deviation  =>', std_dev_deltas, '        **')
print('**                                                    **')
print('********************************************************\n\n')