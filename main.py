import sys
from threading import Thread
import initializer
import angleCalculator
import psiFour
import utils

if(len(sys.argv) != 5):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py proteinName aminoacids_chain numberBitsForRotations method_rotations_generation")
    print ("<!> Example: python main.py Glycylglycine GG 6 random (6 bits for rotations are 64 steps)")
    sys.exit(0)

print('\n###################################################################')
print('##                             QFOLD                             ##')
print('##                                                               ##')
print('## Tool that combines AI and QC to solve protein folding problem ##')
print('###################################################################\n')

proteinName = sys.argv[1].lower()
aminoacids = sys.argv[2]
numberBitsRotation = int(sys.argv[3])
method_rotations_generation = sys.argv[4]
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
    config_variables['basis']
    )

psi = psiFour.PsiFour(config_variables['psi4_path'], config_variables['input_filename_energy_psi4'], config_variables['output_filename_energy_psi4'], config_variables['precalculated_energies_path'], config_variables['energy_method'], config_variables['basis'])

# Number of results is the number of steps multiplied by 2 (one for quantum and other for classical)
total_number_resutls = (config_variables['final_step'] - config_variables['initial_step'])*2
results = results = [{} for x in range(total_number_resutls)]
def angle_calculator_thread(thread_index, option, deltas, step, beta_max, index_min_energy):

        probabilities_matrix = angleCalculator.calculate3DStructure(deltas_dict, step, config_variables['beta_max'], option)

        if option == 0:
            print('\nQuantum probabilities\n')
        else:
            print('\nClassical probabilities\n')

        for key in probabilities_matrix.keys():
            print(key, round(probabilities_matrix[key], 6))

        p_t = 0

        # if the index of min energy calculated by psi 4 is in the results of metropolis, p_t is extracted
        # else, the p_t is set to a very small value close to 0 (not 0 to avoid inf values)
        if index_min_energy in probabilities_matrix.keys():
            p_t = probabilities_matrix[index_min_energy]
        else:
            p_t = 0


        # Result is the calculated TTS

        if p_t >= 1:
            results[thread_index] = 1

        elif p_t == 0:
            results [thread_index] = 9999

        else:
            results[thread_index] = tools.calculateTTS(config_variables['precision_solution'], step, p_t)

#Check if it existes a precalculated energy file with the same parameters, if not call initializer to calculate it
#The format should be energies[proteinName][numberBitsForRotation] ex: energiesGlycylglycine2.json
try:
    f = open(config_variables['precalculated_energies_path']+'delta_energies_'+proteinName+'_'+str(numberBitsRotation)+'_'+method_rotations_generation+'.json')
    f.close()
except IOError:
    print('<!> Info: No precalculated energies file found => Calculating energies\n')
    angleInitializer.calculate_delta_energies(proteinName, numberBitsRotation, method_rotations_generation, aminoacids)

#Create an empty list of enery list
#HARDCODED for proteins with only two aminoacids
#TODO modify to any number of aminoacids (it should a list of list, each position of the list contains a list of phi and psi values of this list position)
[deltas_dict, psi4_min_energy, initial_min_energy, index_min_energy] = psi.readEnergyJson(proteinName, numberBitsRotation, method_rotations_generation)

print('## 3D STRUCTURE CALCULATOR ##\n')

angleCalculator = angleCalculator.AngleCalculator(
    numberBitsRotation, 
    config_variables['ancilla_bits'], 
    config_variables['scaling_factor'], 
    config_variables['number_iterations'],
    len(aminoacids)
    )

q_accumulated_tts = []
c_accumulated_tts = []
x_axis = []

min_q_tts = {'step': 0, 'value': -1}
min_c_tts = {'step': 0, 'value': -1}

threads = []
thread_index = 0
index_to_get_results = []
for step in range(config_variables['initial_step'], config_variables['final_step']):

    print('\nExecuting quantum metropolis with', step, 'steps')

    #Thread for quantum metropolis
    process = Thread(target=angle_calculator_thread, args=[thread_index, 0, deltas_dict, step, config_variables['beta_max'], index_min_energy])
    process.start()
    threads.append(process) 
    index_to_get_results.append(thread_index)
    thread_index += 1


    print('Executing classical metropolis with', step, 'steps\n')

    #Thread for classical metropolis
    process = Thread(target=angle_calculator_thread, args=[thread_index, 1, deltas_dict, step, config_variables['beta_max'], index_min_energy])
    process.start()
    threads.append(process)
    index_to_get_results.append(thread_index)
    thread_index += 1

    if thread_index % config_variables['n_threads_pool'] == 0 or thread_index + config_variables['n_threads_pool'] >= (config_variables['final_step'] - config_variables['initial_step']):

        # It pauses execution until all threads ends
        for process in threads:
            process.join()

        for index in range(index_to_get_results[0], index_to_get_results[-1]+1, 2):

            quantum_TTS = results[index]
            classical_TTS = results[index+1]

            if quantum_TTS < min_q_tts['value'] or min_q_tts['value'] == -1:
                
                min_q_tts['value'] = quantum_TTS
                min_q_tts['step'] = step

            if classical_TTS < min_c_tts['value'] or min_c_tts['value'] == -1:
                
                min_c_tts['value'] = classical_TTS
                min_c_tts['step'] = step

            q_accumulated_tts.append(quantum_TTS)
            c_accumulated_tts.append(classical_TTS)
            x_axis.append(step)

            tools.plot_tts(x_axis, q_accumulated_tts, c_accumulated_tts, proteinName, numberBitsRotation, method_rotations_generation)

        index_to_get_results = []

# Difference between the minimum energy of initializer minus the minimum energy of psi4
min_energy_difference = (1 - (initial_min_energy - psi4_min_energy)) *100
delta_mean = tools.calculate_delta_mean(deltas_dict)
std_dev_deltas = tools.calculate_std_dev_deltas(deltas_dict)

tools.write_tts(config_variables['initial_step'], config_variables['final_step'], q_accumulated_tts, c_accumulated_tts, proteinName, numberBitsRotation, method_rotations_generation)

# Compare the difference between the minimum energy of initializer minus the minimum energy of psi4 with the mean of energy deltas
precision_vs_delta_mean = tools.calculate_diff_vs_mean_diffs(min_energy_difference, delta_mean)

# MERGE RESULTS: if the results generated are comparable with similar results generated previously, it generates the shared plot
# For example, if this execution generates results for minifold 4 bits rotation GG and there are results for random 4 bits GG
# it combines the results into only one plot

results = {}
alternative_results_found = False
for alternative_method in config_variables['methods_initialization']:

    if alternative_method != method_rotations_generation:

        try:
            f = open(config_variables['path_tts_plot']+'tts_results_'+proteinName+'_'+str(numberBitsRotation)+'_'+alternative_method+'.json')
            results[alternative_method] = tools.read_results_file(config_variables['path_tts_plot']+'tts_results_'+proteinName+'_'+str(numberBitsRotation)+'_'+alternative_method+'.json')
            alternative_results_found = True
            f.close()
        except IOError:
            print('<!> Info: No results for method', alternative_method,'found\n')

if alternative_results_found:

    results[method_rotations_generation] = tools.read_results_file(config_variables['path_tts_plot']+'tts_results_'+proteinName+'_'+str(numberBitsRotation)+'_'+method_rotations_generation+'.json') 
    tools.generate_combined_results_plot(results, proteinName, numberBitsRotation)


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