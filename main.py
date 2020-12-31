import sys
import initializer
import angleCalculator
import psiFour
import utils

import time
import datetime

#Read config file with the QFold configuration variables
config_path = './config/config.json'
tools = utils.Utils(config_path)

args = tools.parse_arguments()

print('\n###################################################################')
print('##                             QFOLD                             ##')
print('##                                                               ##')
print('## Tool that combines AI and QC to solve protein folding problem ##')
print('###################################################################\n')

start_time = time.time()

rotationSteps = 2**(int(args.bits))
if args.id == None: args.id = -1

config_variables = tools.get_config_variables()
angleInitializer = initializer.Initializer(
    psi4_path = config_variables['psi4_path'],
    input_file_energies_psi4 = config_variables['input_filename_energy_psi4'], 
    output_file_energies_psi4 = config_variables['output_filename_energy_psi4'],
    energy_method = config_variables['energy_method'],
    precalculated_energies_path = config_variables['precalculated_energies_path'], 
    model_path = config_variables['model_path'], 
    window_size = config_variables['window_size'], 
    max_aa_length = config_variables['maximum_aminoacid_length'],
    initialization_option = config_variables['methods_initialization'],
    n_threads = config_variables['n_threads_pool'],
    basis = config_variables['basis']
    )

psi = psiFour.PsiFour(
    config_variables['psi4_path'], 
    config_variables['input_filename_energy_psi4'], 
    config_variables['output_filename_energy_psi4'], 
    config_variables['precalculated_energies_path'], 
    config_variables['energy_method'], 
    config_variables['n_threads_pool'],
    config_variables['basis'])

#Check if it existes a precalculated energy file with the same parameters, if not call initializer to calculate it
#The format should be energies[args.protein_name][numberBitsForRotation] ex: energiesGlycylglycine2.json
try:
    f = open(config_variables['precalculated_energies_path']+'delta_energies_'+args.protein_name+'_'+str(args.bits)+'_'+args.initialization+'.json')
    f.close()
except IOError:
    print('<!> Info: No precalculated energies file found => Calculating energies\n')
    angleInitializer.calculate_delta_energies(args.protein_name, args.bits, args.initialization, args.aminoacids, args.id)

#Create an empty list of enery list
#HARDCODED for proteins with only two args.aminoacids
#TODO modify to any number of args.aminoacids (it should a list of list, each position of the list contains a list of phi and psi values of this list position)
[deltas_dict, psi4_min_energy, initial_min_energy, index_min_energy, inizialitation_stats] = psi.readEnergyJson(args.protein_name, args.bits, args.initialization)

print('## 3D STRUCTURE CALCULATOR FOR', args.protein_name,'with', args.bits,'bits and', args.initialization,'initialization##\n')

angleCalculator = angleCalculator.AngleCalculator(tools, len(args.aminoacids))

q_accumulated_tts = []
c_accumulated_tts = []
x_axis = []

min_q_tts = {'step': 0, 'value': -1}
min_c_tts = {'step': 0, 'value': -1}

results = []
for step in range(config_variables['initial_step'], config_variables['final_step']):

    # execute for option 0 (quantum) and option 1 (classical)
    for option in ['quantum', 'classical']:

        # calculate the probability matrix of the optimization algorithms
        probabilities_matrix = angleCalculator.calculate3DStructure(deltas_dict, step, option)

        p_t = 0
        # if the index of min energy calculated by psi 4 is in the results of metropolis, p_t is extracted
        # else, the p_t is set to a very small value close to 0 (not 0 to avoid inf values)
        if index_min_energy in probabilities_matrix.keys():
            p_t = probabilities_matrix[index_min_energy]
        else:
            p_t = 0

        result = 0
        # Result is the calculated TTS
        if p_t >= 1:
            result = 1
        elif p_t == 0:
            result = 99999
        else:
            result = tools.calculateTTS(config_variables['precision_solution'], step, p_t)

        
        if option == 'quantum':
            q_accumulated_tts.append(result)
            
            if result < min_q_tts['value'] or min_q_tts['value'] == -1:
        
                min_q_tts['value'] = result
                min_q_tts['step'] = step

        else: 
            c_accumulated_tts.append(result)

            if result < min_c_tts['value'] or min_c_tts['value'] == -1:
        
                min_c_tts['value'] = result
                min_c_tts['step'] = step

    tools.plot_tts(q_accumulated_tts, c_accumulated_tts, args.protein_name, args.aminoacids, args.bits, args.initialization, config_variables['initial_step'])

    final_stats = {'q': min_q_tts, 'c': min_c_tts}

    tools.write_tts(
        config_variables['initial_step'], 
        config_variables['final_step'], 
        q_accumulated_tts, 
        c_accumulated_tts, 
        args.protein_name,
        args.aminoacids,
        args.bits, 
        args.initialization,
        inizialitation_stats,
        final_stats)

# MERGE RESULTS: if the results generated are comparable with similar results generated previously, it generates the shared plot
# For example, if this execution generates results for minifold 4 bits rotation GG and there are results for random 4 bits GG
# it combines the results into only one plot

results = {}
alternative_results_found = False
for alternative_method in config_variables['methods_initialization']:

    if alternative_method != args.initialization:

        try:
            f = open(config_variables['path_tts_plot']+'tts_results_'+args.protein_name+'_'+str(args.bits)+'_'+alternative_method+'.json')
            results[alternative_method] = tools.read_results_file(config_variables['path_tts_plot']+'tts_results_'+args.protein_name+'_'+str(args.bits)+'_'+alternative_method+'.json')
            alternative_results_found = True
            f.close()
        except IOError:
            print('<!> Info: No results for method', alternative_method,'found\n')

if alternative_results_found:

    results[args.initialization] = tools.read_results_file(config_variables['path_tts_plot']+'tts_results_'+args.protein_name+'_'+str(args.bits)+'_'+args.initialization+'.json') 
    tools.generate_combined_results_plot(results, args.protein_name, args.bits)

execution_time = time.time() - start_time

print('\n\n********************************************************')
print('**                       RESULTS                      **')
print('********************************************************')
print('**                                                    **')
print('** Quantum Metropolis   => Min TTS:', '{:.10f}'.format(min_q_tts['value']), 'at step:', min_q_tts['step'], ' **')
print('** Classical Metropolis => Min TTS:', '{:.10f}'.format(min_c_tts['value']), 'at step:', min_c_tts['step'], ' **')
print('**                                                    **')
print('** -------------------------------------------------- **')
print('**                                                    **')
print('** Execution time     =>', str(datetime.timedelta(seconds=execution_time)) ,' in hh:mm:ss  **')
print('********************************************************\n\n')