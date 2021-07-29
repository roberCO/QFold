import os
import time
import datetime
import initializer
import psiFour
import utils
import qms

print('\n###################################################################')
print('##                             QFOLD                             ##')
print('##                                                               ##')
print('## Tool that combines AI and QC to solve protein folding problem ##')
print('###################################################################\n')

start_time = time.time()

#Read config file with the QFold configuration variables
config_path = './config/config.json'

tools = utils.Utils(config_path)
angleInitializer = initializer.Initializer(tools)
psi = psiFour.PsiFour(tools)

args = tools.parse_arguments()

#Check if it existes a precalculated energy file with the same parameters, if not call initializer to calculate it
#The format should be energies[args.protein_name][numberBitsForRotation][initialization] ex: energies_glycylglycine_2_minifold.json
try:
    f = open(tools.config_variables['precalculated_energies_path']+'energies_'+args.protein_name+'_'+str(args.bits)+'_'+args.initialization+'.json')
    f.close()
except IOError:
    print('<!> Info: No precalculated energies file found => Calculating energies\n')
    angleInitializer.calculate_energies(args.protein_name, args.bits, args.initialization, args.aminoacids, args.id)

print('\n## 3D STRUCTURE CALCULATOR FOR', args.protein_name,'with', args.bits,'bits and', args.initialization,'initialization##\n')

file_energies_path_qms = tools.config_variables['precalculated_energies_path'] + 'file_energies_qms_' + args.protein_name + '_' + str(args.bits) + '_' + args.initialization + '.json'

[psi4_min_energy, index_min_energy, initialization_stats] = tools.write_qms_energies_file(args.protein_name, args.bits, args.initialization, file_energies_path_qms)

# call QMS module to calculate results
[min_q_tts, min_c_tts] = qms.executeMetropolis(path=file_energies_path_qms, isCircular=True)

# delete temporal energies file for qms
if os.path.exists(file_energies_path_qms):
    os.remove(file_energies_path_qms)
else:
    print('<*> ERROR: Temporal qms file ' + file_energies_path_qms + 'does not exist! Probably QFold is not working correctly')

execution_time = time.time() - start_time

print('\n********************************************************')
print('**       RESULTS for ', args.protein_name,'with', args.bits,'bits       **')
print('********************************************************')
print('**                                                    **')

if args.mode == 'simulation' or args.mode == 'experiment':
    print('** Quantum Metropolis   => Min TTS:', '{:.4f}'.format(min_q_tts['value']), 'at step:', min_q_tts['step'], ' **')
    print('** Classical Metropolis => Min TTS:', '{:.4f}'.format(min_c_tts['value']), 'at step:', min_c_tts['step'], ' **')

elif args.mode == 'real':

    print('** Quantum Metropolis   => Confidence:', '{:.4f}'.format(min_q_tts['value']*100), '% **')
    print('** Classical Metropolis => Confidence:', '{:.4f}'.format(min_c_tts['value']*100), '% **')
    print('**                                                    **')
    print('**      Quantum success min energy:', min_q_tts['success'],'      **')
    print('**      Classical success min energy:', min_c_tts['success'],'      **')

print('**                                                    **')
print('** -------------------------------------------------- **')
print('**                                                    **')
print('** Execution time     =>', str(datetime.timedelta(seconds=execution_time)) ,' in hh:mm:ss  **')
print('********************************************************\n\n')