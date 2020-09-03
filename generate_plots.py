import utils

import json
import numpy as np

#Read config file with the QFold configuration variables
config_path = './config/config.json'

tools = utils.Utils(config_path)
config_variables = tools.get_config_variables()

# list elements to read
input_files = [
    'glycylglycylglycine_GGG_1_minifold_1_-1'
]

for input_name in input_files:

    data = {}
    # read data
    with open(config_variables['path_tts_plot']+ 'tts_results_'+input_name+'.json') as json_file:
                data[input_name] = json.load(json_file)


# generate plot of minifold vs random inizialization mode
precisions = {}
min_tts = {}

for protein_key in data.keys():

    aas = protein_key.split('_')[1]
    bits = protein_key.split('_')[2]
    init_mode = protein_key.split('_')[3]
    phi_prec = data[protein_key]['initialization_stats']['phis_precision']
    psi_prec = data[protein_key]['initialization_stats']['psis_precision']
    precisions[aas+'_'+bits+'_'+init_mode] = np.mean(np.mean([float(p) for p in phi_prec]) + np.mean([float(p) for p in psi_prec]))
    min_tts[aas+'_'+bits+'_'+init_mode] = min(data[protein_key]['final_stats']['q']['value'], data[protein_key]['final_stats']['c']['value'])


# generate plot of the evolution of tts with different steps comparing classical vs quantum

# generate plot of the comparison between quantum and classical difference of tts  