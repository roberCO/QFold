import json
from os import listdir
from os.path import isfile, join
import re


quantum_folder = "/Users/pabloantoniomorenocasares/Documents/GitHub/QFold/results/"
classical_folder = "/Users/pabloantoniomorenocasares/Documents/GitHub/QFold/results_1000/"
combined_folder = "/Users/pabloantoniomorenocasares/Documents/GitHub/QFold/results_combined/"

input_files = [f for f in listdir(quantum_folder) if isfile(join(quantum_folder, f))]

for input_name in input_files:

    if input_name[:3] == 'tts' and input_name.split('.')[1] == 'json':

        quantum_file = input_name
        params = re.split('\.|\_',input_name)
        params = params[:-1] #remove the file extension
        beta_q = int(params[-1])
        beta_c = 1000

        classical_file = '_'.join(params[:-1]+ ['1000']) + '.json'
        combined_file = '_'.join(params) + '_1000.json'

        with open(quantum_folder + quantum_file) as json_file:
            quantum_data = json.load(json_file)

        try:
            with open(classical_folder + classical_file) as json_file:
                classical_data = json.load(json_file)
        except:
            print("Classical file not found", classical_file)
            continue

        combined_results = {}
        combined_results['initial_step'] = quantum_data['initial_step']
        combined_results['final_step'] = quantum_data['final_step']
        combined_results['quantum_tts'] = quantum_data['quantum_tts']
        combined_results['classical_tts'] = classical_data['classical_tts']
        combined_results['initialization_stats'] = quantum_data['initialization_stats']
        combined_results['final_stats'] = {}
        combined_results['final_stats']['q'] = quantum_data['final_stats']['q']
        combined_results['final_stats']['c'] = classical_data['final_stats']['c']
        with open(combined_folder + combined_file, 'w') as fp:
            json.dump(combined_results, fp)