import numpy as np
from prettytable import PrettyTable

def calculate_stats(data):

    mean_precisions = calculate_mean_precision(data)
    tts_tables = calculate_best_tts_parameters(data)

    return [mean_precisions, tts_tables]

def calculate_mean_precision(data):

    minifold_tts = []
    minifold_precision = []

    random_tts = []
    random_precision = []

    quantum_tts = []
    classical_tts = []

    for key in data.keys():

        quantum_tts.append(data[key]['min_tts_q'])
        classical_tts.append(data[key]['min_tts_c'])

        if 'minifold' in key:

            minifold_tts.append(data[key]['min_tts'])
            minifold_precision.append(data[key]['precision'])

        elif 'random' in key:

            random_tts.append(data[key]['min_tts'])
            random_precision.append(data[key]['precision'])

    mean_stats = {}
    mean_stats['minifold_precision'] = np.mean(minifold_precision)
    mean_stats['minifold_tts'] = np.mean(minifold_tts)
    mean_stats['random_precision'] = np.mean(random_precision)
    mean_stats['random_tts'] = np.mean(random_tts)
    mean_stats['quantum_tts'] = np.mean(quantum_tts)
    mean_stats['classical_tts'] = np.mean(classical_tts)

    return mean_stats
    
def calculate_best_tts_parameters(data):

    structured_data = organize_data(data)
    
    data_tables = []

    for protein_key in structured_data.keys():
        for bits_key in structured_data[protein_key].keys():
            for method_key in structured_data[protein_key][bits_key].keys():

                field_names = []
                q_tts_data = {}
                c_tts_data = {}

                for parameters_key in structured_data[protein_key][bits_key][method_key].keys():

                    field_names.append(parameters_key)
                    q_tts_data[parameters_key] = structured_data[protein_key][bits_key][method_key][parameters_key]['min_tts_q']
                    c_tts_data[parameters_key] = structured_data[protein_key][bits_key][method_key][parameters_key]['min_tts_c']   

                data_tables.append(build_table(protein_key, bits_key, method_key, field_names, q_tts_data, 'quantum'))
                data_tables.append(build_table(protein_key, bits_key, method_key, field_names, c_tts_data, 'classical'))

    return data_tables

def build_table(protein_key, bits_key, method_key, field_names, tts_data, algorithm):

    table = PrettyTable()

    # sort tts data

    tts_data = {k: v for k, v in sorted(tts_data.items(), key=lambda item: item[1])}

    table.field_names = tts_data
    table.add_row([tts_data[tts_key] for tts_key in tts_data.keys()])

    table_title = protein_key + ' - ' + bits_key + '   ' + method_key + '  ' + algorithm
    table = table.get_string(title=table_title)

    return table

def organize_data(data):

    structured_data = {}

    for protein_key in data.keys():
        
        protein_name = protein_key.split('_')[0]
        bits_name = protein_key.split('_')[1]
        method_name = protein_key.split('_')[2]

        beta_type = 'undefined'
        if protein_key.split('_')[4] == '0':
            beta_type = 'beta-fixed'
        elif protein_key.split('_')[4] == '1':
            beta_type = 'beta-variable'

        parameters_name = protein_key.split('_')[3]+'_'+beta_type

        if not protein_name in structured_data.keys():

            parameters = {parameters_name: data[protein_key]}
            method = {method_name: parameters}
            bits = {bits_name: method}

            structured_data[protein_name] = bits

        elif not bits_name in structured_data[protein_name].keys():

            parameters = {parameters_name: data[protein_key]}
            method = {method_name: parameters}

            structured_data[protein_name][bits_name] = method

        elif not method_name in structured_data[protein_name][bits_name].keys():

            parameters = {parameters_name: data[protein_key]}

            structured_data[protein_name][bits_name][method_name] = parameters

        else:

            structured_data[protein_name][bits_name][method_name][parameters_name] = data[protein_key]

    return structured_data