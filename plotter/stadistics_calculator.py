import numpy as np
from prettytable import PrettyTable
from collections import OrderedDict

def calculate_stats(data):

    mean_precisions = calculate_mean_precision(data)
    tts_tables = calculate_best_tts_parameters(data)

    return [mean_precisions, tts_tables]

def calculate_mean_precision(data):

    stats = {}

    for key in data.keys():

        number_aas = len(key.split('_')[0])
        number_bits = key.split('_')[1]
        init = key.split('_')[2]

        if number_aas in stats.keys():

            if number_bits in stats[number_aas].keys():

                if init in stats[number_aas][number_bits].keys():

                    accum_q_tts = stats[number_aas][number_bits][init]['q_tts'] + data[key]['min_tts_q']
                    accum_c_tts = stats[number_aas][number_bits][init]['c_tts'] + data[key]['min_tts_c']

                    accum_precision = stats[number_aas][number_bits][init]['precision'] + data[key]['precision']
                    counter = stats[number_aas][number_bits][init]['counter'] + 1

                    stats[number_aas][number_bits][init].update({'q_tts': accum_q_tts, 'c_tts': accum_c_tts, 'precision': accum_precision, 'counter': counter})

                else:
                    stats[number_aas][number_bits].update({init:{'q_tts': data[key]['min_tts_q'], 'c_tts': data[key]['min_tts_c'], 'precision': data[key]['precision'], 'counter': 1}})

            else:
                stats[number_aas].update({number_bits:{init:{'q_tts': data[key]['min_tts_q'], 'c_tts': data[key]['min_tts_c'], 'precision': data[key]['precision'], 'counter': 1}}})

        else:
            stats.update({number_aas:{number_bits:{init:{'q_tts': data[key]['min_tts_q'], 'c_tts': data[key]['min_tts_c'], 'precision': data[key]['precision'], 'counter': 1}}}})

    for aas in stats.keys():
        for bits in stats[aas].keys():
            for init in stats[aas][bits].keys():

                stats[aas][bits][init]['q_tts'] = stats[aas][bits][init]['q_tts']/stats[aas][bits][init]['counter']
                stats[aas][bits][init]['c_tts'] = stats[aas][bits][init]['c_tts']/stats[aas][bits][init]['counter']
                stats[aas][bits][init]['precision'] = stats[aas][bits][init]['precision']/stats[aas][bits][init]['counter']

    return stats
    
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