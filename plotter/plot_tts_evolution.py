from bokeh.plotting import figure, output_file, show

def plot_tts_ev(data):
    
    for protein_key in data.keys():

        output_file(protein_key + "_tts_evolution.html")

        x = [x_iter for x_iter in range(data[protein_key]['initial_step'], data[protein_key]['final_step'])]

        quantum_tts = data[protein_key]['quantum_tts']
        classical_tts = data[protein_key]['classical_tts']

        plot_tts = figure(title='Evolution of tts with different steps', x_axis_label='STEPS', y_axis_label='TTS', y_axis_type="log", y_range=(1, 10**6))
        plot_tts.line(x, quantum_tts, line_width=2, line_color='green', legend_label='Quantum')
        plot_tts.line(x, classical_tts, line_width=2, line_color='orange', legend_label='Classical')
        plot_tts.title.align = 'center'

        show(plot_tts)

def plot_tts_ev_bits(data):

    data_protein = {}

    for protein_key in data:

        protein_name = protein_key.split('_')[0]
        bits = protein_key.split('_')[1]
        method = protein_key.split('_')[2]

        tts_values = {'q_tts': data[protein_key]['min_tts_q'], 'c_tts': data[protein_key]['min_tts_c']}

        if not protein_name in data_protein.keys():

            bits_dict = {bits: tts_values}
            method_dict = {method: bits_dict}

            data_protein[protein_name] = method_dict

        elif not method in data_protein[protein_name]:

            method_dict = {bits: tts_values}

            data_protein[protein_name][method] = bits_dict

        else:

            data_protein[protein_name][method][bits] = tts_values


    for protein_key in data_protein:

        for method_key in data_protein[protein_key].keys():

            x_labels = []
            quantum_tts = []
            classical_tts = []

            for bits_key in data_protein[protein_key][method_key].keys():

                x_labels.append(int(bits_key))            
                quantum_tts.append(data_protein[protein_key][method_key][bits_key]['q_tts'])
                classical_tts.append(data_protein[protein_key][method_key][bits_key]['c_tts'])

            output_file(protein_key + '_' + method_key + '_' + "_tts_evolution.html")

            plot_tts = figure(title='Evolution of tts with different bits', x_axis_label='Bits', y_axis_label='TTS', y_axis_type="log", y_range=(1, 10**5))
            plot_tts.line(x=x_labels, y=quantum_tts, line_width=2, line_color='green', legend_label='Quantum')
            plot_tts.line(x=x_labels, y=classical_tts, line_width=2, line_color='orange', legend_label='Classical')
            plot_tts.title.align = 'center'

            show(plot_tts)