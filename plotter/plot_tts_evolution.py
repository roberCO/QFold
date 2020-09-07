from bokeh.models import ColumnDataSource
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