from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Turbo256
import random

def plot_q_vs_c(data):
    
    plot_tts = figure(
        title='TTS comparison quantum vs minifold ', 
        x_axis_label='TTS', 
        y_axis_label='relation classical/quantum', 
        x_axis_type="log", 
        y_axis_type="log", 
        y_range=(0, 10**6), 
        x_range=(0, 10**3)
        )

    plot_tts.xaxis.fixed_location = 1

    for protein_key in data:

        min_tts_q = data[protein_key]['min_tts_q']
        min_tts_c = data[protein_key]['min_tts_c']

        relation = min_tts_c / min_tts_q
        min_tts = min(min_tts_q, min_tts_c)

        plot_tts.square(min_tts, relation, line_width=2, line_color=random.sample(Turbo256,15), legend_label=protein_key)


    plot_tts.yaxis.major_label_orientation = "vertical"   
    plot_tts.xaxis.major_label_orientation = "horizontal"

    plot_tts.quad(top=[10**6], bottom=[1], left=[1],right=[10**3], color="green", fill_alpha=0.1)
    plot_tts.quad(top=[1], bottom=[0], left=[0], right=[10**3], color="red", fill_alpha=0.1)
    plot_tts.title.align = 'center'


    show(plot_tts)