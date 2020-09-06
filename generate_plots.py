import utils

import json
import numpy as np
import bokeh
from bokeh.io import export_png

#Read config file with the QFold configuration variables
config_path = './config/config.json'

tools = utils.Utils(config_path)
config_variables = tools.get_config_variables()

# list elements to read
input_files = [
    'glycylglycylglycine_GGG_1_minifold_20_1000',
    'glycylglycylglycine_GGG_1_random_20_1000'

]

for input_name in input_files:

    data = {}
    # read data
    with open(config_variables['path_tts_plot']+ 'tts_results_'+input_name+'.json') as json_file:
                data[input_name] = json.load(json_file)

# prepare the data
precisions = {}
min_tts = {}

for protein_key in data.keys():

    aas = protein_key.split('_')[1]
    bits = protein_key.split('_')[2]
    init_mode = protein_key.split('_')[3]
    phi_prec = data[protein_key]['initialization_stats']['phis_precision']
    psi_prec = data[protein_key]['initialization_stats']['psis_precision']
    # divided by 2 because it is the mean between phi and psi and by 100 to normalize the precision
    precisions[aas+'_'+bits+'_'+init_mode] = np.mean(np.mean(phi_prec) + np.mean(psi_prec))/2/100
    min_tts[aas+'_'+bits+'_'+init_mode] = min(data[protein_key]['final_stats']['q']['value'], data[protein_key]['final_stats']['c']['value'])

# generate plot of minifold vs random inizialization mode
# https://towardsdatascience.com/interactive-bar-charts-with-bokeh-7230e5653ba3

from bokeh.io import output_file, show
from bokeh.models import ColumnDataSource
from bokeh.palettes import Set1
from bokeh.plotting import figure
import random
from bokeh.palettes import Turbo256
from bokeh.models import FactorRange


precision_minifold = [10, 20, 30, 40, 50, 60]
precision_random = [5, 10, 15, 20, 25, 30] 
tts_minifold = [-5, -10, -15, -20, -25, -30]
tts_random =  [-10, -20, -30, -40, -50, -60]


output_file("random_vs_minifold.html")

proteins = ['gg_3', 'gg_4', 'gc_3', 'gc_4', 'aa_3', 'aa_4']
metrics = ["precision_minifold", "precision_random", "tts_minifold", "tts_random"]

y_precision = [ (protein, metric) for protein in proteins for metric in metrics[:2]]
y_tts = [ (protein, metric) for protein in proteins for metric in metrics[2:]]

counts_precision = sum(zip(precision_minifold, precision_random), ())
counts_tts = sum(zip(tts_minifold, tts_random), ())

source_precision = ColumnDataSource(data=dict(y=y_precision, counts_precision=counts_precision, color=random.sample(Turbo256,12)))
source_tts = ColumnDataSource(data=dict(y=y_tts, counts_tts=counts_tts, color=random.sample(Turbo256,12)))

p = figure(x_range=(-100, 100), y_range=FactorRange(*y_precision), plot_height=250, title="Precision and TTS comparison for minifold an random initialization",
            toolbar_location=None)

precision_names = metrics[0]
tts_names = metrics[2]

#p.hbar_stack(metrics, y='proteins', height=0.9, color='#008000', source=ColumnDataSource(precision), legend_label= precision_names)
p.hbar(y='y', right='counts_precision', height=0.9, fill_color='color', source=source_precision)

#p.hbar_stack(metrics, y='proteins', height=0.9, color='#ff0000', source=ColumnDataSource(tts), legend_label=tts_names)
#p.hbar(y='y', left='counts_tts', height=0.9, fill_color='color', source=source_tts)


p.y_range.range_padding = 0.1
p.ygrid.grid_line_color = None
#p.legend.location = "top_left"
p.axis.minor_tick_line_color = None
p.outline_line_color = None

show(p)


# generate plot of the evolution of tts with different steps comparing classical vs quantum

# generate plot of the comparison between quantum and classical difference of tts  