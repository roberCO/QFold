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
from bokeh.layouts import row

precision_minifold = [10, 20, 30, 40, 50, 60]
precision_random = [5, 10, 15, 20, 25, 30] 
tts_minifold = [-5, -10, -15, -20, -25, -30]
tts_random =  [-10, -20, -30, -40, -50, -60]


output_file("random_vs_minifold.html")

proteins = ['gg_3', 'gg_4', 'gc_3', 'gc_4', 'aa_3', 'aa_4']
metrics = ["minifold", "random", "tts_minifold", "tts_random"]

y_precision = [ (protein, metric) for protein in proteins for metric in metrics[:2]]
y_tts = [ (protein, metric) for protein in proteins for metric in metrics[2:]]

counts_precision = sum(zip(precision_minifold, precision_random), ())
counts_tts = sum(zip(tts_minifold, tts_random), ())

color_precision = ['green','green','green','green','green','green','green','green','green','green','green','green']
color_tts = ['red', 'red','red','red','red','red','red','red','red','red','red','red']
source_precision = ColumnDataSource(data=dict(y=y_precision, counts=counts_precision, color=color_precision))
source_tts = ColumnDataSource(data=dict(y=y_tts, counts=counts_tts, color=color_tts))


plot_precision = figure(x_range=(0, 100), y_range=FactorRange(*y_precision), plot_height=250, title="Precision comparison  minifold vs random",
            toolbar_location=None)
plot_tts = figure(x_range=(-100, 0), y_range=FactorRange(*y_tts), plot_height=250, title="TTS comparison minifold vs random",
            toolbar_location=None)

#p.hbar_stack(metrics, y='proteins', height=0.9, color='#008000', source=ColumnDataSource(precision), legend_label= precision_names)
plot_precision.hbar(y='y', right='counts', height=0.9, fill_color='color', source=source_precision)

#p.hbar_stack(metrics, y='proteins', height=0.9, color='#ff0000', source=ColumnDataSource(tts), legend_label=tts_names)
plot_tts.hbar(y='y', right='counts', height=0.9, fill_color='color', source=source_tts)


plot_precision.y_range.range_padding = 0.1
plot_precision.ygrid.grid_line_color = None
#plot_precision.legend.location = "top_left"
plot_precision.axis.minor_tick_line_color = None
plot_precision.outline_line_color = None

plot_tts.y_range.range_padding = 0.1
plot_tts.ygrid.grid_line_color = None
#plot_tts.legend.location = "top_left"
plot_tts.axis.minor_tick_line_color = None
plot_tts.outline_line_color = None
plot_tts.yaxis.visible = False


show(row(plot_tts, plot_precision))


# generate plot of the evolution of tts with different steps comparing classical vs quantum

# generate plot of the comparison between quantum and classical difference of tts  