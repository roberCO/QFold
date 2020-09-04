import utils

import json
import numpy as np

#Read config file with the QFold configuration variables
config_path = './config/config.json'

tools = utils.Utils(config_path)
config_variables = tools.get_config_variables()

# list elements to read
input_files = [
    'glycylglycylglycine_GGG_1_minifold_20_1000'
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
from bokeh.io import output_file, show
from bokeh.models import ColumnDataSource
from bokeh.palettes import GnBu3, OrRd3
from bokeh.plotting import figure

output_file("stacked_split.html")

fruits = ['Apples', 'Pears', 'Nectarines', 'Plums', 'Grapes', 'Strawberries']
years = ["2015", "2016", "2017"]

exports = {'fruits' : fruits,
           '2015'   : [2, 1, 4, 3, 2, 4],
           '2016'   : [5, 3, 4, 2, 4, 6],
           '2017'   : [3, 2, 4, 4, 5, 3]}
imports = {'fruits' : fruits,
           '2015'   : [-1, 0, -1, -3, -2, -1],
           '2016'   : [-2, -1, -3, -1, -2, -2],
           '2017'   : [-1, -2, -1, 0, -2, -2]}

p = figure(y_range=fruits, plot_height=250, x_range=(-16, 16), title="Fruit import/export, by year",
           toolbar_location=None)

p.hbar_stack(years, y='fruits', height=0.9, color=GnBu3, source=ColumnDataSource(exports),
             legend_label=["%s exports" % x for x in years])

p.hbar_stack(years, y='fruits', height=0.9, color=OrRd3, source=ColumnDataSource(imports),
             legend_label=["%s imports" % x for x in years])

p.y_range.range_padding = 0.1
p.ygrid.grid_line_color = None
p.legend.location = "top_left"
p.axis.minor_tick_line_color = None
p.outline_line_color = None

show(p)



# generate plot of the evolution of tts with different steps comparing classical vs quantum

# generate plot of the comparison between quantum and classical difference of tts  