import sys
sys.path.insert(1, '/home/roberto/Escritorio/qfold/QFold/')
import utils

from plot_minifold_vs_random import plot_m_vs_r
from plot_tts_evolution import plot_tts_ev
from plot_quantum_vs_classical import plot_q_vs_c
from stadistics_calculator import calculate_stats

import json
import numpy as np
import bokeh
from bokeh.io import export_png

#Read config file with the QFold configuration variables
config_path = './config/config.json'

tools = utils.Utils(config_path)

# list elements to read
input_files = [
    'glycylglycylglycine_GGG_1_minifold_20_1000',
    'glycylglycylglycine_GGG_1_random_20_1000',
    'glycylglycine_GG_1_minifold_10_1000',
    'glycylglycine_GG_1_random_10_1000'
]

results = {}
for input_name in input_files:
    results.update(tools.read_results_data(input_name))

# generate plot of minifold vs random inizialization mode
plot_m_vs_r(results)

# generate plot of the evolution of tts with different steps comparing classical vs quantum
plot_tts_ev(results)

# generate plot of the comparison between quantum and classical difference of tts
plot_q_vs_c(results)

stats = calculate_stats(results)

print(stats)