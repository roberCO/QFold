import sys
sys.path.insert(1, './')
from os import listdir
from os.path import isfile, join
import utils

from plot_minifold_vs_random import plot_m_vs_r
from plot_tts_evolution import plot_tts_ev, plot_tts_ev_bits
from plot_quantum_vs_classical import plot_q_vs_c, plot_q_vs_c_slope, plot_q_opt_step, plot_q_vs_c_slope_var, plot_q_vs_c_tts_ratio, TTSplotter_matplotlib

from plot_hardware_probability_difs import plot_hardware_prob_difs
from stadistics_calculator import calculate_stats

#Read config file with the QFold configuration variables
config_path = './config/config.json'

tools = utils.Utils(config_path)
config_variables = tools.get_config_variables()

# list elements to read
input_files = [f for f in listdir(config_variables['path_tts_plot']) if isfile(join(config_variables['path_tts_plot'], f))]
results = {}
for input_name in input_files:
    if input_name[:3] == 'tts' and input_name.split('.')[1] == 'json':
        results.update(tools.read_results_data(input_name))

# generate plot of minifold vs random initialization mode
plot_m_vs_r(results)

# generate plot of the evolution of tts with different steps comparing classical vs quantum
plot_tts_ev(results)

# generate plot of the comparison between quantum and classical difference of tts
plot_q_vs_c(results)
plot_q_vs_c_slope(results)
plot_q_vs_c_slope_var(results)
plot_q_vs_c_tts_ratio(results)
plot_q_opt_step(results)
TTSplotter_matplotlib(results, ['fixed'])
TTSplotter_matplotlib(results, ['logarithmic', 'linear', 'geometric', 'exponential'])

# generate plot of evolution quantum and classical with different bits
plot_tts_ev_bits(results)

# generate plot of hardware measurements
plot_hardware_prob_difs(tools)

[stats, tts_tables] = calculate_stats(results)

for table in tts_tables:
    print(table, '\n')

print(stats)