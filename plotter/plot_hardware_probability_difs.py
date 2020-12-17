import random
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import json

def plot_hardware_prob_difs(tools):
    '''
    This function should plot a histogram comparing in two different scales the probabilities of measuring on hardware vs noiseless
    to try to see correlation between them
    '''

    # First load the data
    with open('./results/measurements   .json', 'r') as outfile2: 
        dictionary = json.load(outfile2)

    '''
    For each dipeptide, we want to get the following values:
    
    - For non-zero beta, the average of the measured values of '00'
    - For non-zero beta, the std of the measured values of '00'

    - For zero beta, the average and standard of measured values

    - For non-zero beta, the noiseless value of '00'
    - The zero beta, the noiseless value of '00' is .25
    '''

    # First we recover the data
    y_avgs = []
    y_stds = []
    y_noiseless = []

    aas = []
    betas = tools.config_variables['betas']
    ibmq_shots = tools.config_variables['ibmq_shots']

    beta_zero_meas_avg = np.average(dictionary['--']['0-0']['measurements']['00'])
    beta_zero_meas_std = np.std(dictionary['--']['0-0']['measurements']['00'])

    for aa in dictionary.keys():
        if aa != '--':

            if len(dictionary[aa]) > 1:
                raise Exception('Too many values of beta to choose from. Choose manually in the plotting file')

            # We already introduce the baselines

            print('Noiseless dictionary[aa][str(betas[0]) + - +str(betas[1])][noiseless][00]', dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['noiseless']['00'])

            aas.append(aa)

            y_avgs.append((np.average(dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['measurements']['00']) - beta_zero_meas_avg)/ibmq_shots)
            y_stds.append((np.std(dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['measurements']['00']) + beta_zero_meas_std)/ibmq_shots)
            
            y_noiseless.append(dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['noiseless']['00'] - .25)

    for i in range(len(aas)): 
        print('aas',aas[i])
        print('y_noiseless',y_noiseless[i]) 
        print('(y_avgs,y_stds',y_avgs[i],y_stds[i])

    # Then we plot the results

    # Create some mock data
    fig, ax1 = plt.subplots(1, 1, tight_layout=True)

    color = 'tab:blue'
    ax1.set_xlabel('Dipeptides')
    ax1.set_ylabel('Noiseless probability gap', color=color)

    ax1.scatter(aas, y_noiseless, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:red'
    ax2.set_ylabel('Measured probability gap', color=color)  # we already handled the x-label with ax1
    ax2.errorbar(aas, y_avgs, yerr=y_stds, color=color, marker='s', fmt='.')
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Lower part of the figure

    model = np.polynomial.polynomial.polyfit(y_noiseless, y_avgs, 1)
    (a, b) = model

    predict = np.poly1d(np.flip(model))

    r2 = r2_score(y_avgs, predict(y_noiseless)) # How to write it in the plot?

    max_noiseless = int(np.floor(max(y_noiseless)))
    min_noiseless = int(np.ceil(min(y_noiseless)))

    x_lin_reg = range(0,1)

    y_lin_reg = predict(x_lin_reg)
    '''
    ax1[].set_xlim = ([.2,.3])
    ax1[1].set_ylim = ([.15,.35])

    ax1[1].set_xlabel('Noiseless probability gap')
    ax1[1].set_ylabel('Measured probability gap')

    ax1[1].errorbar(y_noiseless, y_avgs, yerr=y_stds)
    ax1[1].plot(x_lin_reg, y_lin_reg, c = 'r')
    '''

    ax1.annotate("r-squared = {:.3f}".format(r2), (0, 1))

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    #plt.savefig(fname = './results/hardware_measurements')
