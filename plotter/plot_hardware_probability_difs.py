import random
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import json
#import matplotlib

def plot_hardware_prob_difs(tools):
    '''
    This function should plot a histogram comparing in two different scales the probabilities of measuring on hardware vs noiseless
    to try to see correlation between them
    '''

    # First load the data
    with open('./results/measurements.json', 'r') as outfile2: 
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
    y_err = []

    aas = []
    betas = tools.config_variables['betas']
    ibmq_shots = tools.config_variables['ibmq_shots']

    beta_zero_meas_avg = np.average(dictionary['--']['0-0']['measurements']['00'])
    beta_zero_meas_std = np.std(dictionary['--']['0-0']['measurements']['00'])
    beta_zero_runs = len(dictionary['--']['0-0']['measurements']['00'])
    beta_zero_p = beta_zero_meas_avg / ibmq_shots
    beta_zero_err = np.sqrt(beta_zero_p*(1-beta_zero_p)/(ibmq_shots*beta_zero_runs))


    for aa in dictionary.keys():
        if aa != '--':

            if len(dictionary[aa]) > 1:
                raise Exception('Too many values of beta to choose from. Choose manually in the plotting file')

            # We already introduce the baselines

            print('Noiseless dictionary[aa][str(betas[0]) + - +str(betas[1])][noiseless][00]', dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['noiseless']['00'])

            aas.append(aa)

            avg = np.average(dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['measurements']['00'])
            runs = len(dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['measurements']['00'])
            p = avg / ibmq_shots

            y_avgs.append((avg - beta_zero_meas_avg)/ibmq_shots)
            y_stds.append((np.std(dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['measurements']['00']) + beta_zero_meas_std)/ibmq_shots)
            y_err.append(np.sqrt(p*(1-p)/(ibmq_shots*runs))+ beta_zero_err)
            y_noiseless.append(dictionary[aa][str(betas[0]) + '-' +str(betas[1])]['noiseless']['00'])

    for i in range(len(aas)): 
        print('aas',aas[i])
        print('y_noiseless',y_noiseless[i]) 
        print('(y_avgs,y_stds',y_avgs[i],y_stds[i])

    # Then we plot the results

    # Create some mock data
    fig, ax1 = plt.subplots(1, 1, tight_layout=True)

    '''#color = 'tab:blue'
    ax1.set_xlabel('Dipeptides')
    ax1.set_ylabel('Noiseless probability gap')#, color=color)'''

    color = 'tab:blue'
    ax1.set_xlabel('Dipeptides', fontsize = 12)
    ax1.hlines(0, xmin = -.2, xmax = len(aas), color = 'red', linestyles = 'dashed')
    ax1.set_ylabel(r'Measured probability gap: $p(\beta(\vec{t}) = (0.1,1))  - p(\beta(\vec{t}) = (0,0))$', fontsize = 10)  # we already handled the x-label with ax1
    ax1.errorbar(aas, y_avgs, yerr=y_err, marker='d', fmt='.')
    #ax1.tick_params(axis='y')
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)

    # this is an inset axes over the main axes
    a = plt.axes([.55, .3, .4, .2])
    

    xr = np.arange(len(aas))  # the label locations
    width = 0.35  # the width of the bars
    a.hlines(.25, xmin = xr[0]-width, xmax = xr[-1], color = 'red', linestyles = 'dashed')

    bars = a.bar(xr - width/2, y_noiseless, width, label='Noiseless probabilities', color = 'green')
    #a.set_xlabel('Dipeptides')
    #a.set_ylabel('Noiseless probabilities', color=color)
    #a.set_xticklabels(aas)
    a.set_yticklabels([0, 0.25, 0.5, 0.75])

    plt.title('Noiseless probabilities')
    plt.xticks([])
    #plt.yticks([])
    
    # Lower part of the figure
    '''

    model = np.polynomial.polynomial.polyfit(y_noiseless, y_avgs, 1)
    (a, b) = model

    predict = np.poly1d(np.flip(model))

    r2 = r2_score(y_avgs, predict(y_noiseless)) # How to write it in the plot?

    max_noiseless = int(np.floor(max(y_noiseless)))
    min_noiseless = int(np.ceil(min(y_noiseless)))

    x_lin_reg = range(0,1)

    y_lin_reg = predict(x_lin_reg)
    
    ax1[].set_xlim = ([.2,.3])
    ax1[1].set_ylim = ([.15,.35])

    ax1[1].set_xlabel('Noiseless probability gap')
    ax1[1].set_ylabel('Measured probability gap')

    ax1[1].errorbar(y_noiseless, y_avgs, yerr=y_stds)
    ax1[1].plot(x_lin_reg, y_lin_reg, c = 'r')
    

    ax1.annotate("r-squared = {:.3f}".format(r2), (0, 1))
    '''
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    plt.savefig(fname = './results/hardware_measurements')
