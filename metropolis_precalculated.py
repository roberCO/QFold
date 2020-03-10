import sys
import utils
import psiFour
import atom
import numpy as np
import copy
import json

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)


### IMPORTANT: how do we read the energies from the json??? Need to save them in energies_dictionary
with open('energies_'+str(sys.argv[1])+'_'+str(sys.argv[2])+'.json') as f:
    d = json.load(f)
    print(d)

#Global variable
tools = utils.Utils()

proteinName = sys.argv[1]
rotationSteps = pow(2, int(sys.argv[2]))

# Scaling factor

# Get the minimum angle rotation
delta_anglePhi = 1/rotationSteps
delta_anglePsi = 1/rotationSteps

n_iterations = 100000
scaling_factor = 5000 # Modify this parameter to make it reasonable --------

# Create a dictionary of already calculated energies 
energies_dictionary = {}


# Random starting combination of angles
anglePsi_old = np.random.choice(rotationSteps)
anglePhi_old = np.random.choice(rotationSteps)


for iteration in range(n_iterations):

    # First retrieve the present energy
    E_old = energies_dictionary[(anglePhi_old,anglePsi_old)]
    
    # Propose a change
    change_angle = np.random.choice(('phi','psi'))
    change_plus_minus = np.random.choice((1,-1))

    # Calculate the new angles
    if change_angle == 'phi':
        anglePhi_new = anglePhi_old + change_plus_minus * delta_anglePhi
        anglePsi_new = anglePsi_old
    elif change_angle == 'psi':
        anglePhi_new = anglePhi_old
        anglePsi_new = anglePsi_old + change_plus_minus * delta_anglePsi

    # Retrieve the new energy
    E_new = energies_dictionary[(anglePhi_new,anglePsi_new)]   
    
    Delta_E = (E_new - E_old) * scaling_factor

    # Lets use a non_optimal simple schedule
    beta = iteration / n_iterations
    probability_threshold = np.exp(-beta*Delta_E)
    random_number = np.random.random_sample

    # We should accept the change if probability_threshold > 1 (the energy goes down) or if beta is small.
    # If beta small, np.exp(-beta*Delta_E) approx 1.
    if random_number < min(1,probability_threshold): # Accept the change
        anglePhi_old = anglePhi_new
        anglePsi_old = anglePsi_new


