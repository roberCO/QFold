import sys
import utils
import psiFour
import atom
import numpy as np
import copy

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)

#Global variable
tools = utils.Utils()

proteinName = sys.argv[1]
rotationSteps = pow(2, int(sys.argv[2]))

# Scaling factor

# Get the minimum angle rotation
anglePhi = 1/rotationSteps
anglePsi = 1/rotationSteps

n_iterations = 100000
scaling_factor = 5000 # Modify this parameter to make it reasonable --------

#call psi4 to get the atoms of the protein
psi = psiFour.PsiFour()
atoms = psi.getAtomsFromProtein(proteinName)

#Calculate the connection between atoms
atoms = tools.calculateAtomConnection(atoms)

nitroConnections = [['C', 2]]
carboxyConnections = [['C', 1], ['O', 2]]

nitroAtom = tools.findAtom(atoms, 'N', '', nitroConnections)
carboxyAtom = tools.findAtom(atoms, '', 'Carboxy', carboxyConnections)

inputFilenameEnergyPSI4 = 'inputRotations'
outputFilenameEnergyPSI4 = 'outputRotations'

# Create a dictionary of already calculated energies 
energies_dictionary = {}

for Phi in range (rotationSteps):
    for Psi in range (rotationSteps):
        energies_dictionary[(Phi,Psi)] = None


# Random starting combination of angles
anglePsi_old = np.random.choice(rotationSteps)
anglePhi_old = np.random.choice(rotationSteps)

# Calculate energy for the first state:
# From the initial state deepcopy the molecule
copied_atoms = copy.deepcopy(atoms)
copied_nitroAtom = tools.findAtom(copied_atoms, 'N', '', nitroConnections)
copied_carboxyAtom = tools.findAtom(copied_atoms, '', 'Carboxy', carboxyConnections)

#rotate according to the angles anglePhi_old and anglePsi_old. Notice that anglePhi is 1/rotationSteps and anglePhi_old is between (0, rotationSteps). Same for Psi
tools.rotate('phi', anglePhi_old * anglePhi, copied_nitroAtom) 
tools.rotate('psi', anglePsi_old * anglePsi, copied_carboxyAtom)

#Write the file with the actual rotations
psi.writeFileEnergies(copied_atoms, inputFilenameEnergyPSI4)

#Calculate the energy of the actual rotations using PSI4
psi.executePsiCommand(inputFilenameEnergyPSI4, outputFilenameEnergyPSI4)

#Read the PSI4 output file and get the energy
E_old = psi.readEnergyFromFile(outputFilenameEnergyPSI4)

del copied_atoms
del copied_carboxyAtom
del copied_nitroAtom




for iteration in range(n_iterations):

    # First retrieve the present energy
    E_old = energies_dictionary[(anglePhi_old,anglePsi_old)]
    
    # Propose a change
    change_angle = np.random.choice(('phi','psi'))
    change_plus_minus = np.random.choice((1,-1))

    # Calculate the new angles
    if change_angle == 'phi':
        anglePhi_new = anglePhi_new + change_plus_minus * anglePhi
        anglePsi_new = anglePsi_new
    elif change_angle == 'psi':
        anglePhi_new = anglePhi_new
        anglePsi_new = anglePsi_new + change_plus_minus * anglePsi
    
    # Calculate the new energy
    if energies_dictionary((anglePhi_new,anglePsi_new)) != None:
        E_new = energies_dictionary[(anglePhi_new,anglePsi_new)]

    else:
        # From the initial state deepcopy the molecule
        copied_atoms = copy.deepcopy(atoms)
        copied_nitroAtom = tools.findAtom(copied_atoms, 'N', '', nitroConnections)
        copied_carboxyAtom = tools.findAtom(copied_atoms, '', 'Carboxy', carboxyConnections)

        #rotate according to the angles anglePhi_old and anglePsi_old. Notice that anglePhi is 1/rotationSteps and anglePhi_old is between (0, rotationSteps). Same for Psi
        tools.rotate('phi', anglePhi_old * anglePhi, copied_nitroAtom) 
        tools.rotate('psi', anglePsi_old * anglePsi, copied_carboxyAtom)

        #Write the file with the actual rotations
        psi.writeFileEnergies(copied_atoms, inputFilenameEnergyPSI4)

        #Calculate the energy of the actual rotations using PSI4
        psi.executePsiCommand(inputFilenameEnergyPSI4, outputFilenameEnergyPSI4)

        #Read the PSI4 output file and get the energy
        E_new = psi.readEnergyFromFile(outputFilenameEnergyPSI4)
        energies_dictionary[(anglePhi_new,anglePsi_new)] = E_new

        del copied_atoms
        del copied_carboxyAtom
        del copied_nitroAtom

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


