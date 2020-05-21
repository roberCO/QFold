import sys
import initializer
import angleCalculator
import psiFour
import utils

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)

#HARDCODED
aminoacids = 'GG'

proteinName = sys.argv[1].lower()
numberBitsRotation = int(sys.argv[2])
rotationSteps = pow(2, int(numberBitsRotation))

#Global variable
scaling_factor = 80 # Modify this parameter to make it reasonable --------
precision_solution = 0.9 #Calculate the solution for this precision value
n_ancilla_bits = 4 # Number of ancilla bits for the quantum metropolis
t = 100 #Number steps

angleInitializer = initializer.Initializer()
angleCalculator = angleCalculator.AngleCalculator(numberBitsRotation, n_ancilla_bits, scaling_factor)
psi = psiFour.PsiFour()
tools = utils.Utils()

#Check if it existes a precalculated energy file with the same parameters, if not call initializer to calculate it
#The format should be energies[proteinName][numberBitsForRotation] ex: energiesGlycylglycine2.json
try:
    f = open('./precalculated_energies/energies_'+proteinName+'_'+str(numberBitsRotation)+'.json')
    f.close()
except IOError:
    print('<!> Info: No precalculated energies file found => Calculating energies')
    angleInitializer.calculateEnergies(proteinName, numberBitsRotation, aminoacids)

#Create an empty list of enery list
#HARDCODED for proteins with only two aminoacids
#TODO modify to any number of aminoacids (it should a list of list, each position of the list contains a list of phi and psi values of this list position)
[energyList, phi_angle_psi4, psi_angle_psi4] = psi.readEnergyJson(proteinName, numberBitsRotation)

quantum_p_t = angleCalculator.calculate3DStructure(energyList, t, 0)[0][0]
classical_matrix = angleCalculator.calculate3DStructure(energyList, t, 1)
classical_p_t = classical_matrix[0][0]

quantum_TTS = tools.calculateTTS(precision_solution, t, quantum_p_t)
classical_TTS = tools.calculateTTS(precision_solution, t, classical_p_t)

print('Quantum Metropolis => TTS:', quantum_TTS)
print('Classical Metropolis => TTS:', classical_TTS)
