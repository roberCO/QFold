import sys
import initializer
import angleCalculator
import psiFour

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)

#HARDCODED
aminoacids = 'GG'

proteinName = sys.argv[1]
numberBitsRotation = int(sys.argv[2])
rotationSteps = pow(2, int(numberBitsRotation))

#Global variable
beta = 1
scaling_factor = 80 # Modify this parameter to make it reasonable --------
angleInitializer = initializer.Initializer()
angleCalculator = angleCalculator.AngleCalculator(rotationSteps, scaling_factor, beta)
psi = psiFour.PsiFour()

#Check if it existes a precalculated energy file with the same parameters, if not call initializer to calculate it
#The format should be energies[proteinName][numberBitsForRotation] ex: energiesGlycylglycine2.json
try:
    f = open('./precalculated_energies/energies_'+proteinName+'_'+str(numberBitsRotation)+'.json')
    f.close()
except IOError:
    print('<!> Info: No precalculated energies file found => Calculating energies')
    angleInitializer.calculateEnergies(proteinName, numberBitsRotation, aminoacids)

#Create an empty list of enery list
#HARDCODED for proteins with only to aminoacids
#TODO modify to any number of aminoacids (it should a list of list, each position of the list contains a list of phi and psi values of this list position)
energyList = [[0 for x in range(rotationSteps)] for y in range(rotationSteps)] 
energyList = psi.readEnergyJson(proteinName, numberBitsRotation)

angleCalculator.calculate3DStructure(energyList)