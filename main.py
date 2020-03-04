import sys
import utils
import psiFour
import atom
import numpy as np
import quantumProcessor
import copy
import math

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)

#Global variable
tools = utils.Utils()
qProcessor = quantumProcessor.QuantumProcessor()
beta = 1
scaling_factor = 200 # Modify this parameter to make it reasonable --------

proteinName = sys.argv[1]
rotationSteps = pow(2, int(sys.argv[2]))

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

anglePhi = 1/rotationSteps
anglePsi = 1/rotationSteps

energyList = [[0 for x in range(rotationSteps)] for y in range(rotationSteps)] 
anglesEnergy = []
#These two nested loops are hardcoded (it could be n nested loops, 1 per AA) because QFold is going to be used just with two and three aminoacids
#if it is scale to more aminoacids, it should be necessary to implement a recursive function
for x in range(0, rotationSteps):

    for y in range(0, rotationSteps):

        #Perform the rotations over a copy
        copied_atoms = copy.deepcopy(atoms)
        copied_nitroAtom = tools.findAtom(copied_atoms, 'N', '', nitroConnections)
        copied_carboxyAtom = tools.findAtom(copied_atoms, '', 'Carboxy', carboxyConnections)

        #Always rotate from state (0,0)
        tools.rotate('phi', x * anglePhi, copied_nitroAtom) 

        tools.rotate('psi', y * anglePsi, copied_carboxyAtom)
        
        #tools.plotting(atoms, 'phi: ' + str(anglePhi) + ' psi: ' + str(anglePsi))

        #Write the file with the actual rotations
        psi.writeFileEnergies(copied_atoms, inputFilenameEnergyPSI4)

        #Calculate the energy of the actual rotations using PSI4
        psi.executePsiCommand(inputFilenameEnergyPSI4, outputFilenameEnergyPSI4)

        #Read the PSI4 output file and get the energy
        energy = psi.readEnergyFromFile(outputFilenameEnergyPSI4)
        normalizedEnergy = energy*-1
        normalizedEnergy = "%.6f" % normalizedEnergy

        anglesEnergy.append([anglePhi, anglePsi, normalizedEnergy])
        anglePsi += 1/rotationSteps

        energyList[x][y] = energy

        #print ('⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤\n⬤   Phi('+str(x)+'): ' + str(anglePhi) +'\n⬤   Psi('+str(y)+'): '+ str(anglePsi)+ '\n⬤   Energy: ' + str(energy) +'\n⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤\n\n')

        # We eliminate previous copies
        del copied_atoms
        del copied_carboxyAtom
        del copied_nitroAtom

    anglePhi += 1/rotationSteps

for x in range(len(energyList)):
    for y in range(len(energyList[x])):
        print('Phi: ' + str(x) + ' Psi: ' + str(y) + ' Energy: ' + str(energyList[x][y]))

print('⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤')

truthTableList = []
for x in range(len(energyList)):
    for y in range(len(energyList[x])):

        energyReference = energyList[x][y]

        #Phi +1
        phiValue = (x+1) % rotationSteps
        deltaEnergy = (energyList[phiValue][y] - energyReference) * scaling_factor
        print('Delta energy: ' + str(deltaEnergy))
        probability = min(1, math.exp(-1*deltaEnergy*beta))
        truthTableList.append([phiValue, y, 0, +1, probability])

        #Phi -1
        phiValue = (x-1) % rotationSteps
        deltaEnergy = (energyList[phiValue][y] - energyReference) * scaling_factor
        print('Delta energy: ' + str(deltaEnergy))
        probability = min(1, math.exp(-1*deltaEnergy*beta))
        truthTableList.append([phiValue, y, 0, -1, probability])

        #Psi +1
        psiValue = (y+1) % rotationSteps
        deltaEnergy = (energyList[x][psiValue] - energyReference) * scaling_factor
        print('Delta energy: ' + str(deltaEnergy))
        probability = min(1, math.exp(-1*deltaEnergy*beta))
        truthTableList.append([x, psiValue, 1, +1, probability])

        #Psi -1
        psiValue = (y-1) % rotationSteps
        deltaEnergy = (energyList[x][psiValue] - energyReference) * scaling_factor
        print('Delta energy: ' + str(deltaEnergy))
        probability = min(1, math.exp(-1*deltaEnergy*beta))
        truthTableList.append([x, psiValue, 1, -1, probability])

for inputValue in truthTableList:
    print('Phi angle: ' + str(inputValue[0]) + ' psi angle: ' + str(inputValue[1]) + ' rotatedAngle: ' + str(inputValue[2]) + ' rotation value: ' + str(inputValue[3]) + ' probability:' + str(inputValue[4]))


'''
binaryInputOracle = []
for element in anglesEnergy:
    
    phiBinary = tools.floatToBits(element[0])
    psiBinary = tools.floatToBits(element[1])
    energyBinary = tools.floatToBits(float(element[2]))
    elementBinary = str(phiBinary) + str(psiBinary) + str(energyBinary)
    print("Energy: " + str(tools.floatToBits(float(element[2]))) + " phi: " + str(tools.floatToBits(element[0])) + " psi: " + str(tools.floatToBits(element[1])))
    print("Binary code: " + elementBinary)

    binaryInputOracle += [elementBinary]

qProcessor.inputListOracle(binaryInputOracle)
'''