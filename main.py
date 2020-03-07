import sys
import utils
import atom
import numpy as np
import quantumProcessor
import math

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)

#Global variable
tools = utils.Utils()
qProcessor = quantumProcessor.QuantumProcessor()
beta = 1
scaling_factor = 80 # Modify this parameter to make it reasonable --------

proteinName = sys.argv[1]
numberBitsRotation = int(sys.argv[2])
rotationSteps = pow(2, int(numberBitsRotation))

#Check if it existes a precalculated energy file with the same parameters
#The format should be energies[proteinName][numberBitsForRotation] ex: energiesGlycylglycine2.json

try:
    f = open('./precalculated_energies/energies_'+proteinName+'_'+str(numberBitsRotation)+'.json')
    f.close()
except IOError:
    print('No precalculated energies file found. Calculating energies')
    tools.calculateEnergies(proteinName, numberBitsRotation)

energyList = [[0 for x in range(rotationSteps)] for y in range(rotationSteps)] 
energyList = tools.readEnergyJson(proteinName, numberBitsRotation)

'''
for x in range(len(energyList)):
    for y in range(len(energyList[x])):
        print('Phi: ' + str(x) + ' Psi: ' + str(y) + ' Energy: ' + str(energyList[x][y]))

print('⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤')
'''

truthTableList = []
for x in range(len(energyList)):
    for y in range(len(energyList[x])):

        energyReference = energyList[x][y]

        #Phi +1
        phiValue = (x+1) % rotationSteps
        deltaEnergy = (energyList[phiValue][y] - energyReference) * scaling_factor
        probability = min(1, math.exp(-1*deltaEnergy*beta))
        truthTableList.append([phiValue, y, 0, 1, probability])

        #Phi -1
        phiValue = (x-1) % rotationSteps
        deltaEnergy = (energyList[phiValue][y] - energyReference) * scaling_factor
        probability = min(1, math.exp(-1*deltaEnergy*beta))
        truthTableList.append([phiValue, y, 0, 0, probability])

        #Psi +1
        psiValue = (y+1) % rotationSteps
        deltaEnergy = (energyList[x][psiValue] - energyReference) * scaling_factor
        probability = min(1, math.exp(-1*deltaEnergy*beta))
        truthTableList.append([x, psiValue, 1, +1, probability])

        #Psi -1
        psiValue = (y-1) % rotationSteps
        deltaEnergy = (energyList[x][psiValue] - energyReference) * scaling_factor
        probability = min(1, math.exp(-1*deltaEnergy*beta))
        truthTableList.append([x, psiValue, 1, 0, probability])

#sort truthTableList using as key phi, psi and m
tools.sortByAngleMovements(truthTableList)

for inputValue in truthTableList:
    print('Phi angle: ' + tools.number2binary(inputValue[0], numberBitsRotation) + ' psi angle: ' + tools.number2binary(inputValue[1], numberBitsRotation) + ' rotatedAngle: ' + tools.number2binary(inputValue[2], 1) + ' rotation value: ' + tools.number2binary(inputValue[0], 1) + ' probability:' + tools.number2binary(inputValue[4], 10))


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