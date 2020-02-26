import sys
import utils
import psiFour
import atom
import numpy as np
import quantumProcessor

if(len(sys.argv) != 3):
    print ("<*> ERROR: Wrong number of parameters - Usage: python main.py ProteinName numberBitsForRotations")
    print ("<!> Example: python main.py Glycylglycine 6 (6 bits for rotations are 64 steps)")
    sys.exit(0)

#Global variable
tools = utils.Utils()

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
anglePhi = 0

anglesEnergy = []
#These two nested loops are hardcoded (it could be n nested loops, 1 per AA) because QFold is going to be used just with two and three aminoacids
#if it is scale to more aminoacids, it should be necessary to implement a recursive function
for x in range(0, rotationSteps):

    tools.rotate('phi', anglePhi, nitroAtom)
    anglePsi = 0

    for y in range(0, rotationSteps):

        tools.rotate('psi', anglePsi, carboxyAtom)
        #tools.plotting(atoms, 'phi: ' + str(anglePhi) + ' psi: ' + str(anglePsi))


        #Write the file with the actual rotations
        psi.writeFileEnergies(atoms, inputFilenameEnergyPSI4)

        #Calculate the energy of the actual rotations using PSI4
        psi.executePsiCommand(inputFilenameEnergyPSI4, outputFilenameEnergyPSI4)

        #Read the PSI4 output file and get the energy
        energy = psi.readEnergyFromFile(outputFilenameEnergyPSI4)
        normalizedEnergy = energy*-1
        normalizedEnergy = "%.6f" % normalizedEnergy

        anglesEnergy.append([anglePhi, anglePsi, normalizedEnergy])
        anglePsi += 1/rotationSteps

        #print ('⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤\n⬤   Phi: ' + str(anglePhi) +'\n⬤   Psi: '+ str(anglePsi)+ '\n⬤   Energy: ' + str(energy) +'\n⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤\n\n')

    anglePhi += 1/rotationSteps

for element in anglesEnergy:
    
    print("Energy: " + str(tools.floatToBits(float(element[2]))) + " phi: " + str(tools.floatToBits(element[0])) + " psi: " + str(tools.floatToBits(element[1])))
    #print("Energy: " + str(tools.floatToBits(float(element[2]))) + " phi: " + str(tools.floatToBits(element[0])) + " psi: " + str(tools.floatToBits(element[1])))