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

#These two nested loops are hardcoded (it could be n nested loops, 1 per AA) because QFold is going to be used just with two and three aminoacids
#if it is scale to more aminoacids, it should be necessary to implement a recursive function
for x in range(0, rotationSteps):

    print('<!> Rotating phi '+ str(anglePhi) +'!')

    for y in range(0, rotationSteps):

        print('<@> Rotating psi '+str(anglePsi)+'!')

        print('\n------------ DISTANCES BEFORE ROTATION --------------')
        for atomConn in carboxyAtom.linked_to:
            print('Distance between carboxy atom (id: '+str(carboxyAtom.atomId)+') and atom (id: ' + str(atomConn.atomId) + ') : ' + str(tools.distance(carboxyAtom, atomConn)))
        print('------------ ------------------------ --------------')

        #Perform the rotations over a copy
        copied_atoms = copy.deepcopy(atoms)
        copied_nitroAtom = tools.findAtom(copied_atoms, 'N', '', nitroConnections)
        copied_carboxyAtom = tools.findAtom(copied_atoms, '', 'Carboxy', carboxyConnections)

        #Always rotate from state (0,0)
        tools.rotate('phi', x * anglePhi, copied_nitroAtom) 

        tools.rotate('psi', y * anglePsi, copied_carboxyAtom)

        print('\n++++++++++++ DISTANCES AFTER ROTATION ++++++++++++++')
        for atomConn in copied_carboxyAtom.linked_to:
            print('Distance between carboxy atom (id: '+str(copied_carboxyAtom.atomId)+') and atom (id: ' + str(atomConn.atomId) + ') : ' + str(tools.distance(carboxyAtom, atomConn)))
        print('++++++++++++ ++++++++++++++++++++++++ ++++++++++++++')
        
        tools.plotting(atoms, 'phi: ' + str(anglePhi) + ' psi: ' + str(anglePsi))

        #Write the file with the actual rotations
        psi.writeFileEnergies(copied_atoms, inputFilenameEnergyPSI4)

        #Calculate the energy of the actual rotations using PSI4
        psi.executePsiCommand(inputFilenameEnergyPSI4, outputFilenameEnergyPSI4)

        #Read the PSI4 output file and get the energy
        energy = psi.readEnergyFromFile(outputFilenameEnergyPSI4)

        print('<!>  Nitro atom   <!>: ' + str(nitroAtom.element) + ' | x => ' + str(nitroAtom.x) + ' y => ' + str(nitroAtom.y) + ' z => ' + str(nitroAtom.z))
        print('<!> Carboxy atom  <!>: ' + str(carboxyAtom.element) + ' | x => ' + str(carboxyAtom.x) + ' y => ' + str(carboxyAtom.y) + ' z => ' + str(carboxyAtom.z)+'\n')
        print ('Energy => ' + str(energy)+'\n----------------------------------------------------------------\n\n')

        # We eliminate previous copies
        del copied_atoms
        del copied_carboxyAtom
        del copied_nitroAtom
