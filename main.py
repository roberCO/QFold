import sys
import utils
import psiFour
import atom
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def plotting(list_of_atoms, title):
    
    #-----
    VecStart_x = []
    VecStart_y = []
    VecStart_z = []
    VecEnd_x = []
    VecEnd_y = []
    VecEnd_z  = []
    
    #Make list of conections
    list_of_connections = []
    for at1 in list_of_atoms:
        for at2 in list_of_atoms:
            if (at1,at2) not in list_of_connections and (at2,at1) not in list_of_connections and at2 in at1.linked_to:
                list_of_connections += [(at1,at2)]    
                
    for tupl in list_of_connections:
        VecStart_x += [tupl[0].x]
        VecStart_y += [tupl[0].y]
        VecStart_z += [tupl[0].z]
        VecEnd_x += [tupl[1].x]
        VecEnd_y += [tupl[1].y]
        VecEnd_z  += [tupl[1].z]
        
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    ax = fig.add_subplot(111, projection='3d')
    
    for i in range(len(list_of_connections)):
        ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]],color='grey')
    #-----    
    
    xs = []
    ys = []
    zs = []
    c = []
    for at in list_of_atoms:
        xs += [at.x]
        ys += [at.y]
        zs += [at.z]
        if at.element == 'N':
            c += ['blue']
        elif at.element == 'C':
            c += ['black']
        elif at.element == 'O':
            c += ['red']
        elif at.element == 'H':
            c += ['green']    
    
    ax.scatter(xs, ys, zs,c=c,depthshade= False)

    for i in range(len(xs)): 
        ax.text(xs[i],ys[i],zs[i],  '%s' % (str(i)))
    plt.show()
    return ax

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

#These two nested loops are hardcoded (it could be n nested loops, 1 per AA) because QFold is going to be used just with two and three aminoacids
#if it is scale to more aminoacids, it should be necessary to implement a recursive function
for x in range(0, rotationSteps):

    print('<!> Rotating phi '+ str(anglePhi) +'!')
    tools.rotate('phi', anglePhi, nitroAtom)
    anglePsi = 0

    for y in range(0, rotationSteps):

        print('<@> Rotating psi '+str(anglePsi)+'!')
        tools.rotate('psi', anglePsi, carboxyAtom)
        plotting(atoms, 'phi: ' + str(anglePhi) + ' psi: ' + str(anglePsi))

        anglePsi += 1/rotationSteps

        #Write the file with the actual rotations
        psi.writeFileEnergies(atoms, inputFilenameEnergyPSI4)

        #Calculate the energy of the actual rotations using PSI4
        psi.executePsiCommand(inputFilenameEnergyPSI4, outputFilenameEnergyPSI4)

        #Read the PSI4 output file and get the energy
        energy = psi.readEnergyFromFile(outputFilenameEnergyPSI4)

        print('<!>  Nitro atom   <!>: ' + str(nitroAtom.element) + ' | x => ' + str(nitroAtom.x) + ' y => ' + str(nitroAtom.y) + ' z => ' + str(nitroAtom.z))
        print('<!> Carboxy atom  <!>: ' + str(carboxyAtom.element) + ' | x => ' + str(carboxyAtom.x) + ' y => ' + str(carboxyAtom.y) + ' z => ' + str(carboxyAtom.z)+'\n')
        print ('Energy => ' + str(energy)+'\n----------------------------------------------------------------\n\n')

    anglePhi += 1/rotationSteps
