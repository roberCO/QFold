import sys
import utils
import psiFour
import atom
import numpy as np

def distance(atom,atom2):
    return np.sqrt((atom.x-atom2.x)**2+(atom.y-atom2.y)**2+(atom.z-atom2.z)**2)

def calculateAtomConnection(atoms):

    for at1 in atoms:
        for at2 in atoms:
            if at1 != at2:
                if at1.element == 'O' and at2.element == 'C' and distance(at1,at2)<2:
                    at2.c_type = 'Carboxy'
                    at1.linked_to = [at2] + at1.linked_to 
                    at2.linked_to = [at1] + at2.linked_to
                    
                if at1.element == 'N' and at2.element == 'C' and distance(at1,at2)<2:
                    if at2.c_type != 'Carboxy':
                        at2.c_type = 'C_alpha'
                    at1.linked_to = [at2] + at1.linked_to 
                    at2.linked_to = [at1] + at2.linked_to
                    
                if at1.element == 'H' and at2.element == 'C' and distance(at1,at2)<1.3: 
                    at1.linked_to = [at2] + at1.linked_to 
                    at2.linked_to = [at1] + at2.linked_to
                    
                if at1.element == 'C' and at2.element == 'C' and distance(at1,at2)<2  and (at1 not in at2.linked_to) : 
                    at1.linked_to = [at2] + at1.linked_to 
                    at2.linked_to = [at1] + at2.linked_to 

    return atoms



if(len(sys.argv) < 2):
    print ("<*> ERROR: You must specify the number of aminoacids - Usage: python main.py numberOfAminoacids ProteinName [Aminoacid1] [Aminoacid2] [...]")
    print ("<!> Example: python main.py 2 Glycylglycine GLY GLY")
    sys.exit(0)

try:
    numberAA = int(sys.argv[1])
except ValueError:
    print ("<*> ERROR: The number of aminoacids must be an integer and you enter: "+sys.argv[1]+" - Usage: python main.py numberOfAminoacids ProteinName [Aminoacid1] [Aminoacid2] [...]")
    print ("<!> Example: python main.py 2 Glycylglycine GLY GLY")
    sys.exit(0)

if(len(sys.argv) != numberAA+3):
    print ("<*> ERROR: Incorrect number of aminoacids. You specify " + str(numberAA) + " but you enter "+ str(len(sys.argv)-3) +" - Usage: python main.py numberOfAminoacids ProteinName [Aminoacid1] [Aminoacid2] [...]")
    print ("<!> Example: python main.py 2 Glycylglycine GLY GLY")
    sys.exit(0)


proteinName = sys.argv[2]
aminoacids = []

for index in range(0, numberAA):
    aminoacids.append(sys.argv[index+3])


#call psi4 to get the atoms of the protein
psi = psiFour.PsiFour()
atoms = psi.getAtomsFromProtein(proteinName)

#Calculate the connection between atoms
atoms = calculateAtomConnection(atoms)

tools = utils.Utils()

'''
for at in atoms:

    for conn in at.linked_to:

        print("atom: " + at.element + " linked to: " + conn.element + " atom type: " + at.c_type)

    print("------")

#Initial call minifold
initialAngleConfiguration = tools.generateInitialConfig(aminoacids)

#reache initial configuration of angles
atoms = tools.rotateAminoacid(atoms, initialAngleConfiguration)
'''


#HARDCODED
rotationStep = 1/32 #It will multipled by pi (64 steps)


#This is hardcoded because QFold is going to be used just with two and three aminoacids
#if it is scale to more aminoacids, it should be necessary to implement a recursive function

nitroAtom = None
carboxyAtom = None

for at in atoms:

    #The element for the angle phi is nitrogen (N)
    if (at.element == 'N'):

        linkedCarbons = 0
        #The valid nitrogen is the connected with two carbons
        for conn in at.linked_to:

            if(conn.element == 'C'):
                linkedCarbons += 1

        #Atom found and saved in variable
        if (linkedCarbons == 2):
            nitroAtom = at
            #Remove the element that is going to be modified, then the modified element will be added
            atoms.remove(at)



for at in atoms:

    #The element for the angle psi is a carboxy
    if (at.c_type == 'Carboxy'):

        linkedCarbons = 0
        linkedOxygens = 0

        #The valid carboxy is the connected with one carbon (C) and two oxygens
        for conn in at.linked_to:

            if(conn.element == 'C'):
                linkedCarbons += 1

            if(conn.element == 'O'):
                linkedOxygens += 1

        #Atom found and saved in a variable
        if(linkedCarbons == 1 and linkedOxygens == 2):
            carboxyAtom = at
            #Remove the element that is going to be modified, then the modified element will be added
            atoms.remove(at)

anglePhi = 0
for x in range(0, 64):

    tools.rotate('phi', anglePhi, nitroAtom)
    anglePsi = 0

    for y in range(0, 64):

        tools.rotate('psi', anglePsi, carboxyAtom)
        anglePsi += rotationStep

        #Write file with all atoms rotated
        rotationHandle = open('inputRotations.dat', 'w')

        rotationHandle.write('molecule glycylglycine{\n')
        #write input.dat with all rotated atoms
        for at in atoms:
            rotationHandle.write(" " + at.element + " " + str(at.x) + " " + str(at.y) + " " + str(at.z)+'\n')
        
        rotationHandle.write(" " + nitroAtom.element + " " + str(nitroAtom.x) + " " + str(nitroAtom.y) + " " + str(nitroAtom.z)+'\n')
        rotationHandle.write(" " + carboxyAtom.element + " " + str(carboxyAtom.x) + " " + str(carboxyAtom.y) + " " + str(carboxyAtom.z)+'\n')
        rotationHandle.write('}\n\n')
        rotationHandle.write("set basis cc-pvdz\n")
        rotationHandle.write("set reference rhf\n")
        rotationHandle.write("energy('scf')\n")

        rotationHandle.close()

        psi.executePsiCommand('inputRotations', 'energyRotations')

        with open('energyRotations.dat', 'r') as fileHandle:
            for line in fileHandle:
                if 'Final Energy' in line:
                    energy = float(line.split(':')[1])
                    print('energy: ' + str(energy))

    anglePhi += rotationStep

