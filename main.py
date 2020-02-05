import sys
import angleCalculator
import psiFour

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


'''
aC = angleCalculator.AngleCalculator()

#Initial call minifold
initialAngleConfiguration = aC.generateInitialConfig(aminoacids)

#reache initial configuration of angles
atoms = aC.rotateAminoacid(atoms, initialAngleConfiguration)

'''

#HARDCODED
rotationStep = 10


#This is hardcoded because QFold is going to be used just with two and three aminoacids
#if it is scale to more aminoacids, it should be necessary to implement a recursive function



#rotate x axis
for x in range (0, 360, roationStep):

    #reset to initial state to start again the whole rotations cyle
    rotatedAtoms = atoms
    
    #rotate y axis
    for y in range (0, 360, rotationStep):

        aC.rotateAminoacid(aminoacids[1], 30)

        #This rotation is only executed with three AAs
        if(len(aminoacids) == 3)
            #rotate z axis
            for z in range(0, 360, rotationStep):

                rotateAtoms = aC.rotateAminoacid(rotateAtoms, z)


