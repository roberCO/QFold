import utils
import psiFour
import json
import copy
import minifold
import math

class Initializer():

    def __init__(self):

        ## PARAMETERS ##

        #Declare the instances to use the functions of these classes
        self.psi = psiFour.PsiFour()
        self.tools = utils.Utils()

        #HARDCODED. It is assumed that all aminoacids has the nitro and carboxy conexions like that
        #TODO: To study if it is necessary to generalize this assumption
        self.nitroConnections = [['C', 2]]
        self.carboxyConnections = [['C', 1], ['O', 2]]

        #HARDCODED. The name of files would be in a config file
        #TODO: Extract these names to a config file
        self.inputFilenameEnergyPSI4 = 'inputRotations'
        self.outputFilenameEnergyPSI4 = 'outputRotations'

        #It contains the energy of each angle combination
        self.anglesEnergy = []

        #Initialization options (0: random, 1: minifold)
        self.initialization_option = 1

    #Calculate all posible energies for the protein and the number of rotations given
    def calculateEnergies(self, proteinName, numberBitsRotation, aminoacids):

        #Get all atoms from the protein with x/y/z positions and connections
        atoms = self.extractAtoms(proteinName)

        #Identify the nitro and carboxy atoms
        nitroAtom = self.findAtom(atoms, 'N', '', self.nitroConnections)
        carboxyAtom = self.findAtom(atoms, '', 'Carboxy', self.carboxyConnections)

        #Get initial structure of the protein to rotate from it
        predictions = self.calculateInitialStructure(atoms, aminoacids, nitroAtom, carboxyAtom)

        real_angle_phi = self.tools.calculateAngle(atoms[8], atoms[5], atoms[3], atoms[6], 'phi')
        real_angle_psi = self.tools.calculateAngle(atoms[4], atoms[7], atoms[6], atoms[3], 'psi')

        estimation_angle_phi = predictions[0][0]
        estimation_angle_psi = predictions[0][1]

        phi_precision = 0
        psi_precision = 0

        if estimation_angle_phi > real_angle_phi:

            #Calculate the distance if the angles go to zero and starts again
            option_1 = abs(math.pi - estimation_angle_phi) + abs(-math.pi - real_angle_phi)
            option_2 = abs(estimation_angle_phi - real_angle_phi)

            minimum_option = min(option_1, option_2)
            phi_precision = (1-(minimum_option / math.pi))*100

        else:

            #Calculate the distance if the angles go to zero and starts again
            option_1 = abs(math.pi - real_angle_phi) + abs(-math.pi - estimation_angle_phi)
            option_2 = abs(estimation_angle_phi - real_angle_phi)

            minimum_option = min(option_1, option_2)
            phi_precision = (1-(minimum_option / math.pi))*100

        if estimation_angle_psi > real_angle_psi:

            #Calculate the distance if the angles go to zero and starts again
            option_1 = abs(math.pi - estimation_angle_psi) + abs(-math.pi - real_angle_psi)
            option_2 = abs(estimation_angle_psi - real_angle_psi)

            minimum_option = min(option_1, option_2)
            psi_precision = (1-(minimum_option / math.pi))*100

        else:

            #Calculate the distance if the angles go to zero and starts again
            option_1 = abs(math.pi - real_angle_psi) + abs(-math.pi - estimation_angle_psi)
            option_2 = abs(estimation_angle_psi - real_angle_psi)

            minimum_option = min(option_1, option_2)
            psi_precision = (1-(minimum_option / math.pi))*100

        print('PHI precision: ', phi_precision, '%')
        print('PSI precision: ', psi_precision, '%')


        #Calculate all posible energies for the phi and psi angles
        energiesJson = self.calculateAllEnergiesOfRotations(atoms, nitroAtom, carboxyAtom, proteinName, numberBitsRotation)

        self.writeFileEnergies(energiesJson, proteinName, numberBitsRotation)

    #Get the atoms (and the properties) of a protein
    def extractAtoms(self, proteinName):

        #call psi4 to get the atoms of the protein
        atoms = self.psi.getAtomsFromProtein(proteinName)

        #Calculate the connection between atoms
        atoms = self.tools.calculateAtomConnection(atoms)

        return atoms

    #Find the atom that satisfies the parameters
    def findAtom(self, atoms, element, cType, connections):

        #Search in the whole atom list
        for at in atoms:

            #The element for the angle phi is nitrogen (N)
            if ((element != '' and at.element == element) or (cType != '' and at.c_type == cType)):

                #Check all the atom connections
                for conn in connections:

                    linkedElement = conn[0]
                    numberLinkedElements = conn[1]
                    counterLinkedElements = 0

                    #The valid nitrogen is the connected with two carbons
                    for elementConn in at.linked_to:

                        if(elementConn.element == linkedElement):
                            counterLinkedElements += 1

                    #Atom found and saved in variable
                    if (counterLinkedElements == numberLinkedElements):
                        return at

                    else:
                        raise Exception('Element '+at.element+' not found with the proper connections of '+linkedElement+'! '+str(counterLinkedElements)+' found but there should be ' + str(numberLinkedElements))

    def calculateInitialStructure(self, atoms, aminoacids, nitro_atom, carboxy_atom):       
    
        #Set angles to 0. PSI4 returns the optimal angles for the protein, so it is necessary to set these angles to 0
        #Get the value of angles returned by psi4
        angle_phi = self.tools.calculateAngle(atoms[8], atoms[5], atoms[3], atoms[6], 'phi')
        angle_psi = self.tools.calculateAngle(atoms[4], atoms[7], atoms[6], atoms[3], 'psi')


        #Rotate the inverse angles of psi4 to get angles to 0
        self.tools.rotate('phi', angle_phi, nitro_atom) 
        #self.tools.rotate('psi', -angle_psi, carboxy_atom)

        print('Angle phi:\nbefore moved to 0:', angle_phi,'\nmovement:', angle_phi,'\nafter moved to 0:', self.tools.calculateAngle(atoms[8], atoms[5], atoms[3], atoms[6], 'phi'), '\n')
        #print('Angle psi:\nbefore moved to 0:', angle_psi,'\nmovement:', -angle_psi,'\nafter moved to 0:', self.tools.calculateAngle(atoms[4], atoms[7], atoms[6], atoms[3], 'psi'), '\n')

        #Apply the calculated rotations for the angles

        #TODO: random
        if self.initialization_option == 0:
            print('Random angle inizialition')

        #minifold
        elif self.initialization_option == 1:

            mfold = minifold.Minifold()
            angles = mfold.predictAngles(aminoacids)


        return angles

    #This method returns the json with all rotations and energies associated to these rotations
    def calculateAllEnergiesOfRotations(self, atoms, nitroAtom, carboxyAtom, proteinName, numberBitsRotation):

        rotationSteps = pow(2, int(numberBitsRotation))

        #Get the size of each step to rotate
        anglePhi = 1/rotationSteps
        anglePsi = 1/rotationSteps

        #Write the headers of the energies json that is going to be returned
        energiesJson = {}
        energiesJson['protein'] = proteinName
        energiesJson['numberBitsRotation'] = numberBitsRotation
        energiesJson['energies'] = []

        #These two nested loops are hardcoded (it could be n nested loops, 1 per AA) because QFold is going to be used just with two and three aminoacids
        #if it is scale to more aminoacids, it should be necessary to implement a recursive function
        #TODO: To study how to generalize to any number of aminoacids
        for x in range(0, rotationSteps):

            for y in range(0, rotationSteps):

                #Perform the rotations over a copy
                copied_atoms = copy.deepcopy(atoms)
                copied_nitroAtom = self.findAtom(copied_atoms, 'N', '', self.nitroConnections)
                copied_carboxyAtom = self.findAtom(copied_atoms, '', 'Carboxy', self.carboxyConnections)

                #Always rotate from state (0,0)
                self.tools.rotate('phi', x * anglePhi, copied_nitroAtom) 
                self.tools.rotate('psi', y * anglePsi, copied_carboxyAtom)
                
                #Calculate the energy of the protein structure after the previous rotations
                energy = self.calculateEnergyOfRotation(copied_atoms)

                #Normalize the energy
                normalizedEnergy = energy*-1
                normalizedEnergy = "%.6f" % normalizedEnergy

                #Add the angles to list with angles and energies
                self.anglesEnergy.append([anglePhi, anglePsi, normalizedEnergy])

                #Add the values to the file with the precalculated energies
                energiesJson['energies'].append({
                    'phi': x,
                    'psi': y,
                    'energy': energy,
                })

                # We eliminate previous copies
                del copied_atoms
                del copied_carboxyAtom
                del copied_nitroAtom

                #Increment the angle psi counter INSIDE the SECOND loop
                anglePsi += 1/rotationSteps

            #Increment the angle psi counter INSIDE the FIRST loop (OUTSIDE the SECOND)
            anglePhi += 1/rotationSteps

        return energiesJson

    def calculateEnergyOfRotation(self, copied_atoms):

        #Write the file with the actual rotations
        self.psi.writeFileEnergies(copied_atoms, self.inputFilenameEnergyPSI4)

        #Calculate the energy of the actual rotations using PSI4
        self.psi.executePsiCommand(self.inputFilenameEnergyPSI4, self.outputFilenameEnergyPSI4)

        #Read the PSI4 output file and get the energy
        energy = self.psi.readEnergyFromFile(self.outputFilenameEnergyPSI4)

        return energy

    def writeFileEnergies(self, energiesJson, proteinName, numberBitsRotation):

        #Create json with calculated energies
        #HARDCODED path
        #TODO: extract the path to a config file
        with open('./precalculated_energies/energies_'+proteinName+'_'+str(numberBitsRotation)+'.json', 'w') as outfile:
            json.dump(energiesJson, outfile)
