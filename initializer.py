import utils
import psiFour
import json
import copy
import minifold
import math
import random
import sys
import progressbar
import numpy as np

class Initializer():

    def __init__(self, psi4_path, input_file_energies_psi4, output_file_energies_psi4, energy_method, precalculated_energies_path, model_path, window_size, max_aa_length, initialization_option, basis = 'cc-pvdz'):

        ## PARAMETERS ##

        self.model_path = model_path
        self.window_size = window_size
        self.max_aa_length = max_aa_length
        self.initialization_option = initialization_option
        self.precalculated_energies_path = precalculated_energies_path

        #Declare the instances to use the functions of these classes
        self.psi = psiFour.PsiFour(psi4_path, input_file_energies_psi4, output_file_energies_psi4, precalculated_energies_path, energy_method, basis)
        self.tools = utils.Utils()

        #HARDCODED. It is assumed that all aminoacids has the nitro and carboxy conexions like that
        #TODO: To study if it is necessary to generalize this assumption
        self.nitroConnections = [['C', 2]]
        self.carboxyConnections = [['C', 1], ['O', 1], ['N', 1]]

    #Calculate all posible energies for the protein and the number of rotations given
    def calculate_delta_energies(self, proteinName, numberBitsRotation, aminoacids):

        print('## Generating file of energies ##')

        #Get all atoms from the protein with x/y/z positions and connections
        atoms = self.extractAtoms(proteinName)

        #Identify the nitro and carboxy atoms
        nitroAtom = self.findAtom(atoms, 'N', '', self.nitroConnections)
        carboxyAtom = self.findAtom(atoms, '', 'Carboxy', self.carboxyConnections)

        min_energy_psi4 = self.calculateEnergyOfRotation(atoms)

        #Get initial structure of the protein to rotate from it
        atoms = self.calculateInitialStructure(atoms, aminoacids, nitroAtom, carboxyAtom)

        #Calculate all posible energies for the phi and psi angles
        energiesJson = self.calculateAllDeltasOfRotations(atoms, nitroAtom, carboxyAtom, min_energy_psi4, proteinName, numberBitsRotation)

        self.writeFileEnergies(energiesJson, proteinName, numberBitsRotation)

    #Get the atoms (and the properties) of a protein
    def extractAtoms(self, proteinName):

        print('    ⬤ Extracting atoms from proteins')
        #call psi4 to get the atoms of the protein
        atoms = self.psi.getAtomsFromProtein(proteinName)

        print('    ⬤ Calculating connections between atoms')
        #Calculate the connection between atoms
        atoms = self.tools.calculateAtomConnection(atoms)

        return atoms

    #Find the atom that satisfies the parameters
    def findAtom(self, atoms, element, cType, connections):

        found_atoms = []

        # Search in the whole atom list
        for at in atoms:

            # Check if the element or the cType is the same that searched element 
            if ((element != '' and at.element == element) or (cType != '' and at.c_type == cType)):

                # at (atom) has the same element/c_type than the searched atom, now it is necessary to check its connections
                found = True

                #Check all the searched atom connections
                for conn in connections:
                    
                    linked_element = conn[0]
                    number_linked_element = conn[1]
                    counter_at_linked_elements = 0

                    #The valid nitrogen is the connected with two carbons
                    for elementConn in at.linked_to:

                        if(elementConn.element == linked_element):
                            counter_at_linked_elements += 1
                            
                    # It means that the at has different connections that expected
                    if(number_linked_element != counter_at_linked_elements):
                        found = False
                        break

                # If all connections of the atom were found, the atom is added to the found atoms list
                if (found):
                    found_atoms.append(at)

        if len(found_atoms) == 0:
            raise Exception('Element '+element+' not found with the proper connections of '+connections+'!')

        return found_atoms

    def calculateInitialStructure(self, atoms, aminoacids, nitro_atom, carboxy_atom):       
    
        #Set angles to 0. PSI4 returns the optimal angles for the protein, so it is necessary to set these angles to 0
        #Get the value of angles returned by psi4
        [atoms, phi_angles_psi4, psi_angles_psi4] = self.flat_protein(atoms, nitro_atom, carboxy_atom)

        phis_initial_rotation = []
        psis_initial_rotation = []

        #random between -π and π
        if self.initialization_option == 'random':
            print('\n## RANDOM initialization for protein structure ##\n')

            # calculate n random angle values (n is the number of phi/psi angles that is the same than nitro/carboxy atoms)
            for _ in range(len(nitro_atom)):

                phis_initial_rotation.append(random.uniform(-math.pi, math.pi))
                psis_initial_rotation.append(random.uniform(-math.pi, math.pi))

        #minifold
        elif self.initialization_option == 'minifold':
            print('\n## MINIFOLD initialization for protein structure ##\n')
            mfold = minifold.Minifold(self.model_path, self.window_size, self.max_aa_length)
            angles = mfold.predictAngles(aminoacids)

            for angle in angles:

                phis_initial_rotation.append(angle[0])
                psis_initial_rotation.append(angle[1])

        #Rotate all angles to get the initial protein structure
        for index in range(len(phis_initial_rotation)):

                self.tools.rotate('psi', psis_initial_rotation[index], carboxy_atom[index])
                self.tools.rotate('phi', phis_initial_rotation[index], nitro_atom[index]) 


        #Calculate the precision in constrast of the real value calculated by psi4
        self.tools.calculatePrecisionOfAngles(phi_angles_psi4, psi_angles_psi4, phis_initial_rotation, psis_initial_rotation)
        return atoms

    #This method returns the json with all rotations and energies associated to these rotations
    def calculateAllDeltasOfRotations(self, atoms, nitroAtom, carboxyAtom, min_energy_psi4, proteinName, numberBitsRotation):

        rotationSteps = pow(2, int(numberBitsRotation))

        [energies, min_energy, phi_position_min_energy, psi_position_min_energy] = self.calculate_all_energies(rotationSteps, atoms)

        #Write the headers of the energies json that is going to be returned
        energiesJson = {}
        energiesJson['protein'] = proteinName
        energiesJson['numberBitsRotation'] = numberBitsRotation
        energiesJson['psi4_min_energy'] = min_energy_psi4
        energiesJson['initial_min_energy'] = min_energy
        energiesJson['phi_position_min_energy'] = phi_position_min_energy
        energiesJson['psi_position_min_energy'] = psi_position_min_energy
        energiesJson['deltas'] = {}

        print('    ⬤ Calculating deltas for all posible combinations of rotations')
        for index_phi in range(len(energies)):
            for index_psi in range(len(energies[index_phi])):
            
                old_energy = energies[index_phi][index_psi]
            
                for plusminus in [0,1]:
                    pm = 2*plusminus - 1
                    
                    for phipsi in [0,1]:
                        if phipsi == 0:
                            new_phi = (index_phi + pm) % (2**numberBitsRotation)
                            new_psi = index_psi
                            
                        if phipsi == 1:
                            new_phi = index_phi
                            new_psi = (index_psi + pm) % (2**numberBitsRotation)
                            
                        new_energy = energies[new_phi][new_psi]
                        
                        binary_phi = self.tools.angle_to_binary(new_phi, numberBitsRotation)
                        binary_psi = self.tools.angle_to_binary(new_psi, numberBitsRotation)

                        angle_key = binary_phi + binary_psi + str(phipsi) + str(plusminus) 
                        
                        delta = new_energy - old_energy

                        #Add the values to the file with the precalculated energies
                        energiesJson['deltas'][angle_key] = delta

        return energiesJson

    def calculate_all_energies(self, rotationSteps, atoms):

        #These two nested loops are hardcoded (it could be n nested loops, 1 per AA) because QFold is going to be used just with two and three aminoacids
        #if it is scale to more aminoacids, it should be necessary to implement a recursive function
        #TODO: To study how to generalize to any number of aminoacids

        energies = [[0 for x in range(rotationSteps)] for y in range(rotationSteps)]
        anglePhi = 1/rotationSteps
        anglePsi = 1/rotationSteps

        phi_position_min_energy = -1
        psi_position_min_energy = -1
        min_energy = -1

        print('    ⬤ Calculating energies for all posible rotations')
        bar = progressbar.ProgressBar(maxval=2**rotationSteps, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        for x in range(0, rotationSteps):

            for y in range(0, rotationSteps):

                #Perform the rotations over a copy
                copied_atoms = copy.deepcopy(atoms)
                copied_nitroAtom = self.findAtom(copied_atoms, 'N', '', self.nitroConnections)
                copied_carboxyAtom = self.findAtom(copied_atoms, '', 'Carboxy', self.carboxyConnections)

                #Always rotate from state (0,0)
                self.tools.rotate('phi', x * 2*math.pi * anglePhi, copied_nitroAtom) 
                self.tools.rotate('psi', y * 2*math.pi * anglePsi, copied_carboxyAtom)
                
                #Calculate the energy of the protein structure after the previous rotations
                energy = self.calculateEnergyOfRotation(copied_atoms)
                energies[x][y] = energy 

                if energy < min_energy:
                    phi_position_min_energy = x
                    psi_position_min_energy = y
                    min_energy = energy

                # We eliminate previous copies
                del copied_atoms
                del copied_carboxyAtom
                del copied_nitroAtom

                #Increment the angle psi counter INSIDE the SECOND loop
                anglePsi += 1/rotationSteps

                bar.update(x*rotationSteps + y)

            #Increment the angle psi counter INSIDE the FIRST loop (OUTSIDE the SECOND)
            anglePhi += 1/rotationSteps

        bar.finish()

        return [energies, min_energy, phi_position_min_energy, psi_position_min_energy]

    def calculateEnergyOfRotation(self, copied_atoms):

        #Write the file with the actual rotations
        self.psi.writeFileEnergies(copied_atoms)

        #Calculate the energy of the actual rotations using PSI4
        self.psi.executePsiCommand()

        #Read the PSI4 output file and get the energy
        energy = self.psi.readEnergyFromFile()

        return energy

    def flat_protein(self, atoms, nitro_atoms, carboxy_atoms):

        phi_angles_psi4 = []
        psi_angles_psi4 = []

        next_atom = self.get_initial_atom(atoms)

        all_psi = self.get_all_angle_planes(atoms, carboxy_atoms, 'psi')
        all_phi = self.get_all_angle_planes(atoms, nitro_atoms, 'phi')

        end = False
        while(not end):

            end = True
            found = False
            for psi in all_psi:

                for atom in psi:
                    if next_atom.atomId == atom.atomId:
                        found = True
                        break

                if found:

                    end = False

                    angle = self.tools.calculateAngle(psi, 'psi')
                    psi_angles_psi4.append(angle)

                    #Rotate the inverse (*-1) angles of psi4 to get angles to 0
                    c_atom = psi[2] # psi in 2 contains the main atom of the plane (carboxy atom)
                    self.tools.rotate('psi', -angle, c_atom)

                    next_atom = psi[len(psi)-1] 

                    break

            # the end of the chain has been reached
            if end:
                break
            
            found = False
            for phi in all_phi:

                for atom in phi:
                    if next_atom.atomId == atom.atomId:
                        found = True
                        break

                if found: 

                    angle = self.tools.calculateAngle(phi, 'phi')
                    phi_angles_psi4.append(angle)

                    #Rotate the inverse (*-1) angles of psi4 to get angles to 0
                    n_atom = phi[2] # phi in 2 contains the main atom of the plane (nitro atom)
                    self.tools.rotate('phi', -angle, n_atom)

                    next_atom = phi[0] 

                    break
            

        return [atoms, phi_angles_psi4, psi_angles_psi4]

    def get_initial_atom(self, atoms):

        # get initial point to start to flat the protein (initial point is H-N-H)
        initial_atom = atoms[0]
        for at in atoms:
            if at.element == 'N':

                hidrogen_counter = 0
                for conn in at.linked_to:
                    if conn.element == 'H':
                        hidrogen_counter += 1

                if hidrogen_counter >= 2:
                    initial_atom = at
                    break

        return initial_atom

    def get_all_angle_planes(self, atoms, angle_atoms, type_angle):

        all_angle_planes = []

        if type_angle == 'psi':

            atom_1 = 'N'
            atom_2 = 'C_alpha'
            # atom_3 is the at atom
            atom_4 = 'N'

        elif type_angle == 'phi':

            atom_1 = 'Carboxy'
            atom_2 = 'C_alpha'
            # atom_3 is the at atom
            atom_4 = 'Carboxy'

        for atom in angle_atoms:

            all_angle_atoms = []

            #The order in which atoms are added is necessary to calcule correctly the angle
            for at in atom.linked_to:
                if at.element == atom_2 or at.c_type == atom_2:

                    for at2 in at.linked_to:
                        if at2.element == atom_1 or at2.c_type == atom_1:
                            # add nitro
                            all_angle_atoms.append(at2)

                    # add c_alpha
                    all_angle_atoms.append(at)

            # add carboxy
            all_angle_atoms.append(atom)

            # add nitro
            for at in atom.linked_to:
                if at.element == atom_4 or at.c_type == atom_4:
                    all_angle_atoms.append(at)

            all_angle_planes.append(all_angle_atoms)

        return all_angle_planes

    def writeFileEnergies(self, energiesJson, proteinName, numberBitsRotation):

        #Create json with calculated energies
        #TODO: extract the path to a config file
        with open(self.precalculated_energies_path+'energies_'+proteinName+'_'+str(numberBitsRotation)+'.json', 'w') as outfile:
            json.dump(energiesJson, outfile)
