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

    def __init__(self, psi4_path, input_file_energies_psi4, output_file_energies_psi4, energy_method, precalculated_energies_path, model_path, window_size, max_aa_length, initialization_option, n_threads, basis = 'cc-pvdz'):

        ## PARAMETERS ##

        self.model_path = model_path
        self.window_size = window_size
        self.max_aa_length = max_aa_length
        self.initialization_option = initialization_option
        self.precalculated_energies_path = precalculated_energies_path

        #Declare the instances to use the functions of these classes
        self.psi = psiFour.PsiFour(psi4_path, input_file_energies_psi4, output_file_energies_psi4, precalculated_energies_path, energy_method, n_threads, basis)
        self.tools = utils.Utils()

        #HARDCODED. It is assumed that all aminoacids has the nitro and carboxy conexions like that
        #TODO: To study if it is necessary to generalize this assumption
        self.nitroConnections = [['C', 2]]
        self.carboxyConnections = [['C', 1], ['O', 1], ['N', 1]]

    #Calculate all posible energies for the protein and the number of rotations given
    def calculate_delta_energies(self, proteinName, numberBitsRotation, method_rotations_generation, aminoacids):

        print('## Generating file of energies ##')

        #Get all atoms from the protein with x/y/z positions and connections
        atoms = self.extractAtoms(proteinName)

        #Identify the nitro and carboxy atoms
        nitroAtom = self.findAtom(atoms, 'N', '', self.nitroConnections)
        carboxyAtom = self.findAtom(atoms, '', 'Carboxy', self.carboxyConnections)

        min_energy_psi4 = self.calculateEnergyOfRotation(atoms)

        #Get initial structure of the protein to rotate from it
        [atoms, inizialitation_stats] = self.calculateInitialStructure(atoms, aminoacids, nitroAtom, carboxyAtom, method_rotations_generation)

        #Calculate all posible energies for the phi and psi angles
        deltasJson = self.calculateAllDeltasOfRotations(atoms, nitroAtom, carboxyAtom, aminoacids, min_energy_psi4, proteinName, numberBitsRotation, method_rotations_generation)

        # Add the stadistics about the precision of the inizializator
        deltasJson['inizialitation_stats'] = inizialitation_stats
        self.write_json(deltasJson, 'delta_energies', proteinName, numberBitsRotation, method_rotations_generation)

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

    def calculateInitialStructure(self, atoms, aminoacids, nitro_atom, carboxy_atom, method_rotations_generation):       
    
        #Set angles to 0. PSI4 returns the optimal angles for the protein, so it is necessary to set these angles to 0
        #Get the value of angles returned by psi4
        [atoms, phi_angles_psi4, psi_angles_psi4] = self.flat_protein(atoms, nitro_atom, carboxy_atom)

        phis_initial_rotation = []
        psis_initial_rotation = []

        #random between -π and π
        if method_rotations_generation == 'random':
            print('\n## RANDOM initialization for protein structure ##\n')

            # calculate n random angle values (n is the number of phi/psi angles that is the same than nitro/carboxy atoms)
            for _ in range(len(nitro_atom)):

                phis_initial_rotation.append(random.uniform(-math.pi, math.pi))
                psis_initial_rotation.append(random.uniform(-math.pi, math.pi))

        #minifold
        elif method_rotations_generation == 'minifold':
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
        [phis_precision, psis_precision] = self.tools.calculatePrecisionOfAngles(phi_angles_psi4, psi_angles_psi4, phis_initial_rotation, psis_initial_rotation)

        # if it is necessary convert float32 in standard python type (float32 is not serializable by json)
        if type(phis_initial_rotation[0]) is np.float32:
            phis_initial_rotation = [value.item() for value in phis_initial_rotation]
        
        if type(psis_initial_rotation[0]) is np.float32:
            psis_initial_rotation = [value.item() for value in psis_initial_rotation]

        # phis/psis initial rotation is a float 32 and it is not serializable by the json, so it is necessary to convert to a native type of python
        initilization_stats = {
            'phis_precision': phis_precision, 
            'psis_precision': psis_precision, 
            'phi_angles_psi4': phi_angles_psi4, 
            'psi_angles_psi4': psi_angles_psi4, 
            'phis_initial_rotation': phis_initial_rotation,
            'psis_initial_rotation': psis_initial_rotation
            }

        return [atoms, initilization_stats]

    #This method returns the json with all rotations and energies associated to these rotations
    def calculateAllDeltasOfRotations(self, atoms, nitroAtom, carboxyAtom, aminoacids, min_energy_psi4, proteinName, numberBitsRotation, method_rotations_generation):

        rotationSteps = pow(2, int(numberBitsRotation))
        
        # it calculates the number of necessary bits to represent the number of angles
        # example, 4 aminoacids: 3 phis/psis => 2 bits
        bits_number_angles = math.ceil(np.log2(len(aminoacids)-1))

        print('    ⬤ Calculating energies for all posible rotations')
        energies = self.calculate_all_energies(atoms, rotationSteps, 2**(len(aminoacids)-1))

        self.write_json(energies, 'energies', proteinName, numberBitsRotation, method_rotations_generation)

        #Write the headers of the energies json that is going to be returned
        deltasJson = {}
        deltasJson['protein'] = proteinName
        deltasJson['numberBitsRotation'] = numberBitsRotation
        deltasJson['psi4_min_energy'] = min_energy_psi4
        deltasJson['deltas'] = {}

        print('    ⬤ Calculating deltas for all possible combinations of rotations')

        min_energy = 99999
        index_min_energy = -1  
        # iterates over all calculated energies using the keys (contains the values of the phi/psi angles)      
        for e_key in energies.keys():

            old_energy = energies[e_key]

            # check if the energy is lower than the minimum
            if old_energy < min_energy:
                min_energy = old_energy
                index_min_energy = e_key
            
            angle_keys = e_key.split(' ')
            # iterate over all angles keys
            for index_a_key in range(len(angle_keys)):

                # calculate the plus/minus 1 rotation delta
                for plusminus in [0,1]:
                    pm = (-2)*plusminus + 1

                    new_value = (int(angle_keys[index_a_key]) + pm) % (2**numberBitsRotation)

                    # create a key with the values of the angles (if the index is equal to the index of the modified angle, insert the modified one)
                    angle_key = ''
                    binary_key = ''
                    for index_key in range(len(angle_keys)):

                        binary_key += np.binary_repr(int(angle_keys[index_key]), width = numberBitsRotation)

                        if index_key == index_a_key:
                            angle_key += str(new_value)+ ' '
                        
                        else:
                            angle_key += angle_keys[index_key] + ' '

                    new_energy = energies[angle_key.strip()]

                    # add 0/1 for phi/psi
                    if index_a_key % 2 == 0:
                        # if it is even number is phi (add 0)
                        binary_key += str(0)
                    else:
                        # if it is odd number is psi (add 1)
                        binary_key += str(1)
                    
                    # add the index of the phi/psi angle (with more than 2 aminoacids, there are more than one phi/psi)
                    binary_key += np.binary_repr(int(index_a_key/2), width = bits_number_angles)

                    # add 0/1 for plus/minus
                    binary_key += str(plusminus)

                    delta = new_energy - old_energy
                    
                    #Add the values to the file with the precalculated energies
                    deltasJson['deltas'][binary_key] = delta

        deltasJson['initial_min_energy'] = min_energy
        deltasJson['index_min_energy'] = index_min_energy.replace(' ', '')

        return deltasJson

    # RECURSIVE function to calculate all energies of each possible rotation 
    def calculate_all_energies(self, atoms, rotation_steps, protein_sequence_length, index_sequence='', energies = {}):

        # iterate to calculate all possible rotations
        # for example if there are 4 rotation steps, it executes the loop 4 times, but in each iteration, it calls recursively to all rotations starting with 0 (first iteration)  
        for index in range(rotation_steps):

            if protein_sequence_length > 0:
                # returned energy is added to a data structure (this structure is multi-dimensional)
                # index_sequence contains the accumulated index (it helps to know the general index_sequence)
                energies = self.calculate_all_energies(atoms, rotation_steps, protein_sequence_length-1, index_sequence+str(index)+' ', energies)
            
            else:
                
                #Perform the rotations over a copy
                copied_atoms = copy.deepcopy(atoms)
                copied_nitroAtoms = self.findAtom(copied_atoms, 'N', '', self.nitroConnections)
                copied_carboxyAtoms = self.findAtom(copied_atoms, '', 'Carboxy', self.carboxyConnections)

                x_values = []
                y_values = []

                # remove last whitespace
                index_sequence = index_sequence.strip()
                index_values = index_sequence.split(' ')
                for index in range(len(index_values)):

                    if index%2 == 0:
                        # rotation sequence even (0, 2, 4, ...)
                        x_values.append(int(index_values[index])) 

                    if index%2 != 0:
                        # rotation sequence odd (1, 3, 5, ...)
                        y_values.append(int(index_values[index]))

                for index in range(len(x_values)):

                    #Always rotate from state (0,0)
                    self.tools.rotate('phi', (x_values[index]/rotation_steps) * 2*math.pi, copied_nitroAtoms[index]) 
                    self.tools.rotate('psi', (y_values[index]/rotation_steps) * 2*math.pi, copied_carboxyAtoms[index])
                
                #Calculate the energy of the protein structure after the previous rotations
                energies[index_sequence] = self.calculateEnergyOfRotation(copied_atoms)
                print(index_sequence)

                # We eliminate previous copies
                del copied_atoms
                del copied_carboxyAtoms
                del copied_nitroAtoms

                break

        return energies

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

    def write_json(self, json_data, file_name, proteinName, numberBitsRotation, method_rotations_generation):

        #Create json with calculated energies
        #TODO: extract the path to a config file
        with open(self.precalculated_energies_path+file_name+'_'+proteinName+'_'+str(numberBitsRotation)+'_'+method_rotations_generation+'.json', 'w') as outfile:
            json.dump(json_data, outfile)