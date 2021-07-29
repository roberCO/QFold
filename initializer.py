import utils
import psiFour
import json
import copy
import minifold
import math
import random
import sys
import numpy as np
import time
import datetime

class Initializer():

    def __init__(self, tools):

        ## PARAMETERS ##
        self.model_path = tools.config_variables['model_path']
        self.window_size = tools.config_variables['window_size']
        self.max_aa_length = tools.config_variables['maximum_aminoacid_length']
        self.initialization_option = tools.config_variables['methods_initialization']
        self.precalculated_energies_path = tools.config_variables['precalculated_energies_path']
        self.basis = tools.config_variables['basis']

        #Declare the instances to use the functions of these classes
        self.psi = psiFour.PsiFour(tools)
        self.tools = utils.Utils()

    #Calculate all posible energies for the protein and the number of rotations given
    def calculate_energies(self, proteinName, numberBitsRotation, method_rotations_generation, aminoacids, protein_id):

        print('## Generating file of energies ##')

        #Get all atoms from the protein with x/y/z positions and connections
        atoms, backbone = self.extractAtoms(proteinName, aminoacids, protein_id)

        # calculate the energy of the molecule returned by PSI4 (in theory the minimum energy of the molecule)
        min_energy_psi4 = self.calculateEnergyOfRotation(atoms)

        #Get initial structure of the protein to rotate from it
        [atoms, initialization_stats] = self.calculateInitialStructure(atoms, aminoacids, method_rotations_generation, backbone)

        rotationSteps = pow(2, int(numberBitsRotation))
        number_angles = 2*(len(aminoacids)-1)


        #Calculate all posible energies for the phi and psi angles
        energies_json = {}
        [energies_json['energies'], energies_json['index_min_energy']] = self.calculate_all_energies(atoms, rotationSteps, number_angles, number_angles, aminoacids)
        energies_json['initialization_stats'] = initialization_stats
        energies_json['min_energy_psi4'] = min_energy_psi4

        # Add the stadistics about the precision of the initializator
        self.write_json(energies_json, 'energies', proteinName, numberBitsRotation, method_rotations_generation)

    #Get the atoms (and the properties) of a protein
    def extractAtoms(self, proteinName, aminoacids, protein_id):

        print('    ⬤ Extracting atoms from proteins')
        #call psi4 to get the atoms of the protein
        atoms = self.psi.getAtomsFromProtein(proteinName, protein_id)
        
        # if atoms length is 0 means that the proteins was not find in the database
        if len(atoms) == 0:
            raise Exception("Protein name not found. There is no atoms for that protein")

        print('    ⬤ Calculating connections between atoms')
        #Calculate the connection between atoms
        atoms, backbone = self.tools.calculateAtomConnection(atoms,aminoacids)

        return atoms, backbone

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

    def calculateInitialStructure(self, atoms, aminoacids, method_rotations_generation, backbone):

        phis_initial_rotation = []
        psis_initial_rotation = []

        # First we calculate all the angles. Psi uses the first atom from the next aminoacid, whereas phi uses the last from the previous
        psi_angles_psi4 = [self.tools.calculateAngle(backbone[3*j:3*j+4],'psi') for j in range(len(backbone)//3 - 1)]
        phi_angles_psi4 = [self.tools.calculateAngle(backbone[3*j-1:3*j+3],'phi') for j in range(1, len(backbone)//3)]
        
        # Hardware does not modify it
        if method_rotations_generation == 'original':

            phis_initial_rotation = copy.deepcopy(phi_angles_psi4)
            psis_initial_rotation = copy.deepcopy(psi_angles_psi4)

        #random between -π and π
        elif method_rotations_generation == 'random':
            print('\n## RANDOM initialization for protein structure ##\n')

            #Set angles to 0. PSI4 returns the optimal angles for the protein, so it is necessary to set these angles to 0
            atoms = self.flat_protein(atoms, backbone, phi_angles_psi4, psi_angles_psi4)

            # calculate n random angle values (n is the number of phi/psi angles that is the same than nitro/carboxy atoms)
            for _ in range(len(phi_angles_psi4)):

                phis_initial_rotation.append(random.uniform(-math.pi, math.pi))
                psis_initial_rotation.append(random.uniform(-math.pi, math.pi))

                print('Angles', phis_initial_rotation,psis_initial_rotation)

        #minifold
        elif method_rotations_generation == 'minifold':
            print('\n## MINIFOLD initialization for protein structure ##\n')

            #Set angles to 0. PSI4 returns the optimal angles for the protein, so it is necessary to set these angles to 0
            atoms = self.flat_protein(atoms, backbone, phi_angles_psi4, psi_angles_psi4)

            mfold = minifold.Minifold(self.model_path, self.window_size, self.max_aa_length)
            angles = mfold.predictAngles(aminoacids)

            for angle in angles:

                phis_initial_rotation.append(angle[0])
                psis_initial_rotation.append(angle[1])

        #Rotate all angles to get the initial protein structure
        if method_rotations_generation != 'original':

            for index in range(len(phis_initial_rotation)):

                self.tools.rotate(angle_type = 'psi', angle = psis_initial_rotation[index], starting_atom = backbone[3*index+2], backbone = backbone)
                self.tools.rotate(angle_type = 'phi', angle = phis_initial_rotation[index], starting_atom = backbone[3*index+4], backbone = backbone) 

        #Calculate the precision in constrast of the real value calculated by psi4
        [phis_precision, psis_precision] = self.tools.calculatePrecisionOfAngles(phi_angles_psi4, psi_angles_psi4, phis_initial_rotation, psis_initial_rotation)

        # if it is necessary convert float32 in standard python type (float32 is not serializable by json)
        if type(phis_initial_rotation[0]) is np.float32:
            phis_initial_rotation = [value.item() for value in phis_initial_rotation]
        
        if type(psis_initial_rotation[0]) is np.float32:
            psis_initial_rotation = [value.item() for value in psis_initial_rotation]

        # phis/psis initial rotation is a float 32 and it is not serializable by the json, so it is necessary to convert to a native type of python
        initialization_stats = {
            'phis_precision': phis_precision, 
            'psis_precision': psis_precision, 
            'phi_angles_psi4': phi_angles_psi4, 
            'psi_angles_psi4': psi_angles_psi4, 
            'phis_initial_rotation': phis_initial_rotation,
            'psis_initial_rotation': psis_initial_rotation
            }

        return [atoms, initialization_stats]

    # RECURSIVE function to calculate all energies of each possible rotation 
    def calculate_all_energies(self, atoms, rotation_steps, protein_sequence_length, max_lenght, aminoacids, index_sequence='', index_min_energy=-1, energies = {}):

        # iterate to calculate all possible rotations
        # for example if there are 4 rotation steps, it executes the loop 4 times, but in each iteration, it calls recursively to all rotations starting with 0 (first iteration)  
        for index in range(rotation_steps):

            if max_lenght == protein_sequence_length:
                start_time = time.time()

            if protein_sequence_length > 0:
                # returned energy is added to a data structure (this structure is multi-dimensional)
                # index_sequence contains the accumulated index (it helps to know the general index_sequence)
                [energies, index_min_energy] = self.calculate_all_energies(atoms, rotation_steps, protein_sequence_length-1, max_lenght, aminoacids, index_sequence+str(index)+' ', index_min_energy, energies)

            else:
                
                #Perform the rotations over a copy
                copied_atoms = copy.deepcopy(atoms)
                for at in copied_atoms:
                    if at.c_type == 'N_backbone' and ((len(at.linked_to_dict['C']) == 1 and len(at.linked_to_dict['H']) == 2) or self.tools.is_proline_N(at)):
                        nitro_start = at
                        break

                copied_backbone = self.tools.main_chain_builder([nitro_start], aminoacids)

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

                    #Always rotate from state (0,0) angle_type, angle, starting_atom, backbone
                    self.tools.rotate(angle_type='psi', angle=(y_values[index]/rotation_steps) * 2*math.pi, starting_atom = copied_backbone[3*index + 2], backbone = copied_backbone)
                    self.tools.rotate(angle_type='phi', angle=(x_values[index]/rotation_steps) * 2*math.pi, starting_atom = copied_backbone[3*index + 4], backbone = copied_backbone) 
                    
                
                #Calculate the energy of the protein structure after the previous rotations
                energies[index_sequence] = self.calculateEnergyOfRotation(copied_atoms)

                if index_min_energy == -1 or energies[index_sequence] < energies[index_min_energy]: index_min_energy = index_sequence

                # We eliminate previous copies
                del copied_atoms
                del copied_backbone

                break

            if max_lenght == protein_sequence_length:
                total_time = time.time() - start_time
                print("Step", index+1, "of", rotation_steps,"calculated for aminoacids", aminoacids, "in", str(datetime.timedelta(seconds=total_time)), " in hh:mm:ss (", str(datetime.timedelta(seconds=total_time/(rotation_steps**(max_lenght-1)))), "per each)")

        return [energies, index_min_energy]

    def calculateEnergyOfRotation(self, copied_atoms):

        #Write the file with the actual rotations
        self.psi.writeFileEnergies(copied_atoms)

        #Calculate the energy of the actual rotations using PSI4
        self.psi.executePsiCommand()

        #Read the PSI4 output file and get the energy
        energy = self.psi.readEnergyFromFile()

        return energy

    def flat_protein(self, atoms, backbone, phi_angles_psi4, psi_angles_psi4):

        # Next we need to flatten the peptide
        for i in range(len(psi_angles_psi4)):
            # For psi we have to rotate -angle starting in the carboxy of the i-th aminoacid
            self.tools.rotate(angle_type = 'psi', angle = -1*psi_angles_psi4[i], starting_atom = backbone[3*i+2], backbone = backbone)
            # For phi we have to rotate -angle starting in the C_alpha of the (i+1)-th aminoacid
            self.tools.rotate(angle_type ='phi', angle = -1*phi_angles_psi4[i], starting_atom = backbone[3*i+4], backbone = backbone)

        #zeros = [self.tools.calculateAngle(backbone[3*j:3*j+4],'psi') for j in range(len(backbone)//3 - 1)]
        #zeros += [self.tools.calculateAngle(backbone[3*j-1:3*j+3],'phi') for j in range(1, len(backbone)//3)]

        #self.tools.plotting(list_of_atoms = atoms, title = 'Peptide_plot_flattened')
        return atoms

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
    
    def get_energy_configuration_from_position(self, position, initial_args):

        energy = 0

        # calculate the structure (energy and configuration) of the protein from the position calculated by metropolis algorithms
        # it is possible to know the protein structure because it has the initial position and how many degrees was rotated (position * number of rotation bits)
        # First we calculate all the angles. Psi uses the first atom from the next aminoacid, whereas phi uses the last from the previous        

        # first half of position string is phi positions and the other half is psi positions
        phi_positions = position[:int(len(position)/2)]
        psi_positions = position[int(len(position)/2):]

        # get atoms
        atoms = self.psi.getAtomsFromProtein(initial_args.protein_name, initial_args.id)
        atoms, backbone = self.tools.calculateAtomConnection(atoms, initial_args.aminoacids)

        atoms = self.calculate_structure(atoms, initial_args.aminoacids, initial_args.initialization, initial_args.bits, backbone, phi_positions, psi_positions)
        energy = self.calculateEnergyOfRotation(atoms)

        configuration = self.convert_atoms_to_configuration(atoms)

        return [energy, configuration]

    def calculate_structure(self, atoms, aminoacids, init_method, bits, backbone, phi_positions, psi_positions):

        phis_initial_rotation = []
        psis_initial_rotation = []
        rotation_steps = pow(2, int(bits))

        psi_angles_psi4 = [self.tools.calculateAngle(backbone[3*j:3*j+4],'psi') for j in range(len(backbone)//3 - 1)]
        phi_angles_psi4 = [self.tools.calculateAngle(backbone[3*j-1:3*j+3],'phi') for j in range(1, len(backbone)//3)]

        
        #Set angles to 0. PSI4 returns the optimal angles for the protein, so it is necessary to set these angles to 0
        atoms = self.flat_protein(atoms, backbone, phi_angles_psi4, psi_angles_psi4)
        
        #random between -π and π
        if init_method == 'random':
            for _ in range(len(phi_angles_psi4)):

                phis_initial_rotation.append(random.uniform(-math.pi, math.pi))
                psis_initial_rotation.append(random.uniform(-math.pi, math.pi))

        #minifold
        elif init_method == 'minifold':

            mfold = minifold.Minifold(self.model_path, self.window_size, self.max_aa_length)
            angles = mfold.predictAngles(aminoacids)

            for angle in angles:

                phis_initial_rotation.append(angle[0])
                psis_initial_rotation.append(angle[1])


        # rotate to the initial position
        for index in range(len(phis_initial_rotation)):

            self.tools.rotate(angle_type = 'psi', angle = psis_initial_rotation[index], starting_atom = backbone[3*index+2], backbone = backbone)
            self.tools.rotate(angle_type = 'phi', angle = phis_initial_rotation[index], starting_atom = backbone[3*index+4], backbone = backbone)

        # rotate to the selected position
        for index in range(len(phi_positions)):

            self.tools.rotate(angle_type = 'psi', angle = (psis_initial_rotation[index]/rotation_steps) * 2*math.pi, starting_atom = backbone[3*index+2], backbone = backbone)
            self.tools.rotate(angle_type = 'phi', angle = (phis_initial_rotation[index]/rotation_steps) * 2*math.pi, starting_atom = backbone[3*index+4], backbone = backbone)

        return atoms

    def convert_atoms_to_configuration(self, atoms):

        configuration = {}
        for at in atoms: configuration[at.atomId] = dict(element=at.element, c_type=at.c_type, x=at.x, y=at.y, z=at.z)

        return configuration

    def write_json(self, json_data, file_name, proteinName, numberBitsRotation, method_rotations_generation):

        #Create json with calculated energies
        #TODO: extract the path to a config file
        with open(self.precalculated_energies_path+file_name+'_'+proteinName+'_'+str(numberBitsRotation)+'_'+method_rotations_generation+'.json', 'w') as outfile:
            json.dump(json_data, outfile)