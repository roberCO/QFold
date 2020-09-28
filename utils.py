import atom
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import struct
import copy
import math
import json

class Utils():

    def __init__(self, config_path=''):

        if(config_path != ''):
            try:
                f = open(config_path)
                f.close()
            except IOError:
                print('<!> Info: No configuration file')
                raise Exception('It is necessary to create a configuration file (.json) for some variables')

            with open(config_path) as json_file:
                        self.config_variables = json.load(json_file)

    def get_config_variables(self):
        return self.config_variables

    def get_dihedral(self, coords1, coords2, coords3, coords4):
        """Returns the dihedral angle in degrees."""

        a1 = coords2 - coords1
        a2 = coords3 - coords2
        a3 = coords4 - coords3

        v1 = np.cross(a1, a2)
        v1 = v1 / (v1 * v1).sum(-1)**0.5
        v2 = np.cross(a2, a3)
        v2 = v2 / (v2 * v2).sum(-1)**0.5
        porm = np.sign((v1 * a3).sum(-1))
        #Round the value to avoid 1.0000000000002 (python precision error)
        rad = np.arccos(round((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5, 10))
        if not porm == 0:
            rad = rad * porm
        return rad

    def calculateAngle(self, angle_atoms, angle_type):
        'Uses get dihedral to calculate angles between atoms'
        if angle_type == 'phi':
            # For angle phi we take the last atom of the previous aminoacid
            assert(angle_atoms[0].c_type == 'Carboxy' and angle_atoms[1].c_type == 'N_backbone' and angle_atoms[2].c_type =='C_alpha' and angle_atoms[3].c_type == 'Carboxy')
            assert(angle_atoms[0] in angle_atoms[1].linked_to and angle_atoms[1] in angle_atoms[2].linked_to and angle_atoms[2] in angle_atoms[3].linked_to)
            coords1 = np.array([angle_atoms[0].x, angle_atoms[0].y, angle_atoms[0].z])
            coords2 = np.array([angle_atoms[1].x, angle_atoms[1].y, angle_atoms[1].z])
            coords3 = np.array([angle_atoms[2].x, angle_atoms[2].y, angle_atoms[2].z])
            coords4 = np.array([angle_atoms[3].x, angle_atoms[3].y, angle_atoms[3].z])

            return self.get_dihedral(coords1, coords2, coords3, coords4)

        elif angle_type == 'psi':
            # For angle psi we take the first atom of the next aminoacid
            assert(angle_atoms[0].c_type == 'N_backbone' and angle_atoms[1].c_type =='C_alpha' and angle_atoms[2].c_type == 'Carboxy' and angle_atoms[3].c_type == 'N_backbone')
            assert(angle_atoms[0] in angle_atoms[1].linked_to and angle_atoms[1] in angle_atoms[2].linked_to and angle_atoms[2] in angle_atoms[3].linked_to)
            coords1 = np.array([angle_atoms[0].x, angle_atoms[0].y, angle_atoms[0].z])
            coords2 = np.array([angle_atoms[1].x, angle_atoms[1].y, angle_atoms[1].z])
            coords3 = np.array([angle_atoms[2].x, angle_atoms[2].y, angle_atoms[2].z])
            coords4 = np.array([angle_atoms[3].x, angle_atoms[3].y, angle_atoms[3].z])

            return self.get_dihedral(coords1, coords2, coords3, coords4)

        else:
            raise('Angle not recognised!:'+str(angle_type))


    def distance(self, atom, atom2):
        return np.sqrt((atom.x-atom2.x)**2+(atom.y-atom2.y)**2+(atom.z-atom2.z)**2)

    def calculateAtomConnection(self, atoms, aminoacids):

        #Let us first map the topology. Currently cost is O(N^2). Some other algorithm could be desirable
        for at1 in atoms:
            for at2 in atoms:
                if at1 != at2:
                    if at1.element != 'H' and at2.element != 'H' and self.distance(at1,at2)<2  and (at1 not in at2.linked_to): 
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to

                    elif at1.element != 'H' and at2.element == 'H' and self.distance(at1,at2)<1.3  and (at1 not in at2.linked_to):
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to

        # Next we give an structure to each linked_to list
        for at in atoms:
            at.linked_to_dict = {'N': [], 'O': [], 'C': [], 'H': [], 'Other': []}
            for at1 in at.linked_to:
                if at1.element == 'N':
                    at.linked_to_dict['N'].append(at1)
                elif at1.element == 'O':
                    at.linked_to_dict['O'].append(at1)
                elif at1.element == 'C':
                    at.linked_to_dict['C'].append(at1)
                elif at1.element == 'H':
                    at.linked_to_dict['H'].append(at1)
                else:
                    at.linked_to_dict['Other'].append(at1)

        #self.plotting(list_of_atoms = atoms, title = 'Peptide_plot')

        #make a list of nitrogen atoms where one could start the main chain
        nitrogen_starts = []

        # For any aminoacid except proline
        if aminoacids[0] != 'P':
            for at in atoms:
                # This allows to identify any initial N to start except for Proline which has a weird structure
                if at.element == 'N' and len(at.linked_to_dict['C']) == 1 and len(at.linked_to_dict['H'])==2:
                    nitrogen_starts.append(at)

        # For the protein starting at proline
        elif aminoacids[0] == 'P':
            for at in atoms:
                # This allows to identify any initial N to start except for Proline which has a weird structure
                if at.element == 'N' and self.is_proline_N(at):
                    nitrogen_starts.append(at)

        # Find main_chain
        backbone = self.main_chain_builder(nitrogen_starts, aminoacids)

        # Name the atoms
        for (atom,i) in zip(backbone, range(len(backbone))):
            if atom.element == 'N' and (i % 3 == 0):
                atom.c_type = 'N_backbone'
            elif atom.element == 'C' and (i % 3 == 1) and (atom.linked_to_dict['O'] == []):
                atom.c_type = 'C_alpha'
            elif atom.element == 'C' and (i % 3 == 2) and (atom.linked_to_dict['O'] != []):
                atom.c_type = 'Carboxy'
            else:
                raise TypeError('The atom', atom.element, 'does not fulfill the requirements to be part of the backbone')

        return atoms, backbone

    def is_proline_N(self, atom):

        carbon_ring = []

        if atom.element != 'N' or len(atom.linked_to_dict['C']) != 2 or len(atom.linked_to_dict['H'])!=1:
            return False
        else:
            carbons = atom.linked_to_dict['C']
            if len(carbons[0].linked_to_dict['C']) == 1  and len(carbons[0].linked_to_dict['N']) == 1 and len(carbons[1].linked_to_dict['C']) == 2 and len(carbons[1].linked_to_dict['N']) == 1:
                current_carbon = carbons[0]
                ending_carbon = carbons[1]
            elif len(carbons[1].linked_to_dict['C']) == 1  and len(carbons[1].linked_to_dict['N']) == 1 and len(carbons[0].linked_to_dict['C']) == 2 and len(carbons[0].linked_to_dict['N']) == 1:
                current_carbon = carbons[1]
                ending_carbon = carbons[0]
            else:
                return False
            
            for _ in range(2):
                carbon_ring.append(current_carbon)
                current_carbon = (current_carbon.linked_to_dict['C'][0] if current_carbon.linked_to_dict['C'][0] not in carbon_ring else current_carbon.linked_to_dict['C'][1])
                if len(current_carbon.linked_to_dict['C']) != 2 or len(current_carbon.linked_to_dict['N']) != 0 or len(current_carbon.linked_to_dict['O']) != 0 or len(current_carbon.linked_to_dict['H']) != 2:
                    return False
                
            return (True if current_carbon in ending_carbon.linked_to else False)
            


    def main_chain_builder(self, nitrogen_starts, aminoacids):
        '''Takes all the nitrogens that are only connected to a single C and returns the backbone of the protein'''
        best_chains = []
        len_best_chain = 0
        for nitro in nitrogen_starts:
            candidate_chain = []
            nit = nitro
            amino_index = 0
            while True:

                aminolist = []
                aminolist.append(nit)

                # Searching for C-alpha
                carbons = nit.linked_to_dict['C']
                carbons_not_in_chain = [carbon for carbon in carbons if (carbon not in candidate_chain and carbon not in aminolist)]
                if (len(carbons_not_in_chain)==1 and aminoacids[amino_index] != 'P'):
                    car_alpha = carbons_not_in_chain[0]
                    aminolist.append(car_alpha)
                elif (len(carbons_not_in_chain)==2 and aminoacids[amino_index] == 'P'):
                    car_alpha = (carbons_not_in_chain[0] if (len(carbons_not_in_chain[0].linked_to_dict['N']) == 1 and len(carbons_not_in_chain[0].linked_to_dict['C']) == 2 and len(carbons_not_in_chain[0].linked_to_dict['H']) == 1) else carbons_not_in_chain[1])
                    aminolist.append(car_alpha)
                else:
                    break

                # Searching for Carboxy
                carbons = car_alpha.linked_to_dict['C']
                carboxys_not_in_chain = [carbon for carbon in carbons if (carbon not in candidate_chain and carbon not in aminolist and len(carbon.linked_to_dict['O']) > 0)]
                if len(carboxys_not_in_chain)==1:
                    carbox = carboxys_not_in_chain[0]
                    aminolist.append(carbox)
                else:
                    break

                #We have a full aminoacid, so we save it to the candidate list
                candidate_chain += aminolist
                amino_index += 1

                # Searching for next aminoacid Nitrogen
                nitrogens = carbox.linked_to_dict['N']
                nitrogens_not_in_chain = [n for n in nitrogens if (n not in candidate_chain and n not in aminolist)]
                if len(nitrogens_not_in_chain)==1:
                    nit = nitrogens_not_in_chain[0]
                else:
                    break

            # Is the found chain longer than the one we already had?
            if len(candidate_chain) > len_best_chain:
                len_best_chain = len(candidate_chain)
                best_chains = [candidate_chain]
            elif len(candidate_chain) == len_best_chain:
                best_chains.append(candidate_chain)
            else: 
                pass

        if len(best_chains) != 1 or len(best_chains[0])//3 != len(aminoacids):
            raise ValueError('There should be a single lengthy chain!', best_chains, nitrogen_starts)
        else:
            return best_chains[0]


    def rotate(self, angle_type, angle, starting_atom, backbone):

        previous_atom = backbone[backbone.index(starting_atom)-1]

        if angle_type == 'phi':
            if previous_atom.c_type != 'N_backbone' or starting_atom.c_type != 'C_alpha':
                raise Exception('Wrong starting atom for the angle phi:',starting_atom.c_type,'or wrong previous atom',previous_atom.c_type )
                    
        elif angle_type == 'psi':
            if previous_atom.c_type != 'C_alpha' or starting_atom.c_type != 'Carboxy':
                raise Exception('Wrong starting atom for the angle phi:',starting_atom.c_type )

        else:
            raise Exception('Angle not recognised!:',angle_type)
        
        # Define the list of atoms to rotate and then rotate them
        
        backbone2rotate = backbone[backbone.index(starting_atom)+1:]
        ##self.backbone_to_rotate(angle_type,starting_atom, backbone)
        list_of_atoms_to_rotate = self.decorations_to_rotate(backbone2rotate,backbone)

        for atom in list_of_atoms_to_rotate:
            # The axis is defined by the starting atom and the atom prior to the starting atom in the backbone
            atom.rotate(previous_atom, starting_atom, angle, angle_type)  

    #Input: list of backbone atoms. Returns a list of all atoms that must be rotated.
    def decorations_to_rotate(self, backbone2rotate, backbone):
        
        atoms2rotate = backbone2rotate

        newly_added = backbone2rotate

        while newly_added != []:
            previously_added = newly_added
            newly_added = []
            for atom in previously_added:
                for at2 in atom.linked_to:
                    if at2 not in atoms2rotate and at2 not in newly_added and at2 not in backbone:
                        newly_added.append(at2)

            atoms2rotate += newly_added 
        
        return atoms2rotate


    def plotting(self, list_of_atoms, title, plane = False):
    
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

        planePhiPoints = []
        planePsiPoints = []

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

            if at.atomId == 3 or at.atomId == 5 or at.atomId == 6:
                planePhiPoints += [np.array([at.x, at.y, at.z])]
            
            if at.atomId == 4 or at.atomId == 6 or at.atomId == 7:
                planePsiPoints += [np.array([at.x, at.y, at.z])]

        ax.scatter(xs, ys, zs,c=c,depthshade= False)

        if plane:
            XPhi, YPhi, ZPhi = self.calculatePlane(planePhiPoints)
            XPsi, YPsi, ZPsi = self.calculatePlane(planePsiPoints)

            # plot the mesh. Each array is 2D, so we flatten them to 1D arrays
            ax.plot_surface(XPhi, YPhi, ZPhi, color = 'r', alpha = '0.5')
            ax.plot_surface(XPsi, YPsi, ZPsi, color = 'g', alpha = '0.5')


        for i in range(len(xs)): 
            ax.text(xs[i],ys[i],zs[i],  '%s' % (str(i)))
        plt.show()

    def calculatePlane(self, planePoints):

        # These two vectors are in the plane
        v1 = planePoints[2] - planePoints[0]
        v2 = planePoints[1] - planePoints[0]

        # the cross product is a vector normal to the plane
        cp = np.cross(v1, v2)
        a, b, c = cp

        # This evaluates a * x3 + b * y3 + c * z3 which equals d
        d = np.dot(cp, planePoints[2])


        #HARDCODED values
        X, Y = np.meshgrid(planePoints[0], planePoints[0])

        Z = (d - a * X - b * Y) / c

        return X, Y, Z

    def decode_angle_from_index(self, numberBitsRotation, value_angle_psi4, index):

        step = (2*math.pi)/(2**numberBitsRotation)
        calculated_angle_value = value_angle_psi4 + index*step

        if calculated_angle_value > math.pi:
            calculated_angle_value -= 2*math.pi

        return calculated_angle_value

    def calculatePrecisionOfAngles(self, phi_angles_psi4, psi_angles_psi4, phis_initial_rotation, psis_initial_rotation):

        if len(phi_angles_psi4) != len(phis_initial_rotation) or len(psi_angles_psi4) != len(phis_initial_rotation):
            print('<*> ERROR: The number of generated angles (initialization) is different than the number of protein angles')

        phi_precisions = []
        psi_precisions = []

        for index in range(len(phi_angles_psi4)):

            option_1 = 0
            option_2 = 0

            if phis_initial_rotation[index] > phi_angles_psi4[index]:

                #Calculate the distance if the angles go to zero and starts again
                option_1 = abs(math.pi - phis_initial_rotation[index]) + abs(-math.pi - phi_angles_psi4[index])

            else:

                #Calculate the distance if the angles go to zero and starts again
                option_1 = abs(math.pi - phi_angles_psi4[index]) + abs(-math.pi - phis_initial_rotation[index])

            #option_2 is common for both previous cases
            option_2 = abs(phis_initial_rotation[index] -  phi_angles_psi4[index])
            
            minimum_option = min(option_1, option_2)
            phi_precisions.append((1-(minimum_option / (2*math.pi)))*100)

        for index in range(len(psi_angles_psi4)):

            if psis_initial_rotation[index] > psi_angles_psi4[index]:

                #Calculate the distance if the angles go to zero and starts again
                option_1 = abs(math.pi - psis_initial_rotation[index]) + abs(-math.pi - psi_angles_psi4[index])

            else:

                #Calculate the distance if the angles go to zero and starts again
                option_1 = abs(math.pi - psi_angles_psi4[index]) + abs(-math.pi - psis_initial_rotation[index])
                
            option_2 = abs(psis_initial_rotation[index] - psi_angles_psi4[index])
            minimum_option = min(option_1, option_2)
            psi_precisions.append((1-(minimum_option / (2*math.pi)))*100)
        
        print('\nPHI precision: ', np.mean(phi_precisions), '% phi mean real value: ', np.mean(phi_angles_psi4), 'phi mean calculated value:', np.mean(phis_initial_rotation))
        print('PSI precision: ', np.mean(psi_precisions), '% psi mean real value: ', np.mean(psi_angles_psi4), 'psi mean calculated value:', np.mean(psis_initial_rotation), '\n')

        return [phi_precisions, psi_precisions]

    def angle_to_binary(self, angle, number_bits_rotation):

        return ('0'*(number_bits_rotation - len(format(angle,'b'))) + format(angle,'b'))

    def calculateTTS(self, precision_solution, t, p_t):

        return t * (math.log10(1-precision_solution)/(math.log10(1-p_t)))

    def calculate_diff_vs_mean_diffs(self, min_energy_difference, deltas_mean):

        return (1 - (min_energy_difference/deltas_mean)) * 100

    def calculate_delta_mean(self, deltas_dict):

        array = np.array(list(deltas_dict.items()), dtype='float32')
        return array[:,1].mean()

    def calculate_std_dev_deltas(self, deltas_dict):

        array = np.array(list(deltas_dict.items()), dtype='float32')
        return array[:,1].std()

    def plot_tts(self, q_accumulated_tts, c_accumulated_tts, protein_name, aminoacids, number_bits_rotation, method_rotations_generation, initial_step):

        x_axis = [x for x in range(initial_step, initial_step+len(q_accumulated_tts))]

        fig = plt.figure()
        
        ax = fig.add_subplot(111)
        ax.set_title('TTS comparision for Quantum vs Classical Metropolis')

        interval = math.ceil((initial_step+len(q_accumulated_tts)-initial_step)/10)
        plt.xticks(np.arange(min(x_axis), max(x_axis)+1, interval))

        ax.plot(x_axis, q_accumulated_tts, marker='o', markersize=3, color="red", label = 'q_tts')
        ax.plot(x_axis, c_accumulated_tts, marker='o', markersize=3, color="blue", label = 'c_tts')
            
        ax.set_ylabel('TTS')
        ax.set_xlabel('Steps')
        plt.tight_layout()

        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plot_name = self.config_variables['path_tts_plot']+'tts_results_'+protein_name+'_'+aminoacids+'_'+str(number_bits_rotation)+'_'+method_rotations_generation+'_'+str(self.config_variables['beta_max'])+'_'+str(self.config_variables['scaling_factor'])+'.png'

        plt.savefig(plot_name, bbox_inches='tight')
        plt.close()

    def write_tts(self, initial_step, final_step, quantum_tts, classical_tts, protein_name, aminoacids, number_bits_rotation, method_rotations_generation, inizialitation_stats, final_stats):

        tts_json = {}

        tts_json['initial_step'] = initial_step
        tts_json['final_step'] = final_step
        tts_json['quantum_tts'] = quantum_tts
        tts_json['classical_tts'] = classical_tts
        tts_json['initialization_stats'] = inizialitation_stats
        tts_json['final_stats'] = final_stats

        json_name = self.config_variables['path_tts_plot']+'tts_results_'+protein_name+'_'+aminoacids+'_'+str(number_bits_rotation)+'_'+method_rotations_generation+'_'+str(self.config_variables['beta_max'])+'_'+str(self.config_variables['scaling_factor'])+'.json'
        with open(json_name, 'w') as outfile:
            json.dump(tts_json, outfile)

    def read_results_file(self, path_file):

        with open(path_file) as json_file:
            data = json.load(json_file)
        
        return data

    def generate_combined_results_plot(self, data, protein_name, number_bits_rotation):

        fig = plt.figure()

        initial_steps = []
        final_steps = []
        methods = ''
        for key in data:

            initial_steps.append(data[key]['initial_step'])
            final_steps.append(data[key]['final_step'])
            methods += key + ' | '
                
        ax = fig.add_subplot(111)
        # in methods [:-3] remove last ' | '
        ax.set_title('Combined TTS comparision Q vs C for '+ methods[:-3])

        interval = math.ceil((max(final_steps) - min(initial_steps))/10)
        plt.xticks(np.arange(min(initial_steps), max(final_steps), interval))

        for key in data:

            ax.plot([x for x in range(data[key]['initial_step'], data[key]['final_step'])], data[key]['quantum_tts'], marker='o', markersize=3, color=self.config_variables['color_quantum_'+key], label = 'q_tts_'+key)
            ax.plot([x for x in range(data[key]['initial_step'], data[key]['final_step'])], data[key]['classical_tts'], marker='o', markersize=3, color=self.config_variables['color_classical_'+key], label = 'c_tts_'+key)
            
        ax.set_ylabel('TTS')
        ax.set_xlabel('Steps')
        plt.tight_layout()

        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plt.tight_layout()

        plot_name = self.config_variables['path_tts_plot']+'tts_results_'+protein_name+'_'+str(number_bits_rotation)+'_'+'_'+str(self.config_variables['beta_max'])+'_'+str(self.config_variables['scaling_factor'])+'_combined.png'
        plt.savefig(plot_name, bbox_inches='tight')
        plt.close()

    def read_results_data(self, input_name):

        results = {}
        data = {}
        # read data
        path = './'+self.config_variables['path_tts_plot']+ 'tts_results_'+input_name+'.json'
        with open(path) as json_file:
                    data[input_name] = json.load(json_file)

        # prepare the data
        for protein_key in data.keys():

            stats = {}

            aas = protein_key.split('_')[1]
            bits = protein_key.split('_')[2]
            init_mode = protein_key.split('_')[3]
            phi_prec = data[protein_key]['initialization_stats']['phis_precision']
            psi_prec = data[protein_key]['initialization_stats']['psis_precision']
            
            # divided by 2 because it is the mean between phi and psi and by 100 to normalize the precision
            stats['precision'] = np.mean(np.mean(phi_prec) + np.mean(psi_prec))/2/100
            stats['min_tts'] = min(data[protein_key]['final_stats']['q']['value'], data[protein_key]['final_stats']['c']['value'])
            stats['quantum_tts'] = data[protein_key]['quantum_tts']
            stats['classical_tts'] = data[protein_key]['classical_tts']
            stats['initial_step'] = data[protein_key]['initial_step']
            stats['final_step'] = data[protein_key]['final_step']
            stats['min_tts_q'] = data[protein_key]['final_stats']['q']['value']
            stats['min_tts_c'] = data[protein_key]['final_stats']['c']['value']
            
            results[aas+'_'+bits+'_'+init_mode] = stats

        return results
