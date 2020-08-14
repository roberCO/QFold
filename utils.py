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
            assert(angle_atoms[0].c_type == 'Carboxy' and angle_atoms[1].c_type =='C_alpha' and angle_atoms[2].element == 'N' and angle_atoms[3].c_type == 'Carboxy')
            assert(angle_atoms[0] in angle_atoms[1].linked_to and angle_atoms[1] in angle_atoms[2].linked_to and angle_atoms[2] in angle_atoms[3].linked_to)
            coords1 = np.array([angle_atoms[0].x, angle_atoms[0].y, angle_atoms[0].z])
            coords2 = np.array([angle_atoms[1].x, angle_atoms[1].y, angle_atoms[1].z])
            coords3 = np.array([angle_atoms[2].x, angle_atoms[2].y, angle_atoms[2].z])
            coords4 = np.array([angle_atoms[3].x, angle_atoms[3].y, angle_atoms[3].z])

            return self.get_dihedral(coords1, coords2, coords3, coords4)

        elif angle_type == 'psi':
            assert(angle_atoms[0].element == 'N' and angle_atoms[1].c_type =='C_alpha' and angle_atoms[2].c_type == 'Carboxy' and angle_atoms[3].element == 'N')
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

    def calculateAtomConnection(self, atoms):

        for at1 in atoms:
            for at2 in atoms:
                if at1 != at2:
                    if at1.element == 'O' and at2.element == 'C' and self.distance(at1,at2)<2:
                        at2.c_type = 'Carboxy'
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to
                        
                    if at1.element == 'N' and at2.element == 'C' and self.distance(at1,at2)<2:
                        if at2.c_type != 'Carboxy':
                            at2.c_type = 'C_alpha'
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to
                        
                    if at1.element == 'H' and self.distance(at1,at2)<1.3: 
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to
                        
                    if at1.element == 'C' and at2.element == 'C' and self.distance(at1,at2)<2  and (at1 not in at2.linked_to) : 
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to 

        return atoms

    def rotate(self, angle_type, angle, starting_atom):
        if angle_type == 'phi':
            if starting_atom.element != 'N':
                raise Exception('Wrong starting atom for the angle phi:',starting_atom.c_type )
            for atom in starting_atom.linked_to:
                if atom.c_type == 'C_alpha':
                    atom_c_alpha = atom
                    
        elif angle_type == 'psi':
            if starting_atom.c_type != 'Carboxy':
                raise Exception('Wrong starting atom for the angle phi:',starting_atom.c_type )
            for atom in starting_atom.linked_to:
                if atom.c_type == 'C_alpha':
                    atom_c_alpha = atom
        else:
            raise Exception('Angle not recognised!:',angle_type)
        
        # Define the list of atoms to rotate and then rotate them
        
        list_of_atoms_to_rotate = []
        list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type,starting_atom)
        list_of_atoms_to_rotate = self.decorations_to_rotate(list_of_atoms_to_rotate,starting_atom)

        for atom in list_of_atoms_to_rotate:
            atom.rotate(atom_c_alpha, starting_atom, angle, angle_type)  

    def backbone_to_rotate(self, angle_type, starting_atom):
        # Define the list of atoms to rotate and then rotate them
        list_of_atoms_to_rotate = []
        if angle_type == 'phi': # Follows the structure N -> Carboxy -> C_alpha -> N
            for atom in starting_atom.linked_to:
                if starting_atom.element == 'N' and atom.c_type == 'Carboxy':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type, atom)

                elif starting_atom.c_type == 'Carboxy' and atom.c_type == 'C_alpha':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type, atom)

                elif starting_atom.c_type == 'C_alpha' and atom.element == 'N':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type, atom)

                #elif atom.c_type == None:
                #list_of_atoms_to_rotate += additional_atoms()
                ## Leave this until we have backbone

        elif angle_type == 'psi': # Follows the structure N -> C_alpha -> Carboxy -> N
            for atom in starting_atom.linked_to:
                if starting_atom.c_type == 'Carboxy' and atom.element == 'N':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type = angle_type, starting_atom = atom)

                elif starting_atom.c_type == 'C_alpha' and atom.c_type == 'Carboxy':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type = angle_type, starting_atom = atom)

                elif starting_atom.element == 'N' and atom.c_type == 'C_alpha':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type = angle_type, starting_atom = atom)            
        else: 
            raise Exception('Angle_type should be either phi or psi. However it currently is {}'.format(angle_type))
            
        return list_of_atoms_to_rotate


    #Input: list of backbone atoms. Returns a list of all atoms that must be rotated.
    def decorations_to_rotate(self, backbone_list,starting_atom):

        newly_added = backbone_list
        
        while newly_added != []:
            previously_added = newly_added #Don't want the changes made to newly_added to affect to previously added
            newly_added = []
            for prev_atom in previously_added:
                for atom in prev_atom.linked_to:
                    if atom not in backbone_list and atom not in newly_added and atom != starting_atom:
                        newly_added += [atom]
            backbone_list += newly_added
            previously_added = []
        
        for atom in starting_atom.linked_to:
            if atom.element != 'N' and atom.c_type != 'Carboxy' and atom.c_type != 'C_alpha':
                backbone_list += [atom]
        
        #Perhaps use a method to eliminate duplicate references in the list? I think there shouldn't be any.
        return backbone_list

    def plotting(self, list_of_atoms, title):
    
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

    def plot_tts(self, q_accumulated_tts, c_accumulated_tts, protein_name, number_bits_rotation, method_rotations_generation, initial_step):

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
        plot_name = self.config_variables['path_tts_plot']+'tts_results_'+protein_name+'_'+str(number_bits_rotation)+'_'+method_rotations_generation+'_'+str(self.config_variables['beta_max'])+'_'+str(self.config_variables['scaling_factor'])+'.png'

        plt.savefig(plot_name, bbox_inches='tight')
        plt.close()

    def write_tts(self, initial_step, final_step, quantum_tts, classical_tts, protein_name, number_bits_rotation, method_rotations_generation):

        tts_json = {}

        tts_json['initial_step'] = initial_step
        tts_json['final_step'] = final_step
        tts_json['quantum_tts'] = quantum_tts
        tts_json['classical_tts'] = classical_tts

        json_name = self.config_variables['path_tts_plot']+'tts_results_'+protein_name+'_'+str(number_bits_rotation)+'_'+method_rotations_generation+'_'+str(self.config_variables['beta_max'])+'_'+str(self.config_variables['scaling_factor'])+'.json'
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

        plot_name = self.config_variables['path_tts_plot']+'tts_results_'+protein_name+'_'+str(number_bits_rotation)+'_'+'_'+str(self.config_variables['beta_max'])+'_'+str(self.config_variables['scaling_factor'])+'._combined.png'
        plt.savefig(plot_name, bbox_inches='tight')
        plt.close()

