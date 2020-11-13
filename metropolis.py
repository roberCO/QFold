import numpy as np
import matplotlib.pyplot as plt
import utils
import copy
import math

class Metropolis():

    ## TODO: generalise to more than 2 angles

    def __init__(self, number_angles, deltas_dict, tools):

        self.tools = tools
        self.deltas_dict = deltas_dict
        self.number_angles = int(number_angles/2)

        self.bits_rotation = self.tools.args.bits
        
        self.beta = self.tools.config_variables['beta']
        self.beta_type = self.tools.config_variables['beta_type']
        self.n_iterations = self.tools.config_variables['number_iterations']

        self.rotatition_steps = 2**self.bits_rotation
        self.bits_number_angles = math.ceil(np.log2(number_angles))


    def execute_metropolis(self, nW):

        probabilities_matrix = {}
        for _ in range(self.n_iterations):

            if self.tools.args.mode == 'real':
                [phi, psi] = self.calculate_metropolis_result(self.tools.config_variables['w_real_mode'])
            else:
                [phi, psi] = self.calculate_metropolis_result(nW)

            # it is necessary to construct the key from the received phi/psi (from the classical metropolis)
            # the idea is to add 1/n_repetitions to the returned value (to get the normalized number of times that this phi/psi was produced)
            position_angles = ''
            for index in range(len(phi)): position_angles += str(phi[index]) + str(psi[index])

            # if the is already created, sum the entry to the dict, else create the entry
            if position_angles in probabilities_matrix.keys():
                probabilities_matrix[position_angles] += (1/self.n_iterations) 
            else:
                probabilities_matrix[position_angles] = (1/self.n_iterations)

        return probabilities_matrix

    def calculate_metropolis_result(self, nW):

        #Final structure calculated with metropolis. This variable will be returned to angle calculator

        # Data structure with the rotatation (0-rotation steps) of each phi/psi angle
        # for example, if there are 3 aminoacids, there are two phis and two psi
        # the data structure for phis contains two positions the rotation for first phi and for the second phi, etc.
        anglePhi_old = []
        anglePsi_old = []

        for _ in range(self.number_angles):

            # Random initialization of angles
            anglePsi_old.append(np.random.choice(self.rotatition_steps))
            anglePhi_old.append(np.random.choice(self.rotatition_steps))

        for iteration in range(1, nW+1):

            # initially the new angles are equal to the old (then one angle will be randomly modified)
            # deep copy is necessary to avoid two pointer to the same data structure (it is necessary only to modify one of the arrays)
            anglePhi_new = copy.deepcopy(anglePhi_old)
            anglePsi_new = copy.deepcopy(anglePsi_old)
            
            # Propose a change
            # 0 = phi | 1 = psi
            change_angle = np.random.choice((0,1))

            # number of angle (it is possible to have more than one phi/psi)
            position_angle = np.random.choice(self.number_angles)
            position_angle_binary = np.binary_repr(position_angle, width = self.bits_number_angles)

            # 0 = 1 | 1 = -1
            change_plus_minus = np.random.choice((0,1))
            pm = -2*change_plus_minus + 1

            # Calculate the new angles
            if change_angle == 0:
                #Change just +1 or -1 step in the energies dictionary
                anglePhi_new[position_angle] = (anglePhi_old[position_angle] + pm) % self.rotatition_steps
            elif change_angle == 1:
                #Change just +1 or -1 step in the energies dictionary
                anglePsi_new[position_angle] = (anglePsi_old[position_angle] + pm) % self.rotatition_steps
            

            binary_key = ''
            for index in range(len(anglePhi_new)):

                # binary key should contain: phi_1 | psi_1 | phi_2 | psi_2 | ...
                binary_key += np.binary_repr(anglePhi_new[index], width = self.bits_rotation)
                binary_key += np.binary_repr(anglePsi_new[index], width = self.bits_rotation)

            # This choice of Delta_E seems weird.
            # Correspondingly: (state = angle_phi, angle_psi...) +  (move_id = phi/psi+  position_angle_binary) +  move_value

            beta_value = 0
            if self.beta_type == 'fixed':
                beta_value = self.beta
            elif self.beta_type == 'variable':
                beta_value = iteration * (self.beta / nW)
            else:
                print('<*> ERROR: Beta type wrong value. Beta type should be variable or fixed but it is', self.beta_type)

            Delta_E = self.deltas_dict[binary_key + str(change_angle) + position_angle_binary + str(change_plus_minus)]
            if Delta_E >= 0:
                    probability_threshold = math.exp(-beta_value * Delta_E)
            else: 
                probability_threshold = 1

            random_number = np.random.random_sample()

            # We should accept the change if probability_threshold > 1 (the energy goes down) or if beta is small.
            # If beta small, np.exp(-beta*Delta_E) approx 1.
            if random_number < min(1,probability_threshold): # Accept the change
                anglePhi_old = copy.deepcopy(anglePhi_new)
                anglePsi_old = copy.deepcopy(anglePsi_new)

        return [anglePhi_old, anglePsi_old]