import numpy as np
import matplotlib.pyplot as plt
import utils
import copy
import math

class Metropolis():

    def __init__(self, number_angles, deltas_dict, tools):

        self.tools = tools

        self.initialization = tools.args.initialization
        self.bits_rotation = tools.args.bits
        self.mode = tools.args.mode
        self.beta = tools.config_variables['beta']
        self.beta_type = tools.config_variables['beta_type']
        self.kappa = tools.config_variables['kappa']
        self.alpha = tools.config_variables['alpha']
        self.annealing_schedule = tools.config_variables['annealing_schedule']
        self.w_real_mode = tools.config_variables['w_real_mode']

        self.deltas_dict = deltas_dict
        self.number_angles = int(number_angles/2)
        self.rotation_steps = 2**self.bits_rotation
        self.bits_number_angles = math.ceil(np.log2(self.number_angles))

        self.n_iterations = tools.config_variables['number_iterations'] * (self.rotation_steps ** self.number_angles)

    def execute_metropolis(self, nW):

        probabilities_matrix = {}
        for _ in range(self.n_iterations):
            
            [phi, psi] = self.calculate_metropolis_result(nW)
        
            # it is necessary to construct the key from the received phi/psi (from the classical metropolis)
            # the idea is to add 1/n_repetitions to the returned value (to get the normalized number of times that this phi/psi was produced)
            position_angles = ''
            for index in range(len(phi)): position_angles += str(phi[index]) + '-' + str(psi[index]) + "-"
            position_angles = position_angles[:-1]

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

        if self.initialization == 'random':

            for _ in range(self.number_angles):

                # Random initialization of angles
                anglePsi_old.append(np.random.choice(self.rotation_steps))
                anglePhi_old.append(np.random.choice(self.rotation_steps))

        elif self.initialization == 'minifold':

            _, accumulated_probs = self.tools.von_mises_amplitudes(n_qubits = self.bits_rotation, kappa = self.kappa)

            for _ in range(self.number_angles):

                # Random initialization of angles
                r = np.random.random()
                for i in range(self.rotation_steps):
                    if accumulated_probs[i] >= r:
                        anglePhi_old.append(i)
                        break

                r = np.random.random()
                for i in range(self.rotation_steps):
                    if accumulated_probs[i] >= r:
                        anglePsi_old.append(i)
                        break

        for i in range(1, nW+1):

            anglePhi_new, anglePsi_new, change_angle, position_angle, change_plus_minus = self.generate_new_angles(anglePhi_old, anglePsi_old)
            position_angle_binary = np.binary_repr(position_angle, width = self.bits_number_angles)            

            binary_key = ''
            for index in range(len(anglePhi_new)):

                # binary key should contain: phi_1 | psi_1 | phi_2 | psi_2 | ...
                binary_key += np.binary_repr(anglePhi_old[index], width = self.bits_rotation)
                binary_key += np.binary_repr(anglePsi_old[index], width = self.bits_rotation)

            # This choice of Delta_E seems weird.
            # Correspondingly: (state = angle_phi, angle_psi...) +  (move_id = phi/psi+  position_angle_binary) +  move_value
            beta_value = 0
            if self.beta_type == 'fixed':
                beta_value = self.beta
            elif self.beta_type == 'variable':
                if self.annealing_schedule == 'Cauchy' or self.annealing_schedule == 'linear':
                    beta_value = self.beta * i 
                elif self.annealing_schedule == 'Boltzmann' or self.annealing_schedule == 'logarithmic':
                    beta_value = self.beta * np.log(i) + self.beta
                elif self.annealing_schedule == 'geometric':
                    beta_value = self.beta * self.alpha**(-i+1)
                elif self.annealing_schedule == 'exponential': 
                    space_dim = self.number_angles
                    beta_value = self.beta * np.exp( self.alpha * (i-1)**(1/space_dim) )
                else:
                    raise ValueError('<*> ERROR: Annealing Scheduling wrong value. It should be one of [linear, logarithmic, geometric, exponential] but it is', self.annealing_schedule)
            else:
                ValueError('<*> ERROR: Beta type wrong value. Beta type should be variable or fixed but it is', self.beta_type)

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
    
    def generate_new_angles(self, anglePhi_old, anglePsi_old):
        # initially the new angles are equal to the old (then one angle will be randomly modified)
        # deep copy is necessary to avoid two pointer to the same data structure (it is necessary only to modify one of the arrays)
        anglePhi_new = copy.deepcopy(anglePhi_old)
        anglePsi_new = copy.deepcopy(anglePsi_old)
        
        # Propose a change
        # 0 = phi | 1 = psi
        change_angle = np.random.choice((0,1))

        # number of angle (it is possible to have more than one phi/psi)
        position_angle = np.random.choice(self.number_angles)

        # 0 = 1 | 1 = -1
        change_plus_minus = np.random.choice((0,1))
        pm = -2*change_plus_minus + 1

        # Calculate the new angles
        if change_angle == 0:
            #Change just +1 or -1 step in the energies dictionary
            anglePhi_new[position_angle] = (anglePhi_old[position_angle] + pm) % self.rotation_steps
        elif change_angle == 1:
            #Change just +1 or -1 step in the energies dictionary
            anglePsi_new[position_angle] = (anglePsi_old[position_angle] + pm) % self.rotation_steps

        return anglePhi_new, anglePsi_new, change_angle, position_angle, change_plus_minus