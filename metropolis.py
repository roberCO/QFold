import numpy as np
import matplotlib.pyplot as plt
import utils

class Metropolis():

    def __init__(self, bits_rotation, n_iterations, scaling_factor, deltas_dict):

        self.bits_rotation = bits_rotation
        self.n_iterations = n_iterations
        self.scaling_factor = scaling_factor 
        self.deltas_dict = deltas_dict

        self.rotatition_steps = 2**self.bits_rotation

        self.tools = utils.Utils()

    def execute_metropolis(self):

        #Final structure calculated with metropolis. This variable will be returned to angle calculator

        # Random starting combination of angles
        anglePsi_old = np.random.choice(self.rotatition_steps)
        anglePhi_old = np.random.choice(self.rotatition_steps)

        for iteration in range(self.n_iterations):

            # Propose a change
            # 0 = phi | 1 = psi
            change_angle = np.random.choice((0,1))

            # 0 = -1 | 1 = 1
            change_plus_minus = np.random.choice((0,1))
            pm = 2*change_plus_minus - 1

            # Calculate the new angles
            if change_angle == 0:
                #Change just +1 or -1 step in the energies dictionary
                anglePhi_new = (anglePhi_old + pm) % self.rotatition_steps
                anglePsi_new = anglePsi_old
            elif change_angle == 1:
                #Change just +1 or -1 step in the energies dictionary
                anglePhi_new = anglePhi_old
                anglePsi_new = (anglePsi_old + pm) % self.rotatition_steps
            

            phi_binary = self.tools.angle_to_binary(anglePhi_new, self.bits_rotation) 
            psi_binary = self.tools.angle_to_binary(anglePsi_new, self.bits_rotation) 

            # This choice of Delta_E seems weird
            Delta_E = self.deltas_dict[phi_binary + psi_binary + str(change_angle) + str(change_plus_minus)]

            # Lets use a non_optimal simple schedule
            beta = iteration / self.n_iterations
            probability_threshold = np.exp(-beta*Delta_E)
            random_number = np.random.random_sample()

            # We should accept the change if probability_threshold > 1 (the energy goes down) or if beta is small.
            # If beta small, np.exp(-beta*Delta_E) approx 1.
            if random_number < min(1,probability_threshold): # Accept the change
                anglePhi_old = anglePhi_new
                anglePsi_old = anglePsi_new

        return [anglePhi_old, anglePsi_old]