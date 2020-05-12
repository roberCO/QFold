import numpy as np

class Metropolis():

    def __init__(self, n_iterations, scaling_factor, energies):

        self.n_iterations = n_iterations
        self.scaling_factor = scaling_factor 
        self.energies = energies

    def execute_metropolis(self):

        calculated_3d_structure = None

        # Random starting combination of angles
        anglePsi_old = np.random.choice(len(self.energies))
        anglePhi_old = np.random.choice(len(self.energies))

        for iteration in range(self.n_iterations):

            # First retrieve the present energy
            E_old = self.energies[(anglePhi_old,anglePsi_old)]
            
            # Propose a change
            change_angle = np.random.choice(('phi','psi'))
            change_plus_minus = np.random.choice((1,-1))

            # Calculate the new angles
            if change_angle == 'phi':
                #Change just +1 or -1 step in the energies dictionary
                anglePhi_new = anglePhi_old + change_plus_minus * 1
                anglePsi_new = anglePsi_old
            elif change_angle == 'psi':
                #Change just +1 or -1 step in the energies dictionary
                anglePhi_new = anglePhi_old
                anglePsi_new = anglePsi_new + change_plus_minus * 1
            
            # Calculate the new energy
            if self.energies[(anglePhi_new,anglePsi_new)] != None:
                E_new = self.energies[(anglePhi_new,anglePsi_new)]

            else:

                #Read the PSI4 output file and get the energy
                E_new = self.energies()
                self.energies[(anglePhi_new,anglePsi_new)] = E_new

            Delta_E = (E_new - E_old) * self.scaling_factor

            # Lets use a non_optimal simple schedule
            beta = iteration / self.n_iterations
            probability_threshold = np.exp(-beta*Delta_E)
            random_number = np.random.random_sample

            # We should accept the change if probability_threshold > 1 (the energy goes down) or if beta is small.
            # If beta small, np.exp(-beta*Delta_E) approx 1.
            if random_number < min(1,probability_threshold): # Accept the change
                anglePhi_old = anglePhi_new
                anglePsi_old = anglePsi_new

        return calculated_3d_structure