import quantumUtils
import math
import metropolis
import quantumMetropolis
import progressbar

class AngleCalculator():

    def __init__(self, bits_rotation, n_ancilla_bits, scaling_factor, number_iterations):

        self.bits_rotation = bits_rotation
        self.rotationSteps = 2**bits_rotation
        self.n_ancilla_bits = n_ancilla_bits
        self.scaling_factor = scaling_factor
        self.n_iterations = number_iterations
        self.qTools = quantumUtils.QuantumUtils()

    def calculate3DStructure(self, energyList, n_repetitions, beta_max, option=0):

        #Quantum calculation option for 3D structure
        if option == 0: 

            print('\n## Quantum Metropolis ##')
            print('    ⬤  Calculating p_t for', n_repetitions,'steps')

            qMetropolis = quantumMetropolis.QuantumMetropolis(n_repetitions, self.bits_rotation, self.n_ancilla_bits, beta_max, energyList)
            return qMetropolis.execute_quantum_metropolis()

        #Classical calculation option for 3D structure
        elif option == 1:

            print('\n## Classical Metropolis ##')
            print('    ⬤  Calculating p_t for', n_repetitions,'steps')

            probabilities_matrix = [[0]*self.rotationSteps for x in range(self.rotationSteps)]
            classical_metropolis = metropolis.Metropolis(self.n_iterations, self.scaling_factor, energyList)
            
            print('    ⬤  Iterating the classical metropolis')
            bar = progressbar.ProgressBar(maxval=n_repetitions, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
            bar.start()
            for index in range(n_repetitions):
                
                [phi, psi] = classical_metropolis.execute_metropolis()
                probabilities_matrix[phi][psi] += (1/n_repetitions)

                bar.update(index+1)
                
            bar.finish()

            return probabilities_matrix
        