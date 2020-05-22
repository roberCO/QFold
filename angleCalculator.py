import quantumUtils
import math
import metropolis
import quantumMetropolis

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

            print('<i> Quantum Metropolis calculating p_t for', n_repetitions,'steps')

            qMetropolis = quantumMetropolis.QuantumMetropolis(n_repetitions, self.bits_rotation, self.n_ancilla_bits, beta_max, energyList)
            return qMetropolis.execute_quantum_metropolis()

        #Classical calculation option for 3D structure
        elif option == 1:

            print('<i> Classical Metropolis calculating p_t for', n_repetitions,'steps')

            probabilities_matrix = [[0]*self.rotationSteps for x in range(self.rotationSteps)]
            
            for _ in range(n_repetitions):

                classical_metropolis = metropolis.Metropolis(self.n_iterations, self.scaling_factor, energyList)
                
                [phi, psi] = classical_metropolis.execute_metropolis()
                probabilities_matrix[phi][psi] += (1/n_repetitions)
                
            return probabilities_matrix
        