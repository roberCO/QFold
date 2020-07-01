import quantumUtils
import math
import metropolis
import quantumMetropolis

class AngleCalculator():

    def __init__(self, bits_rotation, n_ancilla_bits, scaling_factor, number_iterations, number_aminoacids):

        self.bits_rotation = bits_rotation
        self.rotation_steps = 2**bits_rotation
        self.n_ancilla_bits = n_ancilla_bits
        self.scaling_factor = scaling_factor
        self.n_iterations = number_iterations
        self.number_aminoacids = number_aminoacids
        self.qTools = quantumUtils.QuantumUtils()

    def calculate3DStructure(self, deltas_dict, n_repetitions, beta_max, option=0):

        #Quantum calculation option for 3D structure
        if option == 0: 

            qMetropolis = quantumMetropolis.QuantumMetropolis(n_repetitions, self.bits_rotation, self.n_ancilla_bits, beta_max, deltas_dict)
            return qMetropolis.execute_quantum_metropolis()

        #Classical calculation option for 3D structure
        elif option == 1:

            probabilities_matrix = [[0]*self.rotation_steps for x in range(self.rotation_steps)]
            classical_metropolis = metropolis.Metropolis(self.bits_rotation, self.n_iterations, self.number_aminoacids, self.scaling_factor, deltas_dict)
            
            for _ in range(n_repetitions):
                
                [phi, psi] = classical_metropolis.execute_metropolis()
                probabilities_matrix[phi][psi] += (1/n_repetitions)

            return probabilities_matrix        