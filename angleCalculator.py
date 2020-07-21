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
        self.n_angles = (number_aminoacids -1)*2

    def calculate3DStructure(self, deltas_dict, n_repetitions, beta_max, option=0):

        #Quantum calculation option for 3D structure
        if option == 0: 

            qMetropolis = quantumMetropolis.QuantumMetropolis(n_repetitions, self.bits_rotation, self.n_ancilla_bits, self.n_angles, beta_max, self.scaling_factor, deltas_dict)
            return qMetropolis.execute_quantum_metropolis_n()

        #Classical calculation option for 3D structure
        elif option == 1:

            probabilities_matrix = {}
            classical_metropolis = metropolis.Metropolis(self.bits_rotation, self.n_iterations, self.n_angles/2, self.scaling_factor, deltas_dict)
            
            for _ in range(n_repetitions):
                
                [phi, psi] = classical_metropolis.execute_metropolis()

                # it is necessary to construct the key from the received phi/psi (from the classical metropolis)
                # the idea is to add 1/n_repetitions to the returned value (to get the normalized number of times that this phi/psi was produced)
                position_angles = ''
                for index in range(len(phi)):
                    position_angles += str(phi[index]) + str(psi[index])

                if position_angles in probabilities_matrix.keys():
                    probabilities_matrix[position_angles] += (1/n_repetitions)
                else:
                    # create the new entry for this result
                    probabilities_matrix[position_angles] = (1/n_repetitions)


            print('Classical metropolis calculated for', n_repetitions, 'steps')
            return probabilities_matrix        