import quantumUtils
import math
import metropolis
import quantumMetropolis
import time

class AngleCalculator():

    def __init__(self, initialization, bits_rotation, n_ancilla_bits, number_iterations, oracle_option, number_aminoacids):

        self.initialization = initialization
        self.bits_rotation = bits_rotation
        self.rotation_steps = 2**bits_rotation
        self.n_ancilla_bits = n_ancilla_bits
        self.oracle_option = oracle_option
        self.number_aminoacids = number_aminoacids

        self.qTools = quantumUtils.QuantumUtils()
        self.n_angles = (number_aminoacids -1)*2
        self.n_iterations = number_iterations * self.n_angles * self.rotation_steps

    def calculate3DStructure(self, deltas_dict, n_steps, beta, beta_type, option=0):

        #Quantum calculation option for 3D structure
        if option == 0: 

            qMetropolis = quantumMetropolis.QuantumMetropolis(self.initialization, n_steps, self.bits_rotation, self.n_ancilla_bits, self.n_angles, beta, beta_type, self.oracle_option, deltas_dict)
            
            start_time = time.time()
            [result, time_statevector] = qMetropolis.execute_quantum_metropolis_n()
            q_time = time.time() - start_time
            print("<i> QUANTUM METROPOLIS: Time for", n_steps,"steps:", q_time, "seconds (", time_statevector,"seconds statevector)")

            return result

        #Classical calculation option for 3D structure
        elif option == 1:

            probabilities_matrix = {}
            classical_metropolis = metropolis.Metropolis(self.initialization, self.bits_rotation, n_steps, self.n_angles/2, beta, beta_type, deltas_dict)
            
            start_time = time.time()
            for _ in range(self.n_iterations):

                [phi, psi] = classical_metropolis.execute_metropolis()

                # it is necessary to construct the key from the received phi/psi (from the classical metropolis)
                # the idea is to add 1/n_repetitions to the returned value (to get the normalized number of times that this phi/psi was produced)
                position_angles = ''
                for index in range(len(phi)):

                    if index != 0: position_angles += '-'
                    position_angles += str(phi[index]) + '-' + str(psi[index])

                if position_angles in probabilities_matrix.keys():
                    probabilities_matrix[position_angles] += (1/self.n_iterations)
                else:
                    # create the new entry for this result
                    probabilities_matrix[position_angles] = (1/self.n_iterations)

            print("<i> CLASSICAL METROPOLIS: Time for", n_steps, "steps: %s seconds" % (time.time() - start_time))
            
            return probabilities_matrix        