import quantumUtils
import math
import metropolis
import quantumMetropolis
import time

class AngleCalculator():

    def __init__(self, tools, number_aminoacids):

        self.bits_rotation = tools.args.bits
        self.rotation_steps = 2**self.bits_rotation

        self.n_angles = (number_aminoacids -1)*2
        self.n_iterations = tools.config_variables['number_iterations'] * (self.rotation_steps ** self.n_angles)
        
        self.initialization = tools.args.initialization
        self.n_ancilla_bits = tools.config_variables['ancilla_bits']
        self.kappa = tools.config_variables['kappa']
        self.schedule = tools.config_variables['annealing_schedule']
        self.oracle_option = tools.config_variables['oracle_option']
        self.beta = tools.config_variables['beta']
        self.beta_type = tools.config_variables['beta_type']
        self.alpha = tools.config_variables['alpha']

    def calculate3DStructure(self, deltas_dict, n_steps, option=0):

        #Quantum calculation option for 3D structure
        if option == 'quantum': 

            qMetropolis = quantumMetropolis.QuantumMetropolis(self.initialization, n_steps, self.bits_rotation, self.n_ancilla_bits, self.n_angles, self.beta, self.beta_type, self.kappa, self.alpha, self.schedule, self.oracle_option, deltas_dict)
            
            start_time = time.time()
            [result, time_statevector] = qMetropolis.execute_quantum_metropolis_n()
            q_time = time.time() - start_time
            print("<i> QUANTUM METROPOLIS: Time for", n_steps,"steps:", q_time, "seconds (", time_statevector,"seconds statevector)")

            return result

        #Classical calculation option for 3D structure
        elif option == 'classical':

            probabilities_matrix = {}
            classical_metropolis = metropolis.Metropolis(self.initialization, self.bits_rotation, n_steps, self.n_angles/2, self.beta, self.beta_type, self.alpha, self.schedule, self.kappa, deltas_dict)
            
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