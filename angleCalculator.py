import quantumUtils
import math
import metropolis
import quantumMetropolis
import time
import utils

class AngleCalculator():

    def __init__(self, bits_rotation, n_ancilla_bits, number_iterations, number_aminoacids, initialization_stats, tools, qiskit_api_path, selected_device):

        self.bits_rotation = bits_rotation
        self.rotation_steps = 2**bits_rotation
        self.n_ancilla_bits = n_ancilla_bits
        self.n_iterations = number_iterations
        self.number_aminoacids = number_aminoacids

        self.qTools = quantumUtils.QuantumUtils()
        self.n_angles = (number_aminoacids -1)*2
        self.initialization_stats = initialization_stats

        self.qiskit_api_path = qiskit_api_path
        self.selected_device = selected_device

        self.tools = tools

    def calculate3DStructure(self, deltas_dict, index_min_energy, initial_step, final_step, beta, beta_type, precision_solution):

        q_accumulated_tts = []
        c_accumulated_tts = []

        min_q_tts = {'step': 0, 'value': -1}
        min_c_tts = {'step': 0, 'value': -1}

        for step in range(initial_step, final_step):

            ###### Quantum Metropolis ######
            qMetropolis = quantumMetropolis.QuantumMetropolis(step, self.bits_rotation, self.n_ancilla_bits, self.n_angles, beta, beta_type, deltas_dict, self.tools.args.mode, self.qiskit_api_path, self.selected_device)
        
            start_time = time.time()
            [probabilities_matrix, time_statevector] = qMetropolis.execute_quantum_metropolis_n()
            q_time = time.time() - start_time
            print("<i> QUANTUM METROPOLIS: Time for", step,"steps:", q_time, "seconds (", time_statevector,"seconds statevector)")
            q_tts = self.calculate_tts_from_probability_matrix(probabilities_matrix, index_min_energy, step, precision_solution)

            ###### Classical Metropolis ######
            classical_metropolis = metropolis.Metropolis(self.bits_rotation, step, self.n_angles/2, beta, beta_type, deltas_dict)
            
            start_time = time.time()
            for _ in range(self.n_iterations):

                [phi, psi] = classical_metropolis.execute_metropolis()

                # it is necessary to construct the key from the received phi/psi (from the classical metropolis)
                # the idea is to add 1/n_repetitions to the returned value (to get the normalized number of times that this phi/psi was produced)
                position_angles = ''
                for index in range(len(phi)): position_angles += str(phi[index]) + str(psi[index])

                # if the is already created, sum the entry to the dict, else create the entry
                if position_angles in probabilities_matrix.keys():
                    probabilities_matrix[position_angles] += (1/self.n_iterations) 
                else:
                    probabilities_matrix[position_angles] = (1/self.n_iterations)

            print("<i> CLASSICAL METROPOLIS: Time for", step, "steps: %s seconds" % (time.time() - start_time))
            c_tts = self.calculate_tts_from_probability_matrix(probabilities_matrix, index_min_energy, step, precision_solution)


            ###### Accumulated values Quantum Metropolis ######
            q_accumulated_tts.append(q_tts)     
            if q_tts < min_q_tts['value'] or min_q_tts['value'] == -1: min_q_tts.update(dict(value=q_tts, step=step))
        
            ###### Accumulated values Classical Metropolis ######
            c_accumulated_tts.append(c_tts)
            if c_tts < min_c_tts['value'] or min_c_tts['value'] == -1: min_c_tts.update(dict(value=c_tts, step=step))
            

            # plot data
            self.tools.plot_tts(q_accumulated_tts, c_accumulated_tts, initial_step)

            # generate json data
            final_stats = {'q': min_q_tts, 'c': min_c_tts}
            self.tools.write_tts(initial_step, final_step, q_accumulated_tts, c_accumulated_tts, self.initialization_stats, final_stats)

        return [min_q_tts, min_c_tts]

    def calculate_tts_from_probability_matrix(self, probabilities_matrix, index_min_energy, step, precision_solution):

        p_t = 0
        # if the index of min energy calculated by psi 4 is in the results of metropolis, p_t is extracted
        # else, the p_t is set to a very small value close to 0 (not 0 to avoid inf values)
        if index_min_energy in probabilities_matrix.keys():
            p_t = probabilities_matrix[index_min_energy]
        else:
            p_t = 0

        result = 0
        # Result is the calculated TTS
        if p_t >= 1:
            result = 1
        elif p_t == 0:
            result = 9999
        else:
            result = self.tools.calculateTTS(precision_solution, step, p_t)

        return result