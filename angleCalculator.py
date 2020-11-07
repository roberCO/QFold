import math
import metropolis
import quantumMetropolis
import time
import utils

class AngleCalculator():

    def __init__(self, tools, initialization_stats):

        self.tools = tools
        self.initialization_stats = initialization_stats

        self.n_iterations = self.tools.config_variables['number_iterations']

        self.n_angles = (len(self.tools.args.aminoacids) -1)*2

    def calculate3DStructure(self, deltas_dict, index_min_energy):

        q_accumulated_tts = []
        c_accumulated_tts = []

        min_q_tts = {'step': 0, 'value': -1}
        min_c_tts = {'step': 0, 'value': -1}

        for step in range(self.tools.config_variables['initial_step'], self.tools.config_variables['final_step']):

            ###### Quantum Metropolis ######
            qMetropolis = quantumMetropolis.QuantumMetropolis(step, self.n_angles, deltas_dict, self.tools)
        
            start_time = time.time()

            if self.tools.args.mode == 'simulation':
                [probabilities_matrix, time_statevector] = qMetropolis.execute_quantum_metropolis_n()
            elif self.tools.args.mode == 'experiment':
                self.n_circuit_runs_to_average = self.tools.config_variables['n_circuit_runs_to_average']
                [probabilities_matrix, time_statevector] = qMetropolis.execute_real_hardware(step, self.n_circuit_runs_to_average)
            else:
                print("<*> ERROR!! Quantum execution mode not recognized. The mode selected is ", self.tools.args.mode)

            q_time = time.time() - start_time
            print("<i> QUANTUM METROPOLIS: Time for", step,"steps:", q_time, "seconds (", time_statevector,"seconds statevector)")
            q_tts = self.calculate_tts_from_probability_matrix(probabilities_matrix, index_min_energy, step, self.tools.config_variables['precision_solution'])

            ###### Classical Metropolis ######
            classical_metropolis = metropolis.Metropolis(step, self.n_angles, deltas_dict, self.tools)
            
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
            c_tts = self.calculate_tts_from_probability_matrix(probabilities_matrix, index_min_energy, step, self.tools.config_variables['precision_solution'])


            ###### Accumulated values Quantum Metropolis ######
            q_accumulated_tts.append(q_tts)     
            if q_tts < min_q_tts['value'] or min_q_tts['value'] == -1: min_q_tts.update(dict(value=q_tts, step=step))
        
            ###### Accumulated values Classical Metropolis ######
            c_accumulated_tts.append(c_tts)
            if c_tts < min_c_tts['value'] or min_c_tts['value'] == -1: min_c_tts.update(dict(value=c_tts, step=step))
            

            # plot data
            self.tools.plot_tts(q_accumulated_tts, c_accumulated_tts, self.tools.config_variables['initial_step'])

            # generate json data
            final_stats = {'q': min_q_tts, 'c': min_c_tts}
            self.tools.write_tts(self.tools.config_variables['initial_step'], self.tools.config_variables['final_step'], q_accumulated_tts, c_accumulated_tts, self.initialization_stats, final_stats)

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