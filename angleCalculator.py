import math
import metropolis
import quantumMetropolis
import time
import utils
import collections, functools, operator
import initializer

class AngleCalculator():

    def __init__(self, tools, angle_initializer, initialization_stats):

        self.tools = tools
        self.initializer = angle_initializer
        self.initialization_stats = initialization_stats

        self.mode = tools.args.mode
        self.initial_step = tools.config_variables['initial_step']
        self.final_step = tools.config_variables['final_step']
        self.precision_solution = tools.config_variables['precision_solution']
        self.default_value_tts = tools.config_variables['default_value_tts']

        self.n_angles = (len(self.tools.args.aminoacids) -1)*2

        self.quantum_simulation_activated = tools.config_variables['quantum_simulation_activated']

    def calculate3DStructure(self, deltas_dict, index_min_energy):

        q_accumulated_tts = []
        c_accumulated_tts = []

        min_q_tts = {'step': 0, 'value': -1}
        min_c_tts = {'step': 0, 'value': -1}

        classical_metropolis = metropolis.Metropolis(self.n_angles, deltas_dict, self.tools)

        #'''
        if self.quantum_simulation_activated == True:
            quantum_metropolis = quantumMetropolis.QuantumMetropolis(self.n_angles, deltas_dict, self.tools)
            ###### Quantum Metropolis ######
            start_time = time.time()

            if self.mode == 'simulation':
                [dict_probabilities_matrix, time_statevector] = quantum_metropolis.execute_quantum_metropolis_n(nW = self.final_step)
            elif self.mode == 'experiment': 
                [experiment_result_matrix, time_statevector, execution_stats, measures_dict] = quantum_metropolis.execute_real_hardware(nWs = 2)
            elif self.mode == 'real':
                n_repetitions = self.tools.config_variables['number_repetitions_real_mode']
                accum_probabilities = []
                for _ in range(n_repetitions):
                    [dict_probabilities_matrix, time_statevector] = quantum_metropolis.execute_quantum_metropolis_n(nW = self.tools.config_variables['w_real_mode'] + 1) # The 1 is needed because we implement a range (1, nW)
                    probabilities_matrix = dict_probabilities_matrix[self.tools.config_variables['w_real_mode']]
                    accum_probabilities.append(probabilities_matrix)

                real_q_counts = dict(functools.reduce(operator.add, map(collections.Counter, accum_probabilities)))
                real_q_counts = {k:v/n_repetitions for k,v in real_q_counts.items()}

            else:
                print("<*> ERROR!! Quantum execution mode not recognized. The mode selected is ", self.tools.args.mode)

            q_time = time.time() - start_time
            print("<i> QUANTUM METROPOLIS: Time for", self.initial_step,"steps:", q_time, "seconds (", time_statevector,"seconds statevector)")

            if self.mode == 'simulation':

                for step, probabilities_matrix in dict_probabilities_matrix.items():

                    ###### Accumulated values Quantum Metropolis ######
                    q_tts = self.calculate_tts_from_probability_matrix(probabilities_matrix, index_min_energy, step, self.precision_solution)

                    q_accumulated_tts.append(q_tts)     
                    if q_tts < min_q_tts['value'] or min_q_tts['value'] == -1: min_q_tts.update(dict(value=q_tts, step=step))
        #'''
        ###### Classical Metropolis ######
        for step in range(self.initial_step, self.final_step):  
            start_time = time.time()

            if self.mode == 'real':
                step = self.tools.config_variables['w_real_mode']

            probabilities_matrix = classical_metropolis.execute_metropolis(step)

            if self.mode == 'real':
                # get the key of the maximum value in the probability matrix
                index_min_energy = max(probabilities_matrix.items(), key=operator.itemgetter(1))[0]
                classical_selected_position = index_min_energy
                classical_confidence = probabilities_matrix[index_min_energy]

            print("<i> CLASSICAL METROPOLIS: Time for", step, "steps: %s seconds" % (time.time() - start_time))
            c_tts = self.calculate_tts_from_probability_matrix(probabilities_matrix, index_min_energy, step, self.precision_solution)

            #### create json files ####
            if self.mode == 'simulation':
            
                ###### Accumulated values Classical Metropolis ######
                c_accumulated_tts.append(c_tts)
                if c_tts < min_c_tts['value'] or min_c_tts['value'] == -1: min_c_tts.update(dict(value=c_tts, step=step))
                
                if step == self.final_step - 1: 
                    # plot data
                    #self.tools.plot_tts(q_accumulated_tts, c_accumulated_tts, self.tools.config_variables['initial_step'])

                    # generate json data
                    final_stats = {'q': min_q_tts, 'c': min_c_tts}
                    self.tools.write_tts(q_accumulated_tts, c_accumulated_tts, self.initialization_stats, final_stats)

            elif self.mode == 'experiment':

                p_t = experiment_result_matrix['betas=betas']['raw']['00']/self.tools.config_variables['ibmq_shots']
                q_tts = self.tools.calculateTTS(self.tools.config_variables['precision_solution'], step, p_t)
                
                ###### Accumulated values Quantum Metropolis ######
                q_accumulated_tts.append(q_tts)     
                if q_tts < min_q_tts['value'] or min_q_tts['value'] == -1: min_q_tts.update(dict(value=q_tts, step=step))
            
                ###### Accumulated values Classical Metropolis ######
                c_accumulated_tts.append(c_tts)
                if c_tts < min_c_tts['value'] or min_c_tts['value'] == -1: min_c_tts.update(dict(value=c_tts, step=step))

                self.tools.write_experiment_results(self.initialization_stats, experiment_result_matrix, execution_stats, measures_dict)

            # in real mode it is not necessary to execute the loop (there is only one step/w) so it breaks the loop
            if self.mode == 'real':

                quantum_success = False
                classical_success = False

                [quantum_selected_position, quantum_confidence] = self.get_selected_position_and_confidence(real_q_counts)
                [quantum_energy, quantum_configuration] = self.initializer.get_energy_configuration_from_position(quantum_selected_position, self.tools.args)

                if quantum_selected_position == index_min_energy: quantum_success = True 
                quantum_stats = {'confidence':quantum_confidence, 'success':quantum_success, 'energy':quantum_energy, 'configuration':quantum_configuration}
                
                # if the classical position is equal than the quantum, it is not necessary to recalculate the configuration and energy, it is the same
                if classical_selected_position == quantum_selected_position:
                    classical_energy = quantum_energy
                    classical_configuration = quantum_configuration
                else:
                    [classical_energy, classical_configuration] = self.initializer.get_energy_configuration_from_position(classical_selected_position, self.tools.args)
                
                if classical_selected_position == index_min_energy: classical_success = True 
                classical_stats = {'confidence':classical_confidence, 'success': classical_success, 'energy':classical_energy, 'configuration':classical_configuration}

                self.tools.write_real_results(self.initialization_stats, quantum_stats, classical_stats)

                min_q_tts['value'] = quantum_confidence
                min_q_tts['success'] = quantum_success
                min_c_tts['value'] = classical_confidence
                min_c_tts['success'] = classical_success

                break

        return [min_q_tts, min_c_tts]

    def get_selected_position_and_confidence(self, real_counts):

        max_value = 0
        position = 0
        confidence = 0

        for key in real_counts.keys():
            if real_counts[key] > max_value:
                max_value = real_counts[key]
                position = key
                confidence = real_counts[key]

        return [position, confidence]

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
            result = self.default_value_tts
        else:
            result = self.tools.calculateTTS(precision_solution, step, p_t)

        return result