import math
import metropolis
import quantumMetropolis
import time
import utils
import collections, functools, operator 

class AngleCalculator():

    def __init__(self, tools, initialization_stats):

        self.tools = tools
        self.initialization_stats = initialization_stats

        self.n_angles = (len(self.tools.args.aminoacids) -1)*2

    def calculate3DStructure(self, deltas_dict, index_min_energy):

        q_accumulated_tts = []
        c_accumulated_tts = []

        min_q_tts = {'step': 0, 'value': -1}
        min_c_tts = {'step': 0, 'value': -1}

        quantum_metropolis = quantumMetropolis.QuantumMetropolis(self.n_angles, deltas_dict, self.tools)
        classical_metropolis = metropolis.Metropolis(self.n_angles, deltas_dict, self.tools)
        for step in range(self.tools.config_variables['initial_step'], self.tools.config_variables['final_step']):

            ###### Quantum Metropolis ######
            start_time = time.time()

            if self.tools.args.mode == 'simulation':
                [probabilities_matrix, time_statevector] = quantum_metropolis.execute_quantum_metropolis_n(step)
            elif self.tools.args.mode == 'experiment':
                [experiment_result_matrix, time_statevector, execution_stats] = quantum_metropolis.execute_real_hardware(step)
            elif self.tools.args.mode == 'real':
                n_repetitions = self.tools.config_variables['number_repetitions_real_mode']
                accum_probabilities = []
                for _ in range(n_repetitions):
                    [probabilities_matrix, time_statevector] = quantum_metropolis.execute_quantum_metropolis_n(self.tools.config_variables['w_real_mode'])
                    accum_probabilities.append(probabilities_matrix)

                real_q_counts = dict(functools.reduce(operator.add, map(collections.Counter, accum_probabilities)))
                real_q_counts = {k:v/n_repetitions for k,v in real_q_counts.items()}

            else:
                print("<*> ERROR!! Quantum execution mode not recognized. The mode selected is ", self.tools.args.mode)

            q_time = time.time() - start_time
            print("<i> QUANTUM METROPOLIS: Time for", step,"steps:", q_time, "seconds (", time_statevector,"seconds statevector)")

            if self.tools.args.mode == 'simulation':
                q_tts = self.calculate_tts_from_probability_matrix(probabilities_matrix, index_min_energy, step, self.tools.config_variables['precision_solution'])

            ###### Classical Metropolis ######
            start_time = time.time()

            if self.tools.args.mode == 'real':
                n_repetitions = self.tools.config_variables['number_repetitions_real_mode']
                accum_probabilities = []
                for _ in range(n_repetitions):
                    probabilities_matrix = classical_metropolis.execute_metropolis(self.tools.config_variables['w_real_mode'])
                    accum_probabilities.append(probabilities_matrix)
    
                real_c_counts = dict(functools.reduce(operator.add, map(collections.Counter, accum_probabilities)))
                real_c_counts = {k:v/n_repetitions for k,v in real_c_counts.items()}

            else:
                probabilities_matrix = classical_metropolis.execute_metropolis(step)

            print("<i> CLASSICAL METROPOLIS: Time for", step, "steps: %s seconds" % (time.time() - start_time))
            if self.tools.args.mode == 'simulation' or self.tools.args.mode == 'simulation':
                c_tts = self.calculate_tts_from_probability_matrix(probabilities_matrix, index_min_energy, step, self.tools.config_variables['precision_solution'])


            #### create json files ####
            if self.tools.args.mode == 'simulation':
                
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
                self.tools.write_tts(q_accumulated_tts, c_accumulated_tts, self.initialization_stats, final_stats)

            elif self.tools.args.mode == 'experiment':

                p_t = experiment_result_matrix['betas=betas']['raw']['00']/self.tools.config_variables['ibmq_shots']
                q_tts = self.tools.calculateTTS(self.tools.config_variables['precision_solution'], step, p_t)
                
                ###### Accumulated values Quantum Metropolis ######
                q_accumulated_tts.append(q_tts)     
                if q_tts < min_q_tts['value'] or min_q_tts['value'] == -1: min_q_tts.update(dict(value=q_tts, step=step))
            
                ###### Accumulated values Classical Metropolis ######
                c_accumulated_tts.append(c_tts)
                if c_tts < min_c_tts['value'] or min_c_tts['value'] == -1: min_c_tts.update(dict(value=c_tts, step=step))

                self.tools.write_experiment_results(self.initialization_stats, experiment_result_matrix, execution_stats)

            # in real mode it is not necessary to execute the loop (there is only one step/w) so it breaks the loop
            if self.tools.args.mode == 'real':
                
                self.tools.write_real_results(self.initialization_stats, real_q_counts, real_c_counts)
                break
        
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