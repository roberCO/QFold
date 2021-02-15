import math
import time
import collections, functools, operator

import initializer
import executorIBMQ
import qms

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
        self.classical_simulation_activated = tools.config_variables['classical_simulation_activated']

    def calculate3DStructure(self, energies_dict, index_min_energy):

        q_accumulated_tts = []
        c_accumulated_tts = []

        min_q_tts = {'step': 0, 'value': -1}
        min_c_tts = {'step': 0, 'value': -1}

        # create instance QMS
        qms_solver = qms.QMS(energies_dict, index_min_energy)

        ###### Quantum Metropolis ######
        if self.quantum_simulation_activated == True:
            start_time = time.time()

            if self.mode == 'simulation':
                q_tts_all_steps = qms_solver.execute_quantum_metropolis(nW = self.final_step, 'TTS')
            elif self.mode == 'experiment': 
                ibmq_executor = executorIBMQ.ExecutorIBMQ(self.tools, energies_dict)
                [experiment_result_matrix, time_statevector, execution_stats, measures_dict] = ibmq_executor.execute_real_hardware(nWs = 2)
            elif self.mode == 'real':
                n_repetitions = self.tools.config_variables['number_repetitions_real_mode']
                accum_probabilities = []
                for _ in range(n_repetitions):
                    [dict_probabilities_matrix, time_statevector] = qms_solver.execute_quantum_metropolis(nW = self.tools.config_variables['w_real_mode'] + 1, 'probabities') # The 1 is needed because we implement a range (1, nW)
                    probabilities_matrix = dict_probabilities_matrix[self.tools.config_variables['w_real_mode']]
                    accum_probabilities.append(probabilities_matrix)

                real_q_counts = dict(functools.reduce(operator.add, map(collections.Counter, accum_probabilities)))
                real_q_counts = {k:v/n_repetitions for k,v in real_q_counts.items()}

            else:
                print("<*> ERROR!! Quantum execution mode not recognized. The mode selected is ", self.tools.args.mode)

            for step, q_tts in q_tts_all_steps.items():
                if q_tts < min_q_tts['value'] or min_q_tts['value'] == -1: min_q_tts.update(dict(value=q_tts, step=step))
            
            q_time = time.time() - start_time
            print("<i> QUANTUM METROPOLIS: Time for", self.final_step,"steps:", q_time, "seconds")

        ###### Classical Metropolis ######
        if self.classical_simulation_activated == True:

            start_time = time.time()
            if self.mode == 'real':
                probabilities_dict = qms.execute_classical_metropolis(self.tools.config_variables['w_real_mode'], 'probabilities')
            elif self.mode == 'simulation':
                c_tts_all_steps = qms.execute_classical_metropolis(self.final_step, 'TTS')

            for step, q_tts in q_tts_all_steps.items():
                if q_tts < min_q_tts['value'] or min_q_tts['value'] == -1: min_q_tts.update(dict(value=q_tts, step=step))
            
            c_time = time.time() - start_time
            print("<i> CLASSICAL METROPOLIS: Time for", self.final_step,"steps:", c_time, "seconds")

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