import time
import collections, functools, operator
import json
import numpy as np
import math

from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, Aer, IBMQ
from qiskit.compiler import transpile


# Import measurement calibration functions
import scipy

class ExecutorIBMQ:

    def __init__(self, tools, input_oracle):

        self.tools = tools
        self.input_oracle = input_oracle

        self.qiskit_api_path = self.tools.config_variables['path_qiskit_token']
        self.selected_device = self.tools.config_variables['device_ibm_q']

        self.login_ibmq()

    def execute_real_hardware(self, nWs):

        start_time = time.time()
        shots = self.tools.config_variables['ibmq_shots']
        n_repetitions = self.tools.config_variables['number_repetitions_ibmq']
        n_repetitions_zero_beta = self.tools.config_variables['number_repetitions_ibmq_zero_beta']

        # prepare dictionary with deltas
        deltas_dictionary = collections.OrderedDict(sorted(self.input_oracle.items()))
        deltas = {}
        for (key,value) in deltas_dictionary.items():
            deltas[key[:3]] = value
        
        counts = {}
        measures_dict = {}

        # First we load all the previous results so that for beta = 0 we do not have to recalculate more than necessary
        with open('./results/measurements.json', 'r') as outfile2: 
            dictionary = json.load(outfile2)

        try: # the dictionary has the form dictionary['--']['0-0']['measurements'] = {'00': [1329,3213 ...], '01':...}
            beta0_counts = dictionary['--']['0-0']['measurements'] 
            len_beta0_counts00  = len(beta0_counts['00'])

            if len_beta0_counts00 >= n_repetitions_zero_beta:
                measures_dict['0-0'] = beta0_counts
                runs = [1]

            else:
                beta0_n_repetitions = n_repetitions_zero_beta - len_beta0_counts00
                runs = range(2)

        except:
            beta0_counts = {'00': [], '01':[], '10':[], '11':[]}
            beta0_n_repetitions = n_repetitions_zero_beta
            runs = range(2)

        
        # Then we execute the needed runs
        for index in runs:

            # in the first iteration (index=0) it uses the betas = 0. In the second iteration, it uses the betas of the config file
            if index == 0:
                betas = [1e-10,1e-10]
                key_name_counts = 'betas=0'
                reps = beta0_n_repetitions  
            
            else:
                betas = self.tools.config_variables['betas']
                key_name_counts = 'betas=betas'
                reps = n_repetitions

            counts[key_name_counts] = {}
            # Let us first analyse the noise of the circuit for the ideal case of betas = 0, which should imply .25 chance of success
            qc = self.generate_circ(nWs, deltas, betas)
            

            # get the NOISELESS counts
            counts[key_name_counts]['noiseless'] = self.exe_noiseless(nWs)

            # get the RAW counts
            raw_counts = []
            for i in range(reps):
                print("<i> Waiting to get access to IBMQ processor. Betas = ", betas, ". Iteration = ",i)
                #raw_counts.append(execute(qc, backend, shots=shots).result().get_counts())
                raw_counts.append(execute(qc, Aer.get_backend('qasm_simulator'), shots=shots).result().get_counts())
                print("<i> Circuit in IBMQ executed")

            if index == 0: # Notice that we will add here the measurements for beta =0 already saved in results.json
                measures_dict['0-0']= self.tools.list_of_dict_2_dict_of_lists(raw_counts, beta0_counts = beta0_counts)
            else:
                measures_dict[str(betas[0]) + '-' +str(betas[1])]= self.tools.list_of_dict_2_dict_of_lists(raw_counts)

            # sum all values of the same position and get the mean of each position to store in counts
            raw_counts = dict(functools.reduce(operator.add, map(collections.Counter, raw_counts)))
            raw_counts = {k:v/n_repetitions for k,v in raw_counts.items()}
            counts[key_name_counts]['raw'] = raw_counts


        # In order to see if there is some statistical difference between the two noise circuit (due to the value of beta and the angles)
        # we generate bernouilli distribuitions that follow the same statistics as those that we have measured
        betas = self.tools.config_variables['betas']
        print('measures_dict',measures_dict)
        beta0_bernouilli = self.generate_bernouilli(int(sum(measures_dict['0-0']['00'])), shots*len(measures_dict['0-0']['00']))
        beta1_bernouilli = self.generate_bernouilli(int(sum(measures_dict[str(betas[0]) + '-' +str(betas[1])]['00'])), shots*len(measures_dict[str(betas[0]) + '-' +str(betas[1])]['00']))
        exec_stats, pvalue = scipy.stats.ttest_ind(beta0_bernouilli, beta1_bernouilli, equal_var=False)

        execution_stats = 'The t-test statistic value for there being a significat average difference between measured processes with beta zero and non-zero is ' + str(exec_stats) + ' and the corresponding pvalue is '+ str(pvalue)
        print('<i>', execution_stats)

        time_statevector = time.time() - start_time

        return [counts, time_statevector, execution_stats, measures_dict]

    def login_ibmq(self):

        #read the file that contains the Qiskit user API
        with open(self.qiskit_api_path) as json_file:
            api_token = json.load(json_file)['qiskit_token']

        if api_token == '':
            print('<*> ERROR!! It is necessary to introduce an qiskit API token')

        IBMQ.save_account(api_token, overwrite=True)
        IBMQ.load_account()

        provider = IBMQ.get_provider(hub=self.tools.config_variables['hub'],
                                        group=self.tools.config_variables['group'], 
                                        project=self.tools.config_variables['project'])

        self.backend = provider.get_backend(self.selected_device)
        
        return self.backend

    def generate_circ(self, nWs, deltas, betas):

        assert(len(betas) == nWs)

        move_id  = QuantumRegister(1)
        angle_phi = QuantumRegister(1)
        angle_psi = QuantumRegister(1)
        coin = QuantumRegister(1)
        c_reg = ClassicalRegister(2)
        qc = QuantumCircuit(coin,move_id,angle_psi,angle_phi,c_reg)

        #Circuit ----------
        qc.h(angle_phi)
        qc.h(angle_psi)
        for (i,beta) in zip(range(nWs),betas):
            angles = self.calculate_angles(deltas, beta)
            self.W_step(qc,coin,move_id,angle_psi,angle_phi,angles,nW = i, nWs = nWs)

        # Measure
        qc.measure(angle_phi[0], c_reg[1])
        qc.measure(angle_psi[0], c_reg[0])

        # Transpiling -------

        #layout = {5: angle_phi[0], 6: angle_psi[0], 4: move_id[0], 5: coin[0]}
        layout = {2: angle_psi[0], 3: angle_phi[0], 1: coin[0], 0: move_id[0]} 
        qc = transpile(qc, backend = self.backend, optimization_level=3, 
                    initial_layout=layout, basis_gates = ['u1', 'u2', 'u3', 'cx'], routing_method = 'lookahead')
        
        print('\n⬤⬤⬤⬤  Circuit stadistics after optimization  ⬤⬤⬤⬤\n')
        print('•  Gates = ', qc.count_ops())
        print('•  Depth = ', qc.depth())
        print('\n') 
        
        return qc

    def generate_bernouilli(self, n_0, n):
        array = np.random.binomial(1, n_0/n, n)
        s = np.sum(array)
        while s != n_0:
            i = np.random.randint(n)
            if s<n_0 and array[i] == 0:
                array[i] = 1
            elif s>n_0 and array[i] == 1:
                array[i] = 0
            s = np.sum(array)
        return array

    def exe_noiseless(self, nWs):

        betas = self.tools.config_variables['betas']

        # prepare dictionary with deltas
        deltas_dictionary = collections.OrderedDict(sorted(self.input_oracle.items()))
        deltas = {}
        for (key,value) in deltas_dictionary.items():
            deltas[key[:3]] = value

        print('deltas', deltas)

        aerqc = self.generate_hardware_simulation_circuit(nWs, deltas, betas)

        aerbackend = Aer.get_backend('statevector_simulator')
        backend_options = {"method" : "statevector"}
        experiment = execute(aerqc, aerbackend, backend_options=backend_options)
        state_vector = Statevector(experiment.result().get_statevector(aerqc))

        probabilities = state_vector.probabilities([2,3]) # We are reporting the angles as (psi,phi); since qiskit inverts the reporting order
        print('probabilities',probabilities)
        noiseless_counts = {}
        noiseless_counts['00'] = float(probabilities[0])
        noiseless_counts['01'] = float(probabilities[1])
        noiseless_counts['10'] = float(probabilities[2])
        noiseless_counts['11'] = float(probabilities[3])

        return noiseless_counts

    def generate_hardware_simulation_circuit(self, nWs, deltas, betas):

        assert(len(betas) == nWs)

        move_id  = QuantumRegister(1)
        angle_phi = QuantumRegister(1)
        angle_psi = QuantumRegister(1)
        coin = QuantumRegister(1)
        c_reg = ClassicalRegister(2)
        aerqc = QuantumCircuit(coin,move_id,angle_psi,angle_phi,c_reg)

        #Circuit ----------
        aerqc.h(angle_phi)
        aerqc.h(angle_psi)
        for (i,beta) in zip(range(nWs),betas):
            angles = self.calculate_angles(deltas, beta)
            print('angles step',i,angles)
            self.simulated_W_step(aerqc,coin,move_id,angle_psi,angle_phi,angles,nW = i, nWs = nWs)
        
        return aerqc

    def simulated_W_step(self, qc,coin,move_id,angle_psi,angle_phi,angles,nW,nWs):

        # Perform the preparation of possible moves----
        qc.h(move_id)

        # Prepare the Boltzmann coin ------------------
        self.simulated_hardware_1_coin_flip(qc, coin, move_id, angle_psi, angle_phi, angles, inverse = 1)
        
        # Perform move ---------------------------------
        # For the second angle
        qc.ccx(coin,move_id,angle_psi)

        # For the first angle
        qc.x(move_id)
        qc.ccx(coin,move_id,angle_phi)
        qc.x(move_id)

        if nW < nWs-1: # This happens unless we are in the last step, in which case uncomputing is unnecessary.
            # Unprepare the Boltzmann coin--------------------
            self.simulated_hardware_1_coin_flip(qc, coin, move_id, angle_psi, angle_phi, angles, inverse = -1)

            # Perform the preparation of possible moves ----
            qc.h(move_id)

            #Reflection -------------------------------------
            qc.x(move_id)
            qc.x(coin)

            # Perform a multicontrolled Z
            qc.cz(move_id,coin)

            qc.x(move_id)
            qc.x(coin)

    def W_step(self, qc,coin,move_id,angle_psi,angle_phi,angles,nW,nWs): 
            
        # Perform the preparation of possible moves----
        qc.h(move_id)

        # Prepare the Boltzmann coin ------------------
        self.hardware_1_coin_flip(qc, coin, move_id, angle_psi, angle_phi, angles, inv = 1)
        
        # Perform move ---------------------------------
        # For the second angle
        qc.ccx(coin,move_id,angle_psi)

        # For the first angle
        qc.x(move_id)
        qc.ccx(coin,move_id,angle_phi)
        qc.x(move_id)

        if nW < nWs-1: # This happens unless we are in the last step, in which case uncomputing is unnecessary.
            # Unprepare the Boltzmann coin--------------------
            self.hardware_1_coin_flip(qc, coin, move_id, angle_psi, angle_phi, angles, inv = -1)

            # Perform the preparation of possible moves ----
            qc.h(move_id)

            #Reflection -------------------------------------
            qc.x(move_id)
            qc.x(coin)

            # Perform a multicontrolled Z
            qc.cz(move_id,coin)

            qc.x(move_id)
            qc.x(coin)

    def simulated_hardware_1_coin_flip(self, circuit, coin, move_id, angle_psi, angle_phi, angles, inverse):
        ''' Applies the controlled rotation to the target coin'''
        if inverse == 1:
            circuit.x(coin)
        
        if angles['111'] > .01:
            circuit.mcrx(theta = -inverse * angles['111'], q_controls = [angle_phi[0],angle_psi[0],move_id[0]], q_target = coin[0], use_basis_gates=False)
        circuit.x(move_id)
        
        if angles['110'] > .01:
            circuit.mcrx(theta = -inverse * angles['110'], q_controls = [angle_phi[0],angle_psi[0],move_id[0]], q_target = coin[0], use_basis_gates=False)
        circuit.x(angle_psi)
        
        if angles['100'] > .01:
            circuit.mcrx(theta = -inverse * angles['100'], q_controls = [angle_phi[0],angle_psi[0],move_id[0]], q_target = coin[0], use_basis_gates=False)
        circuit.x(move_id)
        
        if angles['101'] > .01:
            circuit.mcrx(theta = -inverse * angles['101'], q_controls = [angle_phi[0],angle_psi[0],move_id[0]], q_target = coin[0], use_basis_gates=False)
        circuit.x(angle_phi)
        
        if angles['001'] > .01:
            circuit.mcrx(theta = -inverse * angles['001'], q_controls = [angle_phi[0],angle_psi[0],move_id[0]], q_target = coin[0], use_basis_gates=False)
        circuit.x(move_id)
        
        if angles['000'] > .01:
            circuit.mcrx(theta = -inverse * angles['000'], q_controls = [angle_phi[0],angle_psi[0],move_id[0]], q_target = coin[0], use_basis_gates=False)
        circuit.x(angle_psi)
        
        if angles['010'] > .01:
            circuit.mcrx(theta = -inverse * angles['010'], q_controls = [angle_phi[0],angle_psi[0],move_id[0]], q_target = coin[0], use_basis_gates=False)
        circuit.x(move_id)
        
        if angles['011'] > .01:
            circuit.mcrx(theta = -inverse * angles['011'], q_controls = [angle_phi[0],angle_psi[0],move_id[0]], q_target = coin[0], use_basis_gates=False) 
        circuit.x(angle_phi)
        
        if inverse == -1:
            circuit.x(coin)

    def hardware_1_coin_flip(self, circuit, coin, move_id, angle_psi, angle_phi, angles, inv):

        '''Warning! This only works for dipeptide 1 in experiment mode. Do not use elsewhere!'''
        # First we have to identify the non-zero angles. For the rest we accept with probability 1
        circuit.x(coin)
        '''
        Since the angles from 001 and 101 ~= 2.59; and those from 010 and 000 ~= 0.32 (when beta = .1, 
        but they'll always be similar nevertheless), we will perform those rotations together
        '''
        non_zero_angles = {}
        non_zero_angles['0x0'] = (angles['000']+angles['010'])/2
        non_zero_angles['x01'] = (angles['001']+angles['101'])/2
        
        # Let us first perform the first
        circuit.x(angle_phi)
        circuit.x(move_id)
        circuit.mcrx(theta = -inv*non_zero_angles['0x0'],
                    q_controls = [move_id[0],angle_phi[0]], q_target = coin[0], use_basis_gates=True)
        circuit.x(angle_phi)
        circuit.x(move_id)
        
        # Let us perform the second
        circuit.x(angle_psi)
        circuit.mcrx(theta = -inv*non_zero_angles['x01'],
                    q_controls = [move_id[0],angle_psi[0]], q_target = coin[0], use_basis_gates=True)
        circuit.x(angle_psi)

    def calculate_angles(self, deltas_dictionary, beta):

        exact_angles = {}

        for key in deltas_dictionary.keys():

            if deltas_dictionary[key] >= 0:
                probability = math.exp(-beta * deltas_dictionary[key])
            else: 
                probability = 1
            # Instead of encoding the angle corresponding to the probability, we will encode the angle theta such that sin^2(pi/2 - theta) = probability.
            # That way 1 -> 000, but if probability is 0 there is some small probability of acceptance

            # Instead of probability save angles so rotations are easier to perform afterwards sqrt(p) = sin(pi/2-theta/2).
            # The theta/2 is because if you input theta, qiskits rotates theta/2. Also normalised (divided between pi the result)
            exact_angles[key] = math.pi - 2 * math.asin(math.sqrt(probability))

        # Order angles by key
        exact_angles = collections.OrderedDict(sorted(exact_angles.items()))

        return exact_angles

