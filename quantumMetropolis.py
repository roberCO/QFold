import numpy as np
from itertools import product
import time
import math
from math import pi

# Importing standard Qiskit libraries and configuring account
import qiskit
from qiskit import QuantumCircuit, execute, Aer
from qiskit.circuit import QuantumRegister, Qubit, Gate
from qiskit.aqua.components.oracles import Oracle, TruthTableOracle
from qiskit.quantum_info import Statevector

import beta_precalc_TruthTableOracle


class QuantumMetropolis():

    def __init__(self, n_repetitions, angle_precision_bits, probability_bits, n_angles, beta_max, scaling_factor, input_oracle):

        #Global variables

        # Number steps
        self.n_repetitions = n_repetitions
        self.scaling_factor = scaling_factor

        # Number of bits necessary to specify the position of each angle
        self.angle_precision_bits = angle_precision_bits

        #Oracle output ancilla bits
        self.probability_bits = probability_bits

        self.n_angles = n_angles
        self.beta_max = beta_max

        self.move_id_len = int(np.ceil(np.log2(n_angles)))
        self.n_ancilla_bits = self.probability_bits

        if self.n_ancilla_bits < 3:
            raise ValueError('The minimum number of ancilla qubits needed for this algorithm is 3! Currently there are only', self.n_ancilla_bits) 

        if self.n_angles*self.angle_precision_bits + self.move_id_len + self.n_ancilla_bits + 2 > 32:
            raise ValueError('The number of qubits is too large (larger than 32)! Currently there are\n'+
                            str(self.n_angles)+ ' angles, each with '+str(self.angle_precision_bits)+' qubits\n'+
                            'an ancilla with '+str(self.n_ancilla_bits)+' bits\n'+
                            'a move_id register with '+str(self.move_id_len)+' bits\n'+
                            'and finally a single qubit called move_value and another for the coin.\n'+
                            'In total there are '+str(self.n_angles*self.angle_precision_bits + self.move_id_len + self.n_ancilla_bits + 1)+' qubits\n'
                            )

        # The delta E's dictionary
        self.input_oracle = input_oracle

        # For n angles
        [self.move_preparation_gate, self.conditional_move_gate_n, self.reflection_gate] = self.prepare_initial_circuits_n()

    def move_preparation(self, circuit,move_value,move_id):
        '''
        Proposes new moves
        '''
        circuit.h(move_value) #Es un único bit que puede ser 0 para + y 1 para -
        if bin(self.n_angles).count('1') == 1: # if the number of angles is a power of two
            circuit.h(move_id)
        else:
            vector = ([1]*self.n_angles + [0]*(2**(self.move_id_len) - self.n_angles))/np.sqrt(self.n_angles)
            circuit.initialize(vector, [move_id[j] for j in range(self.move_id_len)])

    def conditional_move_npeptide(self,circuit,ancilla,coin,move_value,move_id,angles):
        '''
        Conditioned on coin, perform a move. Tested ok.
        We use a repetitive structure where we perform the conditional sum and subtraction for each angle.
        '''
        # For each angle
        for i in range(self.n_angles):
            angle = angles[i] #Select the angle from the list of registers
            angle_index = np.binary_repr(i, width=self.move_id_len) #convert i to binary

            # Put the given move_id in all 1 to control on it: for instance if we are controling on i=2, move 010 ->111
            for j in range(len(angle_index)):
                if angle_index[j] == '0':
                    circuit.x(move_id[j])

            circuit.mcx(control_qubits= [coin[0]]+[move_id[j] for j in range(move_id.size)], target_qubit = ancilla[0])#create a single control
            self.sumsubtract1(circuit,angle,ancilla[0],ancilla[1],ancilla[2],move_value[0]) #sum or subtract 1 to the angle
            circuit.mcx(control_qubits= [coin[0]]+[move_id[j] for j in range(move_id.size)], target_qubit = ancilla[0])#create a single control        
            
            # Undo the move_id preparation: for instance, if we are controlling on i= 2 move 111->010
            for j in range(len(angle_index)):
                if angle_index[j] == '0':
                    circuit.x(move_id[j])

    def reflection(self, circuit,coin,move_value,move_id):
        '''
        I have to investigate over what is the reflection performed. Is it performed over 000?
        If in state 0000, make it 1111, cccz gate and back to 0000
        '''
        circuit.x(move_id)
        circuit.x(move_value)
        circuit.x(coin)
        
        # Perform a multicontrolled Z
        circuit.h(coin[0])
        circuit.mcx(control_qubits = [move_id[j] for j in range(self.move_id_len)]+ [move_value[0]], target_qubit = coin[0])
        circuit.h(coin[0])
        
        circuit.x(move_id)
        circuit.x(move_value)
        circuit.x(coin)

    def prepare_initial_circuits_n(self):

        # Move preparation gate ---------------------------------------
        s_move_id = QuantumRegister(self.move_id_len) 
        s_move_value = QuantumRegister(1)

        sub_circ = QuantumCircuit(s_move_value,s_move_id)
        self.move_preparation(sub_circ,s_move_value,s_move_id)
        move_preparation_gate = sub_circ.to_instruction()

        # Conditional move gate ----------------------------------------
        s_angles = []
        for i in range(self.n_angles):
            s_angles.append(QuantumRegister(self.angle_precision_bits, name = 'angle' + str(i)))
        s_move_id = QuantumRegister(self.move_id_len)
        s_move_value = QuantumRegister(1)
        s_coin = QuantumRegister(1)
        s_ancilla = QuantumRegister(self.probability_bits)


        # Creates the circuit
        sub_circ = QuantumCircuit(s_ancilla, s_coin, s_move_value, s_move_id)
        for i in range(self.n_angles-1,-1,-1):
            sub_circ = sub_circ + QuantumCircuit(s_angles[i])

        self.conditional_move_npeptide(sub_circ,s_ancilla, s_coin, s_move_value, s_move_id, s_angles)

        # Optimize the circuit

        '''
        print('Before optimization------- conditional_move_npeptide')
        print('gates = ', sub_circ.count_ops())
        print('depth = ', sub_circ.depth())
        sub_circ = transpile(sub_circ, seed_transpiler=1, optimization_level=3)
        print('After optimization--------')
        print('gates = ', sub_circ.count_ops())
        print('depth = ', sub_circ.depth())
        '''
        conditional_move_gate_n = sub_circ.to_instruction()

        # Reflection gate --------------------------------------------------
        s_move_id = QuantumRegister(self.move_id_len) 
        s_move_value = QuantumRegister(1)
        s_coin = QuantumRegister(1)

        sub_circ = QuantumCircuit(s_coin, s_move_value, s_move_id)
        self.reflection(sub_circ,s_coin, s_move_value, s_move_id)

        # We could optimize the circuit, but it will probably not be worth it (fairly small)

        reflection_gate = sub_circ.to_instruction()

        return [move_preparation_gate, conditional_move_gate_n, reflection_gate]

    def sum1(self, circuit,qubit_string,control,start,end):
        
        circuit.cx(control,end) # iff control = 1, end = 1
        circuit.x(start)
        circuit.cx(control,start) # iff control = 1, start = 0

        for i in range(qubit_string.size,-1,-1):
            '''
            Next thing we analise if all qubits to the right have value 1, 
            and save it in the current qubit and start.
            Don't need to add control, since end already does that work
            '''
            if i < qubit_string.size:
                # For i = 0, there is only the start to worry about
                circuit.mcx(control_qubits= [qubit_string[j] for j in range(i-1,-1,-1)]+[end], target_qubit = qubit_string[i])
            circuit.mcx(control_qubits = [qubit_string[j] for j in range(i-1,-1,-1)]+[end], target_qubit = start)

            '''
            Next, controlling on the current qubit and start, we change all the following qubits to 0.
            We have to control with qubit_string[n_qubit], start and control because start could be at state 1 without control also being in that state.
            '''
            if i == qubit_string.size:
                for j in range(i-1,-1,-1):
                    circuit.ccx(control,start,qubit_string[j])
                circuit.ccx(control,start,end)
            elif i == 0:
                circuit.mcx(control_qubits = [control,qubit_string[i],start], target_qubit = end)
            else:
                for j in range(i-1,-1,-1):            
                    circuit.mcx(control_qubits = [control,qubit_string[i],start], target_qubit = qubit_string[j])
                circuit.mcx(control_qubits = [control,qubit_string[i],start], target_qubit = end)
        circuit.x(start)

    def subtract1(self, circuit,qubit_string,control,start,end):
        '''
        Outputs:
        subtracts register 2 (1 qubit) from register 1 in register 1. Tested ok.
        
        Input:
        circuit: QuantumCircuit with registers qubit_string, control, ancilla
        
        qubit_string: QuantumRegister
        
        control: Qubit. Use ancilla[0] or similar
        
        start: Qubit. Use ancilla[1] or similar
        end: Qubit. Use ancilla[2] or similar
        
        Comments: In binary, subtracting is the same procedure as summing when we exchange 0s and 1s
        '''
        circuit.x(qubit_string)

        self.sum1(circuit,qubit_string,control,start,end)
        
        circuit.x(qubit_string)

    def sumsubtract1(self,circuit,qubit_string,control,start,end,move_value):
        '''
        Outputs:
        Sum/subtracts register 2 (control, 1 qubit) from register 1 (qubit_string) in register 1. Tested ok.

        Input:
        circuit: QuantumCircuit with registers qubit_string, control, ancilla and move_value
        qubit_string: QuantumRegister where the sum/subtraction is performed
        control: Qubit. Use ancilla[0] or similar. It encodes the probability of change.
        start: Qubit. Use ancilla[1] or similar
        end: Qubit. Use ancilla[2] or similar
        move_value: 1 to subtract, 0 to sum

        Comments: In binary, subtracting is the same procedure as summing when we exchange 0s and 1s
        '''
        circuit.cx(move_value,qubit_string)

        self.sum1(circuit,qubit_string,control,start,end)

        circuit.cx(move_value,qubit_string)

    def coin_flip(self, circuit,ancilla,coin):
        '''
        Prepares the coin with the probability encoded in the ancilla.
        The important thing to notice is that we are using the same convention as qiskit: littleendian.
        That means that the larger the index of the ancilla bit, the more significant it is, and the larger the rotation
        '''
        
        #Necesitamos usar el número guardado en las ancillas para realizar rotaciones controladas.  
        #Notice that ancilla encodes 1-probability, rather than probability.
        #Notice also that cu3(theta) rotates theta/2. As the first angle to rotate is pi/4 we need to start in theta = pi/2

        circuit.x(coin) # Start in 1 and decrease it, since we encoded the angle corresponding 1-probability
        for i in range(ancilla.size-1,-1,-1): # See how to perform an rx rotation in https://qiskit.org/documentation/stubs/qiskit.circuit.library.U3Gate.html
            circuit.cu3(theta = -math.pi*2**(i-ancilla.size), phi  = 0, lam = 0, control_qubit = ancilla[i], target_qubit = coin)
    
    def coin_flip_func_n(self, oracle):
        
        '''
        Defines de coin_flip_gate using the oracle that is provided on the moment.
        Notice that oracle gate has registers oracle.variable_register and oracle.output_register in that order
        oracle.variable_register should have size angle_phi.size + angle_psi.size + move_id.size + move_value.size
        oracle.output_register should have size self.probability_bits
        '''

        # Construct an instruction for the oracle
        oracle.construct_circuit()
        oracle_circuit = oracle.circuit

        #print(oracle_circuit)
        oracle_gate = oracle_circuit.to_instruction()

        # Let us create a circuit for coin_flip
        cf_angles = []
        for i in range(self.n_angles):
            cf_angles.append(QuantumRegister(self.angle_precision_bits, name = 'angle' + str(i)))
        cf_move_id = QuantumRegister(self.move_id_len) 
        cf_move_value = QuantumRegister(1)
        cf_coin = QuantumRegister(1)
        cf_ancilla = QuantumRegister(self.probability_bits)

        cf_circ = QuantumCircuit(cf_ancilla,cf_coin,cf_move_value,cf_move_id)
        for i in range(self.n_angles-1,-1,-1):
            cf_circ = cf_circ + QuantumCircuit(cf_angles[i])


        # Main operations
        cf_circ.append(oracle_gate, [cf_move_value[0]]+[cf_move_id[j] for j in range(cf_move_id.size)] + [cf_angles[k][j] for (k,j) in product(range(self.n_angles-1,-1,-1), range(self.angle_precision_bits))] + [cf_ancilla[j] for j in range(self.n_ancilla_bits)])
        self.coin_flip(cf_circ,cf_ancilla,cf_coin)
        cf_circ.append(oracle_gate.inverse(), [cf_move_value[0]]+[cf_move_id[j] for j in range(cf_move_id.size)]+[cf_angles[k][j] for (k,j) in product(range(self.n_angles-1,-1,-1), range(self.angle_precision_bits))]+ [cf_ancilla[j] for j in range(self.n_ancilla_bits)])
        
        coin_flip_gate = cf_circ.to_instruction()
        
        return coin_flip_gate

    def W_func_n(self, oracle):
        
        '''This defines the parametrised gate W using the oracle that is provided to it, and we can reuse its inverse too.'''

        # State definition. All angles range from 0 to 2pi
        w_angles = []
        for i in range(self.n_angles):
            w_angles.append(QuantumRegister(self.angle_precision_bits, name = 'angle' + str(i)))

        # Move proposal
        w_move_id = QuantumRegister(self.move_id_len, name = 'move_id') #Which angle are we modifying
        w_move_value = QuantumRegister(1, name = 'move_value') #0 -> decrease the angle. 1-> increase it

        # Coin
        w_coin = QuantumRegister(1, name = 'coin')

        # Ancillas
        w_ancilla = QuantumRegister(self.probability_bits, name = 'ancilla')

        # Circuit
        qc = QuantumCircuit(w_ancilla,w_coin,w_move_value,w_move_id)
        for i in range(self.n_angles-1,-1,-1):
            qc = qc + QuantumCircuit(w_angles[i])

        
        # Define the coin_flip_gate
        coin_flip_gate = self.coin_flip_func_n(oracle)

        # Move preparation
        qc.append(self.move_preparation_gate, [w_move_value[0]]+[w_move_id[j] for j in range(self.move_id_len)])
        
        # Coin flip    
        qc.append(coin_flip_gate,  [w_ancilla[j] for j in range(self.n_ancilla_bits)]+[w_coin[0],w_move_value[0]]+ [w_move_id[j] for j in range(self.move_id_len)]+[w_angles[k][j] for (k,j) in product(range(self.n_angles-1,-1,-1), range(self.angle_precision_bits))])

        # Conditional move
        qc.append(self.conditional_move_gate_n, [w_ancilla[j] for j in range(self.n_ancilla_bits)]+[w_coin[0],w_move_value[0]]+ [w_move_id[j] for j in range(self.move_id_len)]+[w_angles[k][j] for (k,j) in product(range(self.n_angles-1,-1,-1), range(self.angle_precision_bits))])

        # Inverse coin flip
        qc.append(coin_flip_gate.inverse(),[w_ancilla[j] for j in range(self.n_ancilla_bits)]+[w_coin[0],w_move_value[0],]+ [w_move_id[j] for j in range(self.move_id_len)]+[w_angles[k][j] for (k,j) in product(range(self.n_angles-1,-1,-1), range(self.angle_precision_bits))])

        # Inverse move preparation
        qc.append(self.move_preparation_gate.inverse(), [w_move_value[0]]+[w_move_id[j] for j in range(self.move_id_len)])

        # Reflection
        qc.append(self.reflection_gate, [w_coin[0],w_move_value[0]]+[w_move_id[j] for j in range(self.move_id_len)])

        W_gate = qc.to_instruction()
        
        return W_gate

    def U_func_n(self):
        
        '''
        This defines the gate U that initially spreads the output of minifold, and we can reuse its inverse too.
        It is basically the gate W but with the coin flip being sin^2 (theta = pi/6) = 1/4 probability of acceptance
        '''

        # State definition. All angles range from 0 to 2pi
        u_angles = []
        for i in range(self.n_angles):
            u_angles.append(QuantumRegister(self.angle_precision_bits, name = 'angle' + str(i)))

        # Move proposal
        u_move_id = QuantumRegister(self.move_id_len, name = 'move_id') #Which angle are we modifying
        u_move_value = QuantumRegister(1, name = 'move_value') #0 -> decrease the angle. 1-> increase it

        # Coin
        u_coin = QuantumRegister(1, name = 'coin')

        # Ancillas
        u_ancilla = QuantumRegister(self.probability_bits, name = 'ancilla')

        # Circuit
        qc = QuantumCircuit(u_ancilla, u_coin, u_move_value,u_move_id)
        for i in range(self.n_angles-1,-1,-1):
            qc = qc + QuantumCircuit(u_angles[i])


        # Move preparation
        qc.append(self.move_preparation_gate, [u_move_value[0]]+ [u_move_id[j] for j in range(self.move_id_len)])
        
        # Coin flip: equivalent to rx: https://qiskit.org/documentation/stubs/qiskit.circuit.library.U3Gate.html
        qc.u3( theta =  math.pi/3, phi = 0, lam = 0, qubit=u_coin)

        # Conditional move
        qc.append(self.conditional_move_gate_n, [u_ancilla[j] for j in range(self.n_ancilla_bits)]+[u_coin[0],u_move_value[0]]+ [u_move_id[j] for j in range(self.move_id_len)]+[u_angles[k][j] for (k,j) in product(range(self.n_angles-1,-1,-1), range(self.angle_precision_bits))])

        # Inverse coin flip
        qc.u3( theta = math.pi/3, phi = 0, lam = 0, qubit=u_coin).inverse()

        # Inverse move preparation
        qc.append(self.move_preparation_gate.inverse(), [u_move_value[0]]+ [u_move_id[j] for j in range(self.move_id_len)])

        # Reflection
        qc.append(self.reflection_gate, [u_coin[0],u_move_value[0]]+[u_move_id[j] for j in range(self.move_id_len)])

        U_gate = qc.to_instruction()
        
        return U_gate

    def execute_quantum_metropolis_n(self):

        # State definition. All angles range from 0 to 2pi
        # State definition. All angles range from 0 to 2pi

        g_angles = []
        for i in range(self.n_angles):
            g_angles.append(QuantumRegister(self.angle_precision_bits, name = 'angle' + str(i)))

        # Move proposal
        g_move_id = QuantumRegister(self.move_id_len, name = 'move_id') #Which angle are we modifying
        g_move_value = QuantumRegister(1, name = 'move_value') #0 -> decrease the angle. 1-> increase it

        # Coin
        g_coin = QuantumRegister(1, name = 'coin')

        # Ancillas
        g_ancilla = QuantumRegister(self.probability_bits, name = 'ancilla')

        # Circuit
        qc = QuantumCircuit(g_ancilla,g_coin,g_move_value,g_move_id)
        for i in range(self.n_angles-1,-1,-1):
            qc = qc + QuantumCircuit(g_angles[i])

        # Metropolis algorithm (create one input oracle for each beta)
        #list_gates = []

        # If initialization is totally mixed use
        for g_angle in g_angles:
            qc.h(g_angle)

        beta = self.beta_max
        #It creates one different oracle for each beta
        oracle = beta_precalc_TruthTableOracle.Beta_precalc_TruthTableOracle(self.input_oracle, beta, in_bits = self.n_angles*self.angle_precision_bits + self.move_id_len + 1,out_bits = self.probability_bits, scaling_factor = self.scaling_factor)
        
        W_gate = self.W_func_n(oracle)
        
        #list_gates.append(W_gate) # We deepcopy W_gate to not interfere with other calls

        for i in range(self.n_repetitions):
            #list_gates[i].params[0]= beta
            qc.append(W_gate,  [g_ancilla[j] for j in range(self.n_ancilla_bits)] + [g_coin[0],g_move_value[0]]+ [g_move_id[j] for j in range(self.move_id_len)] +[g_angles[k][j] for (k,j) in product(range(self.n_angles-1,-1,-1), range(self.angle_precision_bits))])

        start_time = time.time()
        
        backend = Aer.get_backend('statevector_simulator')
        backend_options = {"method" : "statevector"}
        experiment = execute(qc, backend, backend_options=backend_options)
        state_vector = Statevector(experiment.result().get_statevector(qc))

        # calculate the indices of the angles
        # the angles qubits are at the end of the statevector
        number_bits_angles = self.angle_precision_bits * self.n_angles
        total_bits = state_vector.num_qubits
        angle_qubits = [qubit_index for qubit_index in range ((total_bits - number_bits_angles), total_bits)]

        probabilities = state_vector.probabilities(angle_qubits)

        #state = qi.Statevector.from_instruction(qc)
        time_statevector = time.time() - start_time

        # Extract probabilities in the measurement of the angles phi and psi
        #probabilities = state.probabilities([j+self.n_ancilla_bits+2+self.move_id_len for j in range(self.angle_precision_bits * self.n_angles)])

        probs = {}
        for index_probabilites in range(2**(self.angle_precision_bits *self.n_angles)):

            key = self.convert_index_to_key(index_probabilites, self.angle_precision_bits, self.n_angles)
            probs[key] = probabilities[index_probabilites]#.as_integer

        return [probs, time_statevector]

# TO DO: initizialisation. 

    # this method converts the index returned by statevector into a string key. 
    # for example: key 10 is converted to 22 if there are two angles and two precision bits
    # for example: key 8 is converted to 0010 if there are four angles and three precision bits
    def convert_index_to_key(self, key_int, precision_bits, n_angles):

        key_str = ''

        # iterate over the number of angles
        for index_angle in range(n_angles):

            # generate a denominator to divide the key_int (integer key)
            # this denominator is equivalent to the 'weight' of this angle position
            # for example, if there are 4 angles, it goes to the first angle (from left) and calculate the denominator
            # then it goes to the next angle and calculate the new denominator
            denominator = 2**(precision_bits*((n_angles-1) - index_angle))

            result = int(key_int/denominator)
            key_str += str(result)

            # the key_int value is necessary to be updated
            key_int -= result * denominator

        return key_str