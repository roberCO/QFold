# Importing standard Qiskit libraries and configuring account
from qiskit import QuantumCircuit, execute, Aer, IBMQ
from qiskit.compiler import transpile, assemble
import qiskit.quantum_info as qi

# Import libraries
import qiskit
from qiskit.circuit import QuantumRegister, ClassicalRegister, QuantumCircuit, Qubit, Clbit, Gate, Parameter, InstructionSet
from qiskit.aqua.components.oracles import Oracle, TruthTableOracle
from qiskit import BasicAer
import numpy as np
import math
from copy import deepcopy

import logging
import beta_precalc_TruthTableOracle

from math import pi

from qiskit.aqua.utils.controlled_circuit import apply_cu3
from qiskit.aqua import AquaError
from itertools import product
import time


class QuantumMetropolis():

    def __init__(self, n_repetitions, angle_precision_bits, probability_bits, n_angles, beta_max, input_oracle):

        #Global variables

        # Number steps
        self.n_repetitions = n_repetitions

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

    def move_preparation(self, circuit,move_id,move_value):
        '''
        Proposes new moves
        '''
        circuit.h(move_value) #Es un único bit que puede ser 0 para + y 1 para -
        if bin(self.n_angles).count('1') == 1: # if the number of angles is a power of two
            circuit.h(move_id)
        else:
            vector = ([1]*self.n_angles + [0]*(2**(self.move_id_len) - self.n_angles))/np.sqrt(self.n_angles)
            circuit.initialize(vector, [move_id[j] for j in range(self.move_id_len)])

    def conditional_move_npeptide(self,circuit,angles,move_id,move_value,coin,ancilla):
        '''
        Conditioned on coin, perform a move. Tested ok.
        We use a repetitive structure where we perform the conditional sum and substraction for each angle.
        '''
        # For each angle
        for i in range(self.n_angles):
            angle = angles[i] #Select the angle from the list of registers
            angle_index = np.binary_repr(i, width=self.move_id_len) #convert i to binary

            # Put the given move_id in all 1 to control on it: for instance if we are controling on i=2, move 010 ->111
            for j in range(len(angle_index)):
                if angle_index[j] == '0':
                    circuit.x(move_id[j])

            circuit.mcx(control_qubits= [move_id[j] for j in range(move_id.size)] +[coin[0]], target_qubit = ancilla[0])#create a single control
            self.sumsubstract1(circuit,angle,ancilla[0],ancilla[1],ancilla[2],move_value[0]) #sum or substract 1 to the angle
            circuit.mcx(control_qubits= [move_id[j] for j in range(move_id.size)] +[coin[0]], target_qubit = ancilla[0])#create a single control        
            
            # Undo the move_id preparation: for instance, if we are controlling on i= 2 move 111->010
            for j in range(len(angle_index)):
                if angle_index[j] == '0':
                    circuit.x(move_id[j])

    def reflection(self, circuit,move_id,move_value,coin):
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

        sub_circ = QuantumCircuit(s_move_id, s_move_value)
        self.move_preparation(sub_circ,s_move_id,s_move_value)
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
        sub_circ = QuantumCircuit(s_angles[0])
        for i in range(1, self.n_angles):
            sub_circ = sub_circ + QuantumCircuit(s_angles[i])
        sub_circ = sub_circ + QuantumCircuit(s_move_id, s_move_value, s_coin, s_ancilla)

        self.conditional_move_npeptide(sub_circ,s_angles,s_move_id,s_move_value,s_coin,s_ancilla)

        # Optimize the circuit

        print('Before optimization------- conditional_move_npeptide')
        print('gates = ', sub_circ.count_ops())
        print('depth = ', sub_circ.depth())
        sub_circ = transpile(sub_circ, seed_transpiler=1, optimization_level=3)
        print('After optimization--------')
        print('gates = ', sub_circ.count_ops())
        print('depth = ', sub_circ.depth())

        conditional_move_gate_n = sub_circ.to_instruction()

        # Reflection gate --------------------------------------------------
        s_move_id = QuantumRegister(self.move_id_len) 
        s_move_value = QuantumRegister(1)
        s_coin = QuantumRegister(1)

        sub_circ = QuantumCircuit(s_move_id,s_move_value,s_coin)
        self.reflection(sub_circ,s_move_id,s_move_value,s_coin)

        # We could optimize the circuit, but it will probably not be worth it (fairly small)

        reflection_gate = sub_circ.to_instruction()

        return [move_preparation_gate, conditional_move_gate_n, reflection_gate]

    def sum1(self, circuit,qubit_string,control,start,end):
        
        n_qubits = qubit_string.size     # calculate n_qubits
        circuit.cx(control,end) # iff control = 1, end = 1
        circuit.x(start)
        circuit.cx(control,start) # iff control = 1, start = 0
        
        for i in range(n_qubits+1): #Don't need to add control, since start already does that work
            '''
            Next thing we analise if all qubits to the right have value 1, 
            and save it in the current qubit and start
            '''
            if i > 0:
                # For i = 0, there is only the start to worry about
                circuit.mcx(control_qubits= [qubit_string[j] for j in range(n_qubits-i)]+[end], target_qubit = qubit_string[n_qubits-i])
            circuit.mcx(control_qubits = [qubit_string[j] for j in range(n_qubits-i)]+[end], target_qubit = start)

            '''
            Next, controlling on the current qubit and start, we change all the following qubits to 0.
            We have to control with qubit_string[n_qubit]
            '''
            if i == 0:
                for j in range(n_qubits-i):
                    circuit.ccx(control,start,qubit_string[j])
                circuit.ccx(control,start,end)
            elif i == n_qubits:
                circuit.mcx(control_qubits = [control,qubit_string[n_qubits-i],start], target_qubit = end)
            else:
                for j in range(n_qubits-i):            
                    circuit.mcx(control_qubits = [control,qubit_string[n_qubits-i],start], target_qubit = qubit_string[j])
                circuit.mcx(control_qubits = [control,qubit_string[n_qubits-i],start], target_qubit = end)
        circuit.x(start)

    def substract1(self, circuit,qubit_string,control,start,end):
        '''
        Outputs:
        Substracts register 2 (1 qubit) from register 1 in register 1. Tested ok.
        
        Input:
        circuit: QuantumCircuit with registers qubit_string, control, ancilla
        
        qubit_string: QuantumRegister
        
        control: Qubit. Use ancilla[0] or similar
        
        start: Qubit. Use ancilla[1] or similar
        end: Qubit. Use ancilla[2] or similar
        
        Comments: In binary, substracting is the same procedure as summing when we exchange 0s and 1s
        '''
        circuit.x(qubit_string)

        self.sum1(circuit,qubit_string,control,start,end)
        
        circuit.x(qubit_string)

    def sumsubstract1(self,circuit,qubit_string,control,start,end,move_value):
        '''
        Outputs:
        Sum/Substracts register 2 (control, 1 qubit) from register 1 (qubit_string) in register 1. Tested ok.

        Input:
        circuit: QuantumCircuit with registers qubit_string, control, ancilla and move_value
        qubit_string: QuantumRegister where the sum/substraction is performed
        control: Qubit. Use ancilla[0] or similar. It encodes the probability of change.
        start: Qubit. Use ancilla[1] or similar
        end: Qubit. Use ancilla[2] or similar
        move_value: 1 to substract, 0 to sum

        Comments: In binary, substracting is the same procedure as summing when we exchange 0s and 1s
        '''
        circuit.cx(move_value,qubit_string)

        self.sum1(circuit,qubit_string,control,start,end)

        circuit.cx(move_value,qubit_string)

    def coin_flip(self, circuit,coin,ancilla):
        '''
        Prepares the coin with the probability encoded in the ancilla.
        coin: register[0]
        ancilla: register[1:]
        '''
        
        #Necesitamos usar el número guardado en las ancillas para realizar rotaciones controladas.  
        #Notice that ancilla encodes 1-probability, rather than probability.
        #Notice also that cu3(theta) rotates theta/2. As the first angle to rotate is pi/4 we need to start in theta = pi/2

        circuit.x(coin) # Start in 1 and decrease it, since we encoded the angle corresponding 1-probability
        for i in range(ancilla.size): # See how to perform an rx rotation in https://qiskit.org/documentation/stubs/qiskit.circuit.library.U3Gate.html
            circuit.cu3(theta = -math.pi/(2**(i+1)), phi  = -math.pi/2, lam = math.pi/2, control_qubit = ancilla[i], target_qubit = coin)
        
    def coin_flip_func_n(self, oracle):
        
        '''
        Defines de coin_flip_gate using the oracle that is provided on the moment.
        Notice that oracle gate has registers oracle.variable_register and oracle.output_register in that order
        oracle.variable_register should have size angle_phi.size + angle_psi.size + move_id.size + move_value.size
        oracle.output_register should have size ancilla.size
        '''

        # Construct an instruction for the oracle
        oracle.construct_circuit()
        oracle_circuit = oracle.circuit

        # Optimize the circuit of the oracle. First try it with printout
        print('<i> Oracle gate and depth count')
        print('gates = ', oracle_circuit.count_ops())
        print('depth = ', oracle_circuit.depth())

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

        cf_circ = QuantumCircuit(cf_angles[0])
        for i in range(1, self.n_angles):
            cf_circ = cf_circ + QuantumCircuit(cf_angles[i])
        cf_circ = cf_circ + QuantumCircuit(cf_move_id, cf_move_value ,cf_coin , cf_ancilla)

        # Main operations
        cf_circ.append(oracle_gate, [cf_angles[k][j] for (k,j) in product(range(self.n_angles), range(self.angle_precision_bits))]+ [cf_move_id[j] for j in range(cf_move_id.size)] + [cf_move_value[0]] +[cf_ancilla[j] for j in range(cf_ancilla.size)])
        self.coin_flip(cf_circ,cf_coin,cf_ancilla)
        cf_circ.append(oracle_gate.inverse(), [cf_angles[k][j] for (k,j) in product(range(self.n_angles), range(self.angle_precision_bits))]+ [cf_move_id[j] for j in range(cf_move_id.size)] + [cf_move_value[0]] +[cf_ancilla[j] for j in range(cf_ancilla.size)])
        
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
        qc = QuantumCircuit(w_angles[0])
        for i in range(1, self.n_angles):
            qc = qc + QuantumCircuit(w_angles[i])
        qc = qc + QuantumCircuit(w_move_id, w_move_value ,w_coin , w_ancilla)
        
        # Define the coin_flip_gate
        self.coin_flip_gate = self.coin_flip_func_n(oracle)

        # Move preparation
        qc.append(self.move_preparation_gate, [w_move_id[j] for j in range(self.move_id_len)]+[w_move_value[0]])
        
        # Coin flip    
        qc.append(self.coin_flip_gate, [w_angles[k][j] for (k,j) in product(range(self.n_angles), range(self.angle_precision_bits))] + [w_move_id[j] for j in range(self.move_id_len)] +[ w_move_value[0],w_coin[0]] + [w_ancilla[j] for j in range(w_ancilla.size)])

        # Conditional move
        qc.append(self.conditional_move_gate_n, [w_angles[k][j] for (k,j) in product(range(self.n_angles), range(self.angle_precision_bits))] + [w_move_id[j] for j in range(self.move_id_len)] +[ w_move_value[0],w_coin[0]] + [w_ancilla[j] for j in range(w_ancilla.size)])

        # Inverse coin flip
        qc.append(self.coin_flip_gate.inverse(), [w_angles[k][j] for (k,j) in product(range(self.n_angles), range(self.angle_precision_bits))] + [w_move_id[j] for j in range(self.move_id_len)] +[ w_move_value[0],w_coin[0]] + [w_ancilla[j] for j in range(w_ancilla.size)])

        # Inverse move preparation
        qc.append(self.move_preparation_gate.inverse(), [w_move_id[j] for j in range(self.move_id_len)]+[w_move_value[0]])

        # Reflection
        qc.append(self.reflection_gate, [w_move_id[j] for j in range(self.move_id_len)]+[w_move_value[0],w_coin[0]])

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
        qc = QuantumCircuit(u_angles[0])
        for i in range(1, self.n_angles):
            qc = qc + QuantumCircuit(u_angles[i])
        qc = qc + QuantumCircuit(u_move_id, u_move_value , u_coin , u_ancilla)

        # Move preparation
        qc.append(self.move_preparation_gate, [u_move_id[j] for j in range(self.move_id_len)]+[u_move_value[0]])
        
        # Coin flip: equivalent to rx: https://qiskit.org/documentation/stubs/qiskit.circuit.library.U3Gate.html
        qc.u3( theta =  math.pi/6, phi = -math.pi/2, lam = math.pi/2, qubit=u_coin)

        # Conditional move
        qc.append(self.conditional_move_gate_n, [u_angles[k][j] for (k,j) in product(range(self.n_angles), range(self.angle_precision_bits))] + [u_move_id[j] for j in range(self.move_id_len)] +[ u_move_value[0],u_coin[0]] + [u_ancilla[j] for j in range(self.probability_bits)])

        # Inverse coin flip
        qc.u3( theta = -math.pi/6, phi = -math.pi/2, lam = math.pi/2, qubit=u_coin)
        # Inverse move preparation
        qc.append(self.move_preparation_gate.inverse(), [u_move_id[j] for j in range(self.move_id_len)]+[u_move_value[0]])

        # Reflection
        qc.append(self.reflection_gate, [u_move_id[j] for j in range(self.move_id_len)]+[u_move_value[0],u_coin[0]])

        U_gate = qc.to_instruction()
        
        return U_gate

    def execute_quantum_metropolis_n(self):

        # State definition. All angles range from 0 to 2pi
        # State definition. All angles range from 0 to 2pi

        print('<i> Execute n quantum metropolis')
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
        qc = QuantumCircuit(g_angles[0])
        for i in range(1, self.n_angles):
            qc = qc + QuantumCircuit(g_angles[i])
        qc = qc + QuantumCircuit(g_move_id, g_move_value ,g_coin , g_ancilla)

        # Metropolis algorithm (create one input oracle for each beta)
        list_gates = []

        # If initialization is totally mixed use
        # for g_angle in g_angles:
        #   qc.h(g_angle)

        # If initialization is from minifold

        # qc.x (SELECT THE BITS THAT HAVE TO BE PUT AT STATE 1) TO BE DONE!!!

        print('<i> Calculating gates')

        U_gate = self.U_func_n()

        for i in range(3):
            qc.append(U_gate, [g_angles[k][j] for (k,j) in product(range(self.n_angles), range(self.angle_precision_bits))] + [g_move_id[j] for j in range(self.move_id_len)] + [g_move_value[0],g_coin[0]] + [g_ancilla[j] for j in range(g_ancilla.size)])

        # End of initialization

        for i in range(self.n_repetitions):

            print('<i> Quantum metropolis ->', i,'repetition')
            
            beta = ((1+i)/self.n_repetitions)*self.beta_max
            
            #It creates one different oracle for each beta
            oracle = beta_precalc_TruthTableOracle.Beta_precalc_TruthTableOracle(self.input_oracle, beta, in_bits = self.n_angles*self.angle_precision_bits + self.move_id_len + 1,out_bits = self.probability_bits)
            print('<i> Oracle created')
            
            W_gate = self.W_func_n(oracle)
            print('<i> w gate created')
            
            list_gates.append(W_gate) # We deepcopy W_gate to not interfere with other calls
            print('<i> w_gate added')

            #list_gates[i].params[0]= beta
            qc.append(W_gate, [g_angles[k][j] for (k,j) in product(range(self.n_angles), range(self.angle_precision_bits))] + [g_move_id[j] for j in range(self.move_id_len)] + [g_move_value[0],g_coin[0]] + [g_ancilla[j] for j in range(g_ancilla.size)])
            print('<i> q circuit created')

            print('\n')

        start_time = time.time()
        
        print('<i> Calculating statevector')
        state = qi.Statevector.from_instruction(qc)
        print("<i>QUANTUM METROPOLIS: Time to calculate statevector: %s seconds" % (time.time() - start_time))

        # Extract probabilities in the measurement of the angles phi and psi
        probabilities = state.probabilities([j for j in range(self.angle_precision_bits * self.n_angles)])

        probs = {}
        for i in range(2**(self.angle_precision_bits *self.n_angles)):

            b = str(np.binary_repr(i, width = self.angle_precision_bits *self.n_angles))
            probs[b] = probabilities[i]

        print('<i> Final probabilities', probs)
        return probs

# TO DO: initizialisation. 