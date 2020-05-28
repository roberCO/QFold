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

from sympy.combinatorics.graycode import GrayCode
from math import pi

from qiskit.aqua.utils.controlled_circuit import apply_cu3
from qiskit.aqua import AquaError
from itertools import product
import time
import progressbar


class QuantumMetropolis():

    def __init__(self, n_repetitions, n_precision_bits, n_ancilla_bits, beta_max, input_oracle):

        #Global variables

        # Number steps
        self.n_repetitions = n_repetitions
        # Dipeptide: Indicate number of precision bits
        self.n_precision_bits = n_precision_bits
        #Oracle ancilla bits
        self.n_ancilla_bits = n_ancilla_bits
        self.beta_max = beta_max

        self.input_oracle = input_oracle

        [self.move_preparation_gate, self.conditional_move_gate, self.reflection_gate] = self.prepare_initial_circuits()

    # This is the move preparation gate
    def move_preparation(self, circuit,move_id,move_value):
        '''
        Proposes new moves
        '''
        circuit.h(move_id)
        circuit.h(move_value)

    # Use as qc.append(move_preparation, [move_id[0], move_value[0]])
    def conditional_move_dipeptide(self, circuit,coin,move_id,move_value,angle_phi,angle_psi,ancilla):
        '''
        Conditioned on coin, perform a move. For a dipeptide!
        We use a repetitive structure where we perform the conditional sum and substraction for each angle.
        Checked ok
        '''
        # For angle_phi = angle_id = 0 ----------------------------------------------
        circuit.x(move_id) # Put move_id in 1
        
        # Conditional on move_id = 0, move_value = 1 and coin = 1, increase angle_phi by on
        circuit.mcrx(theta = pi, q_controls = [move_id[0], move_value[0], coin[0]], q_target = ancilla[0])
        #circuit.append(Gate("mcx", 5, []), [move_id[0], move_id[1], move_value, coin, ancilla[0]]) #create a single control
        self.sum1(circuit,angle_phi,ancilla[0],ancilla[1],ancilla[2]) #calc_anc_size(len_angle_value) =6 in this cased
        circuit.mcrx(theta = pi, q_controls = [move_id[0], move_value[0], coin[0]], q_target = ancilla[0])
        
        # Conditional on move_id = 0, move_value = 0 and coin = 1, decrease angle_ psi by on
        # Put move_id in 11
        circuit.x(move_value)
        circuit.mcrx(theta = pi, q_controls = [move_id[0], move_value[0], coin[0]], q_target = ancilla[0])#create a single control
        self.substract1(circuit,angle_phi,ancilla[0],ancilla[1],ancilla[2]) #calc_anc_size(len_angle_value) =6 in this cased
        circuit.mcrx(theta = pi, q_controls = [move_id[0], move_value[0], coin[0]], q_target = ancilla[0])#create a single control
        circuit.x(move_value)
        
        circuit.x(move_id)
        
        # For angle_psi = angle_id = 1 ----------------------------------------------
        # Put move_id in 11
        
        # Conditional on move_id = 1, move_value = 1 and coin = 1, increase angle_psi by one
        circuit.mcrx(theta = pi, q_controls = [move_id[0],move_value[0], coin[0]], q_target = ancilla[0])#create a single control
        self.sum1(circuit,angle_psi,ancilla[0],ancilla[1],ancilla[2]) #calc_anc_size(len_angle_value) =6 in this cased
        circuit.mcrx(theta = pi, q_controls = [move_id[0],move_value[0], coin[0]], q_target = ancilla[0])#create a single control
        
        # Conditional on move_id = 1, move_value = 0 and coin = 1, decrease angle_psi by one
        # Put move_id in 11
        circuit.x(move_value)
        circuit.mcrx(theta = pi, q_controls = [move_id[0], move_value[0], coin[0]], q_target = ancilla[0])#create a single control
        self.substract1(circuit,angle_psi,ancilla[0],ancilla[1],ancilla[2]) #calc_anc_size(len_angle_value) =6 in this cased
        circuit.mcrx(theta = pi, q_controls = [move_id[0], move_value[0], coin[0]], q_target = ancilla[0])#create a single control
        circuit.x(move_value)

    # To use it or its inverse
    # qc.append(conditional_move.inverse(), [(angle_phi[j] for j in range(angle_phi.size)),(angle_psi[j] for j in range(angle_psi.size)),move_id[0], move_value[0], coin[0], ancilla[0], ancilla[1], ancilla[2]])
    def reflection(self, circuit,move_id,move_value,coin):
        '''
        I have to investigate over what is the reflection performed. Is it performed over 000?
        If in state 0000, make it 1111, cccz gate and back to 0000
        '''
        circuit.x(move_id)
        circuit.x(move_value)
        circuit.x(coin)
        
        circuit.mcrz(lam = pi, q_controls = [move_id[0]] + [move_value[0]], q_target = coin[0]) #For dipeptide
        # circuit.mcrz(lam = pi, q_controls = [move_id[0]] + [move_id[1]]+ [move_value[0]], q_target = coin[0]) #For tripeptide
        
        circuit.x(move_id)
        circuit.x(move_value)
        circuit.x(coin)

    def prepare_initial_circuits(self):

        # Code the moves as a gate
        s_move_id = QuantumRegister(1) 
        s_move_value = QuantumRegister(1)

        sub_circ = QuantumCircuit(s_move_id, s_move_value)
        self.move_preparation(sub_circ,s_move_id,s_move_value)
        move_preparation_gate = sub_circ.to_instruction()


        # Let us define it as a portable gate
        s_move_id = QuantumRegister(1) 
        s_move_value = QuantumRegister(1)
        s_coin = QuantumRegister(1)
        s_ancilla = QuantumRegister(self.n_ancilla_bits)
        s_angle_phi = QuantumRegister(self.n_precision_bits, name = 'angle_phi')
        s_angle_psi = QuantumRegister(self.n_precision_bits, name = 'angle_psi') 

        sub_circ = QuantumCircuit(s_angle_phi,s_angle_psi,s_move_id,s_move_value,s_coin,s_ancilla)

        self.conditional_move_dipeptide(sub_circ,s_coin,s_move_id,s_move_value,s_angle_phi,s_angle_psi,s_ancilla)
        conditional_move_gate = sub_circ.to_instruction()


        # Let us define it as a portable gate
        s_move_id = QuantumRegister(1) 
        s_move_value = QuantumRegister(1)
        s_coin = QuantumRegister(1)

        sub_circ = QuantumCircuit(s_move_id,s_move_value,s_coin)
        self.reflection(sub_circ,s_move_id,s_move_value,s_coin)
        reflection_gate = sub_circ.to_instruction()

        return [move_preparation_gate, conditional_move_gate, reflection_gate]


    # Add one to the circuit passed as parameter
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
                circuit.mcrx(theta = pi, q_controls = [qubit_string[j] for j in range(n_qubits-i)]+[end], q_target = qubit_string[n_qubits-i])
            circuit.mcrx(theta = pi, q_controls = [qubit_string[j] for j in range(n_qubits-i)]+[end], q_target = start)

            '''
            Next, controlling on the current qubit and start, we change all the following qubits to 0.
            We have to control with qubit_string[n_qubit]
            '''
            if i == 0:
                for j in range(n_qubits-i):
                    circuit.ccx(control,start,qubit_string[j])
                circuit.ccx(control,start,end)
            elif i == n_qubits:
                circuit.mcrx(theta = pi, q_controls = [control,qubit_string[n_qubits-i],start], q_target = end)
            else:
                for j in range(n_qubits-i):            
                    circuit.mcrx(theta = pi, q_controls = [control,qubit_string[n_qubits-i],start], q_target = qubit_string[j])
                circuit.mcrx(theta = pi, q_controls = [control,qubit_string[n_qubits-i],start], q_target = end)
        circuit.x(start)

    #Substract one to the circuit passed as parameter
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
        for i in range(ancilla.size):
            circuit.cu3(theta = -math.pi/(2**(i+1)), phi  = 0, lam = 0, control_qubit = ancilla[i], target_qubit = coin)
        
    def coin_flip_func(self, oracle):
        
        '''
        Defines de coin_flip_gate using the oracle that is provided on the moment.
        Notice that oracle gate has registers oracle.variable_register and oracle.output_register in that order
        oracle.variable_register should have size angle_phi.size + angle_psi.size + move_id.size + move_value.size
        oracle.output_register should have size ancilla.size
        '''

        # Construct an instruction for the oracle
        oracle.construct_circuit()
        oracle_circuit = oracle.circuit
        #print(oracle_circuit)
        oracle_gate = oracle_circuit.to_instruction()


        # Let us create a circuit for coin_flip
        cf_move_id = QuantumRegister(1) 
        cf_move_value = QuantumRegister(1)
        cf_coin = QuantumRegister(1)
        cf_ancilla = QuantumRegister(self.n_ancilla_bits)
        cf_angle_phi = QuantumRegister(self.n_precision_bits, name = 'angle_phi')
        cf_angle_psi = QuantumRegister(self.n_precision_bits, name = 'angle_psi') 

        cf_circ = QuantumCircuit(cf_angle_phi, cf_angle_psi, cf_move_id, cf_move_value ,cf_coin , cf_ancilla)

        # Main operations
        cf_circ.append(oracle_gate, [cf_angle_phi[j] for j in range(cf_angle_phi.size)]+[cf_angle_psi[j] for j in range(cf_angle_psi.size)]+ [cf_move_id[0], cf_move_value[0]] +[cf_ancilla[j] for j in range(cf_ancilla.size)])
        self.coin_flip(cf_circ,cf_coin,cf_ancilla)
        cf_circ.append(oracle_gate.inverse(), [cf_angle_phi[j] for j in range(cf_angle_phi.size)]+[cf_angle_psi[j] for j in range(cf_angle_psi.size)]+ [cf_move_id[0], cf_move_value[0]] +[cf_ancilla[j] for j in range(cf_ancilla.size)])

        coin_flip_gate = cf_circ.to_instruction()
        
        return coin_flip_gate

    '''
    Use as:
        #coin_flip_gate.params[0]= a_given_beta
        qc.append(coin_flip_gate.inverse(), [angle_phi[j] for j in range(angle_phi.size)]+[angle_psi[j] for j in range(angle_psi.size)] + [move_id[0], move_value[0],coin[0]] + [ancilla[j] for j in range(ancilla.size)])
    '''
    def W_func(self, oracle):
        
        '''This defines the parametrised gate W using the oracle that is provided to it, and we can reuse its inverse too.'''

        # State definition. All angles range from 0 to 2pi
        w_angle_phi = QuantumRegister(self.n_precision_bits, name = 'angle_phi')
        w_angle_psi = QuantumRegister(self.n_precision_bits, name = 'angle_psi') 

        # Move proposal
        w_move_id = QuantumRegister(1, name = 'move_id') #Which angle are we modifying
        w_move_value = QuantumRegister(1, name = 'move_value') #0 -> decrease the angle. 1-> increase it

        # Coin
        w_coin = QuantumRegister(1, name = 'coin')

        # Ancillas
        w_ancilla = QuantumRegister(self.n_ancilla_bits, name = 'ancilla')

        # Circuit
        qc = QuantumCircuit(w_angle_phi,w_angle_psi,w_move_id,w_move_value,w_coin,w_ancilla)

        #D beta = Parameter('β') ---- Deprecated
        #D coin_flip_gate.params[0]= beta ---- Deprecated
        
        # Define the coin_flip_gate
        coin_flip_gate = self.coin_flip_func(oracle)

        # Move preparation
        qc.append(self.move_preparation_gate, [w_move_id[0], w_move_value[0]])
        
        # Coin flip    
        qc.append(coin_flip_gate, [w_angle_phi[j] for j in range(w_angle_phi.size)]+[w_angle_psi[j] for j in range(w_angle_psi.size)] + [w_move_id[0], w_move_value[0],w_coin[0]] + [w_ancilla[j] for j in range(w_ancilla.size)])

        # Conditional move
        qc.append(self.conditional_move_gate, [w_angle_phi[j] for j in range(w_angle_phi.size)]+[w_angle_psi[j] for j in range(w_angle_psi.size)] + [w_move_id[0], w_move_value[0], w_coin[0]] + [w_ancilla[j] for j in range(w_ancilla.size)])

        # Inverse coin flip
        qc.append(coin_flip_gate.inverse(), [w_angle_phi[j] for j in range(w_angle_phi.size)]+[w_angle_psi[j] for j in range(w_angle_psi.size)] + [w_move_id[0], w_move_value[0],w_coin[0]] + [w_ancilla[j] for j in range(w_ancilla.size)])

        # Inverse move preparation
        qc.append(self.move_preparation_gate, [w_move_id[0], w_move_value[0]])

        # Reflection
        qc.append(self.reflection_gate, [w_move_id[0], w_move_value[0],w_coin[0]])

        W_gate = qc.to_instruction()
        
        return W_gate

    '''
    Use as:
        #W_gate.params[0]= a_given_beta
        qc.append(W_gate.inverse(), [(angle_psi[j] for j in range(angle_psi.size)),angle_phi[j] for j in range(angle_phi.size)),move_id[0], move_value[0],coin[0],(ancilla[j] for j in range(ancilla.size))])
    '''

    def execute_quantum_metropolis(self):

        # State definition. All angles range from 0 to 2pi
        g_angle_phi = QuantumRegister(self.n_precision_bits, name = 'angle_phi')
        g_angle_psi = QuantumRegister(self.n_precision_bits, name = 'angle_psi') 

        # Move proposal
        g_move_id = QuantumRegister(1, name = 'move_id') #Which angle are we modifying
        g_move_value = QuantumRegister(1, name = 'move_value') #0 -> decrease the angle. 1-> increase it

        # Coin
        g_coin = QuantumRegister(1, name = 'coin')

        # Ancillas
        g_ancilla = QuantumRegister(self.n_ancilla_bits, name = 'ancilla')

        # Circuit
        qc = QuantumCircuit(g_angle_phi,g_angle_psi,g_move_id,g_move_value,g_coin,g_ancilla)

        # Metropolis algorithm (create one input oracle for each beta)
        print('    ⬤  Generating one oracle per each beta')
        list_gates = []
        bar = progressbar.ProgressBar(maxval=self.n_repetitions, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        for i in range(self.n_repetitions):
            
            beta = (1+i)/self.n_repetitions*self.beta_max
            
            #It creates one different oracle for each beta
            oracle = beta_precalc_TruthTableOracle.Beta_precalc_TruthTableOracle(self.input_oracle, beta, out_bits = self.n_ancilla_bits)
            
            W_gate = self.W_func(oracle)
            
            list_gates.append(W_gate) # We deepcopy W_gate to not interfere with other calls
            #list_gates[i].params[0]= beta
            qc.append(W_gate, [g_angle_phi[j] for j in range(g_angle_phi.size)] + [g_angle_psi[j] for j in range(g_angle_psi.size)] + [g_move_id[0], g_move_value[0],g_coin[0]] + [g_ancilla[j] for j in range(g_ancilla.size)])
            bar.update(i+1)

        bar.finish()

        # If instead we want to return the statevector
        print('    ⬤  Calculating statevector')
        start_time = time.time()
        state = qi.Statevector.from_instruction(qc)
        print("--- %s seconds ---" % (time.time() - start_time))

        print('    ⬤  Calculating probabilities')
        # Extract probabilities in the measurement of the angles phi and psi
        probabilities = state.probabilities([j for j in range(self.n_precision_bits * 2)])

        print('    ⬤  Calculating relevant probabilities')
        relevant_probabilities = []
        probs = []
        for i in range(2**(self.n_precision_bits *2)):

            probs.append(probabilities[i])
            if (i+1) % 2**self.n_precision_bits == 0:
                relevant_probabilities.append(probs)
                probs = []

        print('    ⬤  Returning quantum probabilities')
        return relevant_probabilities