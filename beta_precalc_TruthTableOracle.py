from qiskit.aqua.components.oracles import Oracle, TruthTableOracle
from sympy.combinatorics.graycode import GrayCode, gray_to_bin, bin_to_gray
from qiskit.circuit import QuantumRegister, ClassicalRegister, QuantumCircuit, Qubit
import utils

import math

class Beta_precalc_TruthTableOracle(TruthTableOracle):
    '''Outputs the binary angle of rotation to get the correct probability. Tested ok'''
    def __init__(self, deltas_dictionary, beta, in_bits, out_bits, optimization=True, mct_mode='noancilla'):

        self.tools = utils.Utils()

        self.beta = beta
        self.in_bits = in_bits
        self.out_bits = out_bits
        self.deltas_dictionary = deltas_dictionary

        #self.circuit = self.create_circuit()
        
        self.bitmap = self.calculate_bitmap()
        for i in range(out_bits):
            print('In column',i,'there are', self.bitmap[i].count('1'), '1s and ', self.bitmap[i].count('0'), '0s.')
        print('\n')
        for bm in self.bitmap:
            print(len(bm))
        print('\n')

        super().__init__(self.bitmap, optimization, mct_mode)

    def calculate_bitmap(self):
        new_bitmap = []
        angles = {}

        for key in self.deltas_dictionary.keys():

            if self.deltas_dictionary[key] >= 0:
                probability = math.exp(-self.beta * self.deltas_dictionary[key])
            else: 
                probability = 1
                
            # Instead of encoding the probability, we will encode 1-probability. That way 1 -> 000, 
            # but if probability is 0 there is some small probability of acceptance
            probability = 1 - probability
            
            # Instead of probability save angles so rotations are easier to perform afterwards sqrt(p) = sin(theta)
            angle = math.asin(math.sqrt(probability))
            
            # Make the angle be between [0,1]. Since the maximum is pi/2
            angle /= (math.pi/2) 
            
            # Convert it into an integer and a string
            str_angle2 = self.int_angle_func(angle,self.out_bits)
            angles[key] = str_angle2

        # Printout
        for i in range(2**self.out_bits):
            st = self.tools.angle_to_binary(i, self.out_bits)
            print(st, 'appears', list(angles.values()).count(st), 'times in the angles dictionary')
        
        # Encoding the new bitmap
        new_bitmap = []
        for o in range(self.out_bits):
            string = ''
            for key in self.deltas_dictionary.keys():
                string += str(angles[key])[o]
            new_bitmap += [string]

        return new_bitmap

    #Method to convert angles int to binary
    def int_angle_func(self, angle,out_bits):
        out_str = ''
        a = angle
        for bits in range(1,out_bits+1):

            if a + 1e-10 > 1/(2**(bits)):
                out_str += '1'
                a -= 1/(2**bits)
            else:
                out_str += '0'
        
        return out_str

    '''
    def create_circuit(self):
        angles = {}

        for key in self.deltas_dictionary.keys():

            if self.deltas_dictionary[key] >= 0:
                probability = math.exp(-self.beta * self.deltas_dictionary[key])
            else: 
                probability = 1
                
            # Instead of encoding the probability, we will encode 1-probability. That way 1 -> 000, 
            # but if probability is 0 there is some small probability of acceptance
            probability = 1 - probability
            
            # Instead of probability save angles so rotations are easier to perform afterwards sqrt(p) = sin(theta)
            angle = math.asin(math.sqrt(probability))
            
            # Make the angle be between [0,1]. Since the maximum is pi/2
            angle /= (math.pi/2) 
            
            # Convert it into an integer and a string
            str_angle2 = self.int_angle_func(angle,self.out_bits)
            angles[key] = str_angle2

        # Printout
        for i in range(2**self.out_bits):
            st = self.tools.angle_to_binary(i, self.out_bits)
            print(st, 'appears', list(angles.values()).count(st), 'times in the angles dictionary')

        # Now prepare the circuit
        in_reg = QuantumRegister(self.in_bits)
        out_reg = QuantumRegister(self.out_bits)
        circuit = QuantumCircuit(in_reg,out_reg)

        circuit.x(in_reg)
        prev_state = '0'*self.in_bits
        indices = GrayCode(self.in_bits)

        for index in indices:
            out = angles[index]
            if out != '0'*self.out_bits:

                
                for i in range(self.out_bits):
                    if out[i] == '1':
                        circuit.mcx(control_qubits= in_reg[:], target_qubit= out_reg[i])

                # Save which was the previous controlled state        
                prev_state = index
    '''