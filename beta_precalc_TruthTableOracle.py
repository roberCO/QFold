from qiskit.aqua.components.oracles import Oracle, TruthTableOracle
from sympy.combinatorics.graycode import GrayCode, gray_to_bin, bin_to_gray
from qiskit.circuit import QuantumRegister, ClassicalRegister, QuantumCircuit, Qubit
import utils
from collections import OrderedDict

import math
import numpy as np

class Beta_precalc_TruthTableOracle(TruthTableOracle):
    '''Outputs the binary angle of rotation to get the correct probability. Tested ok'''
    def __init__(self, deltas_dictionary, beta, in_bits, out_bits, optimization=True, mct_mode='noancilla'):

        self.tools = utils.Utils()

        self.beta = beta
        self.in_bits = in_bits
        self.out_bits = out_bits
        self.deltas_dictionary = OrderedDict(sorted(deltas_dictionary.items()))
        
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
                
            # Instead of encoding the angle corresponding to the probability, we will encode the angle theta such that sin^2(pi/2 - theta) = probability.
            # That way 1 -> 000, but if probability is 0 there is some small probability of acceptance
            
            # Instead of probability save angles so rotations are easier to perform afterwards sqrt(p) = sin(pi/2-theta/2).
            # The theta/2 is because if you input theta, qiskits rotates theta/2. Also normalised (divided between pi the result)
            angle = 1 - 2/math.pi * math.asin(math.sqrt(probability))
            
            # Convert it into an integer and a string
            if angle == 1.:
                raise ValueError('Warning: angle seems to be pi/2, and that should not be possible')
            
            # angle will be between 0 and 1, so we move it to between 0 and 2^out_bits. Then calculate the integer and the binary representation
            angles[key] = np.binary_repr(int(angle*2**self.out_bits), width= self.out_bits)

        # Order angles by key
        angles = OrderedDict(sorted(angles.items()))

        # Printout
        for i in range(2**self.out_bits):
            st = np.binary_repr(i, width = self.out_bits)
            print(st, 'appears', list(angles.values()).count(st), 'times in the angles dictionary')
        
        # Encoding the new bitmap
        new_bitmap = []
        for o in range(self.out_bits):
            string = ''
            for key in angles.keys():
                string += str(angles[key])[o]
            new_bitmap += [string]

        return new_bitmap