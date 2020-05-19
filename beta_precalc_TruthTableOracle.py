from qiskit.aqua.components.oracles import Oracle, TruthTableOracle

import math

class Beta_precalc_TruthTableOracle(TruthTableOracle):
    '''Outputs the binary angle of rotation to get the correct probability. Tested ok'''
    def __init__(self, energies_dictionary, beta, out_bits, optimization=True, mct_mode='noancilla'):
        self.beta = beta
        self.out_bits = out_bits
        self.in_bits = len(list(energies_dictionary.keys())[0]) # The key of energies_dictionary is the input bits
        self.energies_dictionary = energies_dictionary
        #print(self.energies_dictionary)
        # print('Is parameter',self.parameters) --- Deprecated
        self.calculate_bitmap()
        super().__init__(self.bitmap, optimization, mct_mode)
        
    def calculate_bitmap(self):
        new_bitmap = []
        angles = {}
        for i in range(int(2**self.in_bits)):
            st = '0'*(self.in_bits - len(str(format(i,'b')))) + str(format(i,'b'))
            #print(st)
            if self.energies_dictionary[st] >= 0:
                #print(type(np))
                probability = math.exp(-self.beta * self.energies_dictionary[st])
                #print('probability: ',probability)
                #print('energy: ', self.energies_dictionary[st])
                #print('beta:', self.beta)
            else: 
                probability = 1
                #print('probability = 1')
                
            # Instead of encoding the probability, we will encode 1-probability. That way 1 -> 000, 
            #but if probability is 0 there is some small probability of acceptance
            probability = 1 - probability
            
            
            # Instead of probability save angles so rotations are easier to perform afterwards sqrt(p) = sin(theta)
            angle = math.asin(math.sqrt(probability))
            
            # Make the angle be between [0,1]. Since the maximum is pi/2
            angle /= (math.pi/2) 
            
            # Convert it into an integer and a string
            #int_angle = format(int(angle*2**self.out_bits), 'b')
            #str_angle = str(int_angle) 
            str_angle2 = self.int_angle_func(angle,self.out_bits)
            # Convert it to binary
            #int_angle = format(int(angle*2**out_bits), 'b')
            '''
            if int_angle == '1' + '0'*out_bits:
                angles[st] = '1'*out_bits # As we only have out_bits, the 10000 is substituted by 1111
            else:
                str_angle = str(int_angle)
                angles[st] = '0'*(out_bits - len(str_angle)) + str_angle
            ''' 
            angles[st] = str_angle2
            #print('dict_key',st)
            #print('energies',self.energies_dictionary[st])
            #print('probability',probability)
            #print('angle',angle)
            #print('angles', angles[st])
        
        # Encoding the new bitmap
        new_bitmap = []
        for o in range(self.out_bits):
            string = ''
            for i in range(int(2**self.in_bits)):
                st = '0'*(self.in_bits - len(str(format(i,'b')))) + str(format(i,'b'))
                string += str(angles[st])[o]
            new_bitmap += [string]
        #print('angles',angles)
        #print('bitmap',new_bitmap)
        self.bitmap = new_bitmap

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