import quantumUtils
import math
import metropolis
import quantumMetropolis

class AngleCalculator():

    def __init__(self, rotationSteps, scaling_factor):

        self.rotationSteps = rotationSteps
        self.scaling_factor = scaling_factor
        self.qTools = quantumUtils.QuantumUtils()

    def calculate3DStructure(self, energyList, option=0):

        #Quantum calculation option for 3D structure
        if option == 0: 

            #HARDCODED
            n_precision_bits = 2
            n_ancilla_bits = 4
            qMetropolis = quantumMetropolis.QuantumMetropolis(n_precision_bits, n_ancilla_bits, energyList)
            return qMetropolis.execute_quantum_metropolis()

        #Classical calculation option for 3D structure
        elif option == 1:

            n_repetitions = 100
            probabilities_matrix = [[0]*self.rotationSteps for x in range(self.rotationSteps)]
            
            for _ in range(n_repetitions):

                n_iterations = 10000
                scaling_factor = 5000 
                classical_metropolis = metropolis.Metropolis(n_iterations, scaling_factor, energyList)
                
                [phi, psi] = classical_metropolis.execute_metropolis()
                probabilities_matrix[phi][psi] += (1/n_repetitions)
                
            return probabilities_matrix
        