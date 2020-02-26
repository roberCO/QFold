import qiskit

class QuantumProcessor():
    def __init__(self):
        print('Init qProcessor')

    def inputListOracle(self, inputList):

        truthTableOracle = qiskit.TruthTableOracle(inputList)

