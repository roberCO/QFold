import qiskit.aqua.components.oracles as oracle

class QuantumProcessor():
    def __init__(self):
        print('Init qProcessor')

    def inputListOracle(self, inputList):

        truthTableOracle = oracle.TruthTableOracle(inputList)

