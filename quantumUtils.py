import qiskit.aqua.components.oracles as oracle

class QuantumUtils():

    def constructBitMapFromList(self, values):

        #Convert the received values to binary
        binaryValues = []
        for val in values:
            binaryValues.append(self.number2binary(val, 10))

        #Iterate over the binary values to build the bitmap
        bitmap = []
        for x in range(0, 10):

            bitmapEntry = ""
            for y in range(0, len(binaryValues)):
                bitmapEntry += binaryValues[y][x]

            bitmap.append(bitmapEntry)

        return bitmap

    def inputListOracle(self, inputList):

        return oracle.TruthTableOracle(inputList)

    def sortAngleMovement(self, value):

        #The number to sort is composed by values of phi, psi, m (angle 0->phi/1->psi), m (rotation 0->-1/1->1)
        composedNumber = str(value[0]) + str(value[1]) + str(value[2]) + str(value[3])
        return int(composedNumber)

    def sortByAngleMovements(self, listToSort):
        
        listToSort.sort(key=self.sortAngleMovement)

    def number2binary(self, number, numberBits):

        if number < 1 and number != 0:

            #Round HARDCODED
            number = round(number, 3)
            while(number < 1000):
                number *= 10

        numberBinary = "{0:b}".format(int(number))
        numberBinary = numberBinary.zfill(numberBits)

        return numberBinary

