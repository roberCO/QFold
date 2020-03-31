import quantumUtils
import math

class AngleCalculator():

    def __init__(self, rotationSteps, scaling_factor, beta):

        self.rotationSteps = rotationSteps
        self.scaling_factor = scaling_factor
        self.beta = beta
        self.qTools = quantumUtils.QuantumUtils()

    def calculate3DStructure(self, energyList):

        truthTableList = self.createTruthTableList(energyList)

        #Construct the bitmap
        energyValues = []
        for values in truthTableList:
            energyValues.append(values[4])
        
        binaryInputOracle = self.qTools.constructBitMapFromList(energyValues)

        
        self.qTools.inputListOracle(binaryInputOracle)

    def createTruthTableList(self, energyList):

        truthTableList = []
        for x in range(len(energyList)):
            for y in range(len(energyList[x])):

                energyReference = energyList[x][y]

                #Phi +1
                phiValue = (x+1) % self.rotationSteps
                deltaEnergy = (energyList[phiValue][y] - energyReference) * self.scaling_factor
                probability = min(1, math.exp(-1*deltaEnergy*self.beta))
                truthTableList.append([phiValue, y, 0, 1, probability])

                #Phi -1
                phiValue = (x-1) % self.rotationSteps
                deltaEnergy = (energyList[phiValue][y] - energyReference) * self.scaling_factor
                probability = min(1, math.exp(-1*deltaEnergy*self.beta))
                truthTableList.append([phiValue, y, 0, 0, probability])

                #Psi +1
                psiValue = (y+1) % self.rotationSteps
                deltaEnergy = (energyList[x][psiValue] - energyReference) * self.scaling_factor
                probability = min(1, math.exp(-1*deltaEnergy*self.beta))
                truthTableList.append([x, psiValue, 1, +1, probability])

                #Psi -1
                psiValue = (y-1) % self.rotationSteps
                deltaEnergy = (energyList[x][psiValue] - energyReference) * self.scaling_factor
                probability = min(1, math.exp(-1*deltaEnergy*self.beta))
                truthTableList.append([x, psiValue, 1, 0, probability])

        #sort truthTableList using as key phi, psi and m
        self.qTools.sortByAngleMovements(truthTableList)

        return truthTableList


        