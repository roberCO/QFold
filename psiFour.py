import subprocess
import pandas as pd
import re
import atom

class PsiFour():

    def getAtomsFromProtein(self, protein):
        
        fileName = 'inputFile'
        outputFileName = 'outputFile'

        #create input file
        self.createInputFile(protein, fileName)

        #execute psi4
        self.executePsiCommand(fileName, outputFileName)

        #read/parse outputfile
        return self.parsePsiOutputFile(protein)


    def createInputFile(self, protein, filename):

        inputFile = open(filename+'.dat', 'w')

        inputFile.write('molecule ' + protein + '{\n')
        inputFile.write(' pubchem: '+ protein+'\n')
        inputFile.write('}\n\n')

        inputFile.write('set basis cc-pvdz\n')
        inputFile.write('set reference rhf\n')
        inputFile.write('energy(\'scf\')\n')

        inputFile.close()

    def executePsiCommand(self, inputFileName, outputFileName):

        #execute psi4 by command line (it generates the file output.dat with the information)
        subprocess.run(["psi4", inputFileName+".dat", outputFileName+".dat"], stdout=subprocess.DEVNULL)

    def writeFileEnergies(self, atoms, nitroAtom, carboxyAtom, inputFilenameEnergyPSI4):

        #Write file with all atoms rotated
        rotationHandle = open(inputFilenameEnergyPSI4+'.dat', 'w')

        rotationHandle.write('molecule glycylglycine{\n')
        #write input.dat with all rotated atoms
        for at in atoms:
            rotationHandle.write(" " + at.element + " " + str(at.x) + " " + str(at.y) + " " + str(at.z)+'\n')
        
        rotationHandle.write(" " + nitroAtom.element + " " + str(nitroAtom.x) + " " + str(nitroAtom.y) + " " + str(nitroAtom.z)+'\n')
        rotationHandle.write(" " + carboxyAtom.element + " " + str(carboxyAtom.x) + " " + str(carboxyAtom.y) + " " + str(carboxyAtom.z)+'\n')
        rotationHandle.write('}\n\n')
        rotationHandle.write("set basis cc-pvdz\n")
        rotationHandle.write("set reference rhf\n")
        rotationHandle.write("energy('scf')\n")

        rotationHandle.close()

    def readEnergyFromFile(self, outputFilenameEnergyPSI4):

        energy = 0
        with open(outputFilenameEnergyPSI4+'.dat', 'r') as fileHandle:
            for line in fileHandle:

                #If the PSI4 algorithm converges
                if 'Final Energy' in line:
                    energy = float(line.split(':')[1])
                
                #If the PSI4 algorithm does not converge, the energy used is the calculated in the last iteration (iteration 100)
                if 'iter 100:' in line:
                    energy = float(line.split()[3])

        return energy

    def parsePsiOutputFile(self, protein):

        with open("outputFile.dat","r") as filehandle:

            isDataLine = False
            atoms = []
            for line in filehandle:

                #if line is an empty string after reading data
                if (isDataLine and line.isspace()):
                    break
                
                # Data has ------ and it is necessary to avoid it
                if(isDataLine and not '--' in line):

                    lineChunks = line.split()
                    atoms += [atom.Atom(lineChunks[0], float(lineChunks[1]), float(lineChunks[2]), float(lineChunks[3]), float(lineChunks[4]))]

                if 'Center' in line:
                    isDataLine = True
                        
        return atoms