import subprocess
import pandas as pd
import re
import atom

class PsiFour():

    def getAtomsFromProtein(self, protein):
        
        #create input file
        self.createInputFile(protein)

        #execute psi4
        self.executePsiCommand()

        #read/parse outputfile
        return self.parsePsiOutputFile(protein)


    def createInputFile(self, protein):

        inputFile = open('inputFile.dat', 'w')

        inputFile.write('molecule ' + protein + '{\n')
        inputFile.write(' pubchem: '+ protein+'\n')
        inputFile.write('}\n\n')

        inputFile.write('set basis cc-pvdz\n')
        inputFile.write('set reference rhf\n')
        inputFile.write('energy(\'scf\')\n')

        inputFile.close()

    def executePsiCommand(self):

        #execute psi4 by command line (it generates the file output.dat with the information)
        subprocess.run(["psi4", "inputFile.dat", "outputFile.dat"])

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