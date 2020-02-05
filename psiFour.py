import subprocess
import pandas as pd
import re

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

        #Read the output file and write csv with the data
        #file_csv = open(protein + ".csv","w")
        with open("outputFile.dat","r") as filehandle:

            isDataLine = False
            data = []
            for line in filehandle:

                #if line is an empty string after reading data
                if (isDataLine and line.isspace()):
                    break
                
                # Data has ------ and it is necessary to avoid it
                if(isDataLine and not '--' in line):

                    lineChunks = line.split()
                    data.append(lineChunks)

                    #lineFormated = re.sub('\s+', ' ', line).strip()
                    #file_csv.write(lineFormated+'\n') 

                if 'Center' in line:
                    isDataLine = True
                

        #file_csv.close()

        atoms = pd.DataFrame(data, columns = ["Element", "X", "Y", "Z", "Mass"])
        
        return atoms

        #data = pd.read_csv(protein + ".csv", header=None, delim_whitespace=True)
        #data.columns = ["Element", "X", "Y", "Z", "Mass"]
        #print(data.shape)
        #print(type(data))
        
