#####################################################################################################################################
##  Minifold is a project developed by Eric Alcaide (https://github.com/EricAlcaide/MiniFold)                                      ##
##  Thanks to the autor for his awesome project similar to AlphaFold but in miniature                                              ##
##                                                                                                                                 ##
##  Initial project uses Jupyter Notebooks, but in this project everything was translated to plain python for performance reasons  ## 
#####################################################################################################################################

import os
import numpy as np
import matplotlib.pyplot as plt
import math

# Import libraries
import keras
import keras.backend as K
from keras.models import Model
# Optimizer and regularization
from keras.regularizers import l2
from keras.losses import mean_squared_error, mean_absolute_error
# Keras layers
from keras.layers.convolutional import Conv1D
from keras.layers import Dense, Dropout, Flatten, Input, BatchNormalization, Activation
from keras.layers.pooling import MaxPooling1D, AveragePooling1D, MaxPooling2D, AveragePooling2D
from keras.optimizers import Adam

# Model architecture
from resnet_1d_angles import resnet_v2, custom_mse_mae

class MinifoldTrainer():
    
    def __init__(self, inputPath):

        self.inputPath = inputPath

        #HARDCODED paths
        self.extracted_aminoacids_path = './data/full_under_200.txt'
        self.full_extracted_aminoacids_path = './data/full_angles_under_200.txt'
        self.output_path = './data/outputs.txt'
        self.input_AA_path = './data/input_aa.txt'
        self.input_PSSM_path = './data/input_pssm.txt'

        #HARDCODED parameters
        self.maximum_aminoacid_length = 200

        if not os.path.isfile(inputPath):
            raise IOError('<!> ERROR: %s does not exist!' %inputPath)

        #Create the path to store the data
        if not os.path.exists('./data/'):
            os.mkdir('./data/')
    
    def train(self):

        #Get protein under defined maximum aminoacid length
        self.getProteinsFromRaw(self.maximum_aminoacid_length)

        #Get angles froom coords
        self.getAnglesFromCoords()

        #Angle data preparation
        self.angleDataPreparation()

        #Generate model (using resnet_1d_angles)
        self.generateModel()

    def getProteinsFromRaw(self, max_aminoacid_length):

        # Scan first n proteins
        names = []
        seqs = []
        coords = []
        pssms = []

        with open(self.inputPath) as f:
            lines = f.readlines()

            for index in range(1, len(lines)):
                
                if len(coords) == 995:
                    break
                
                # Start recording
                if lines[index] == "[ID]\n":
                    names.append(lines[index+1])
                
                elif lines[index] == "[PRIMARY]\n":
                    seqs.append(lines[index+1])

                elif lines[index] == "[TERTIARY]\n":
                    coords.append(self.coords_split(lines[index+1:index+3], "\t"))

                elif lines[index] == "[EVOLUTIONARY]\n":
                    pssms.append(self.coords_split(lines[index+1:index+21], "\t"))

        under = []
        for index in range (1, len(seqs)):
            if len(seqs[index]) < max_aminoacid_length:
                under.append(index)

        dists = []
        # Get distances btwn pairs of AAs - only for prots under 200
        for k in under:
            # Get distances from coordinates
            dist = []
            for i in range(1, len(coords[k][1])):
                # Only pick coords for C-alpha carbons! - position (1/3 of total data)
                # i%3 == 2 Because juia arrays start at 1 - Python: i%3 == 1
                if i%3 == 2:
                    aad = [] # Distance to every AA from a given AA
                    for j in range(1, len(coords[k][1])):
                        if j%3 == 2:
                            aad.append(self.norm([coords[k][1][i],coords[k][2][i],coords[k][3][i]]-[coords[k][1][j],coords[k][2][j],coords[k][3][j]]))
                    
                    dist.append(aad)
            
            dists.append(dist)
            
        # Data is OK. Save it to a file.       
        with open(self.extracted_aminoacids_path, 'w+') as f:
            aux = [0]
            for k in under:
                aux.append(aux[len(aux)]+1)
                # ID
                f.write("\n[ID]\n")
                f.write(names[k])
                # Seq
                f.write("\n[PRIMARY]\n")
                f.write(seqs[k])
                # PSSMS
                f.write("\n[EVOLUTIONARY]\n")
                f.writedlm(pssms[k])
                # Coords
                f.write("\n[TERTIARY]\n")
                f.writedlm(coords[k])
                # Dists
                f.write("\n[DIST]\n")
                # Check that saved proteins are less than 200 AAs
                if len(dists[aux[len(aux)]][1])>200:
                    print("error when checking protein in dists n: ", aux[len(aux)], " length: ", len(dists[aux[len(aux)]][1]))
                    break
                else:
                    f.writedlm(dists[aux[len(aux)]])

            print('<i> File ', self.extracted_aminoacids_path, 'written')

    def getAnglesFromCoords(self):

        names = []
        seqs = []
        psis = []
        phis = []
        pssms = []
        coords = []

        # Opn file and read text
        with open(self.extracted_aminoacids_path, "r") as f:
            lines = f.read().split('\n')

        # Extract numeric data from text
        for i,line in enumerate(lines):
            if len(names) == 601:
                break

            # Read each protein separately
            if line == "[ID]":
                names.append(lines[i+1])
            elif line == "[PRIMARY]":
                seqs.append(lines[i+1])
            elif line == "[EVOLUTIONARY]":
                pssms.append(self.parse_line_angle_from_coords(lines[i+1:i+22]))
            elif line == "[TERTIARY]":
                coords.append(self.parse_line_angle_from_coords(lines[i+1:i+3+1]))

        # Organize by atom type
        coords_nterm = [self.separate_coords(full_coords, 0) for full_coords in coords]
        coords_calpha = [self.separate_coords(full_coords, 1) for full_coords in coords]
        coords_cterm = [self.separate_coords(full_coords, 2) for full_coords in coords]

        # Compute angles for a protein
        phis, psis = [], [] # phi always starts with a 0 and psi ends with a 0
        for k in range(len(coords)):
            phi, psi = [0.0], []
            # Use our own functions inspired from bioPython
            for i in range(len(coords_calpha[k])):
                # Calculate phi, psi
                # CALCULATE PHI - Can't calculate for first residue
                if i>0:
                    phi.append(self.get_dihedral(coords_cterm[k][i-1], coords_nterm[k][i], coords_calpha[k][i], coords_cterm[k][i])) # my_calc
                    
                # CALCULATE PSI - Can't calculate for last residue
                if i<len(coords_calpha[k])-1: 
                    psi.append(self.get_dihedral(coords_nterm[k][i], coords_calpha[k][i], coords_cterm[k][i], coords_nterm[k][i+1])) # my_calc
                
            # Add an extra 0 to psi (unable to claculate angle with next aa)
            psi.append(0)
            # Add protein info to register
            phis.append(phi)
            psis.append(psi)

            # Data is OK. Can save it to file.
            with open(self.full_extracted_aminoacids_path, 'w+') as f:
                for k in range(len(names)-1):
                    # ID
                    f.write("\n[ID]\n")
                    f.write(names[k])
                    # Seq
                    f.write("\n[PRIMARY]\n")
                    f.write(seqs[k])
                    # PSSMS
                    f.write("\n[EVOLUTIONARY]\n")
                    for j in range(len(pssms[k])):
                        f.write(self.stringify_angle_from_coords(pssms[k][j])+"\n")
                    # PHI
                    f.write("\n[PHI]\n")
                    f.write(self.stringify_angle_from_coords(phis[k]))
                    # PSI
                    f.write("\n[PSI]\n")
                    f.write(self.stringify_angle_from_coords(psis[k]))

                print('<i> File ', self.full_extracted_aminoacids_path, 'written')
        


    def angleDataPreparation(self):

        # Opn file and read text
        with open(self.full_extracted_aminoacids_path, "r") as f:
            lines = f.read().split('\n')

        # Scan first n proteins
        names = []
        seqs = []
        psis = []
        phis = []
        pssms = []

        # Extract numeric data from text
        for i,line in enumerate(lines):
            if len(names) == 601:
                break
            # Read each protein separately
            if line == "[ID]":
                names.append(lines[i+1])
            elif line == "[PRIMARY]":
                seqs.append(lines[i+1])
            elif line == "[EVOLUTIONARY]":
                pssms.append(self.parse_lines(lines[i+1:i+22]))
            elif lines[i] == "[PHI]":
                phis.append(self.parse_line_angle_data_preparation(lines[i+1]))
            elif lines[i] == "[PSI]":
                psis.append(self.parse_line_angle_data_preparation(lines[i+1]))

        input_aa = []
        input_pssm = []
        outputs = []

        long = 0 # Counter to ensure everythings fine

        for i in range(len(seqs)): 
            if len(seqs[i])>17*2:
                long += len(seqs[i])-17*2
                for j in range(17,len(seqs[i])-17):
                # Padd sequence
                    input_aa.append(self.onehotter_aa(seqs[i], j))
                    input_pssm.append(self.pssm_cropper(pssms[i], j))
                    outputs.append([phis[i][j], psis[i][j]])
                    # break

        input_aa = np.array(input_aa).reshape(len(input_aa), 17*2, 22)
        input_pssm = np.array(input_pssm).reshape(len(input_pssm), 17*2, 21)

        # Save outputs to txt file
        with open(self.output_path, "w+") as f:
            for o in outputs:
                f.write(self.stringify_angle_data_preparation(o)+"\n")
            
            print('<i> File ', self.output_path, 'written')


        # Save AAs & PSSMs data to different files (together makes a 3dims tensor)
        # Will concat later
        with open(self.input_AA_path, "w+") as f:
            for aas in input_aa:
                f.write("\nNEW\n")
                for j in range(len(aas)):
                    f.write(self.stringify_angle_data_preparation(aas[j])+"\n")
            
            print('<i> File ', self.input_AA_path, 'written')

        with open(self.input_PSSM_path, "w+") as f:
            for k in range(len(input_pssm)):
                f.write("\nNEW\n")
                for j in range(len(input_pssm[k])):
                    f.write(self.stringify_angle_data_preparation(input_pssm[k][j])+"\n")
            
            print('<i> File ', self.input_PSSM_path, 'written')

    def generateModel(self):

        ## LOAD DATASET ##

        # Load outputs/labels from file
        outputs = np.genfromtxt(self.output_path)

        # Get inputs data
        aas = self.get_ins()
        pssms = self.get_ins(pssm=True)

        ## REMOVE nan VALUES FROM DATASET
        #If there is a nan in the outputs list is necessary to discard the row in the aa and pssm list
        outputs_no_nan = []
        aas_no_nan = []
        pssms_no_nan = []

        for index in range(0, len(outputs)):
            if(str(outputs[index][0]) != 'nan'):
                outputs_no_nan.append(outputs[index])
                aas_no_nan.append(aas[index])
                pssms_no_nan.append(pssms[index])
            
        outputs = np.array(outputs_no_nan)
        aas = np.array(aas_no_nan)
        pssms = np.array(pssms_no_nan)

        #######

        out = []
        out.append(np.sin(outputs[:,0]))
        out.append(np.cos(outputs[:,0]))
        out.append(np.sin(outputs[:,1]))
        out.append(np.cos(outputs[:,1]))
        out = np.array(out).T
        print('out shape: ', out.shape)


        # Concatenate input features
        inputs = np.concatenate((aas[:, :, :20], pssms[:, :, :20], aas[:, :, 20:]), axis=2) 
        np.set_printoptions(threshold=np.inf)

        # Separate data between training and testing
        split = 38700
        x_train, x_test = inputs[:split], inputs[split:]
        y_train, y_test = out[:split], out[split:]

        ## LOADING MODEL ##

        # Using AMSGrad optimizer for speed 
        adam = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, decay=0.0, amsgrad=True)
        # Create model
        model = resnet_v2(input_shape=(17*2,42), depth=20, num_classes=4, conv_first=True)
        model.compile(optimizer=adam, loss=custom_mse_mae, metrics=["mean_absolute_error", "mean_squared_error"])

        # Resnet (pre-act structure) with 34*42 columns as inputs - leaving a subset for validation
        model.fit(x_train, y_train, epochs=5, batch_size=16, verbose=1, shuffle=True, validation_data=(x_test, y_test))


    # Helper function to save data to a .txt file
    def stringify_angle_data_preparation(self, vec):
        return "".join(str(v)+" " for v in vec)

    # Helper functions to extract numeric data from text
    def parse_lines(self, raw):
        return np.array([[float(x) for x in line.split(" ") if x != ""] for line in raw])

    def parse_line_angle_data_preparation(self, line):
        return np.array([float(x) for x in line.split(" ") if x != ""])

    def coords_split(self, lister, splice):
        # Split all passed sequences by "splice" and return an array of them
        # Convert string fragments to float 
        coords = []
        for c in lister:
            coords.append([float(a) for a in c.split(splice)])
        
        return coords

    # Length of masking - 17x2 AAs
    def onehotter_aa(self, seq, pos):
        pad = 17
        # Pad sequence
        key = "HRKDENQSYTCPAVLIGFWM"
        # Van der Waals radius
        vdw_radius = {"H": 118, "R": 148, "K": 135, "D": 91, "E": 109, "N": 96, "Q": 114,
                    "S": 73, "Y": 141, "T": 93, "C": 86, "P": 90, "A": 67, "V": 105,
                    "L": 124, "I": 124, "G": 48, "F": 135, "W": 163, "M": 124}
        radius_rel = vdw_radius.values()
        basis = min(radius_rel)/max(radius_rel)
        # Surface exposure 
        surface = {"H": 151, "R": 196, "K": 167, "D": 106, "E": 138, "N": 113, "Q": 144,
                    "S": 80, "Y": 187, "T": 102, "C": 104, "P": 105, "A": 67, "V": 117,
                    "L": 137, "I": 140, "G": 0, "F": 175, "W": 217, "M": 160}
        surface_rel = surface.values()
        surface_basis = min(surface_rel)/max(surface_rel)
        # One-hot encoding
        one_hot = []
        for i in range(pos-pad, pos+pad): # alponer los guiones ya tiramos la seq para un lado
            vec = [0 for i in range(22)]
            # mark as 1 the corresponding indexes
            for j in range(len(key)):
                if seq[i] == key[j]:
                    vec[j] = 1
                    # Add Van der Waals relative radius
                    vec[-2] = vdw_radius[key[j]]/max(radius_rel)-basis
                    vec[-1] = surface[key[j]]/max(surface_rel)-surface_basis

                    print('AA: ', key[j], ' VDW radius: ', vec[-2], ' a: ', vdw_radius[key[j]], ' b: ', max(radius_rel), ' c: ', basis)
                    print('AA: ', key[j], ' surface:    ', vec[-1], ' a: ', surface[key[j]], ' b: ', max(surface_rel), ' c: ', surface_basis, '\n')

            
            one_hot.append(vec) 
        
        return np.array(one_hot)

    #Crops the PSSM matrix
    def pssm_cropper(self, pssm, pos):
        pssm_out = []
        pad = 17
        for row in enumerate(pssm):
            pssm_out.append(row[pos-pad:pos+pad])
        # PSSM is Lx21 - solution: transpose
        return np.array(pssm_out)

    # Could use "Using LinearAlgebra + built-in norm()" but gotta learn Julia
    def norm(self, vector):
        return math.sqrt(sum([v*v for v in vector]))

    def parse_line_angle_from_coords(self, raw):
        return np.array([[float(x) for x in line.split("\t") if x != ""] for line in raw])

    #Get the coordinates for 1 atom type
    def separate_coords(self, full_coords, pos): # pos can be either 0(n_term), 1(calpha), 2(cterm)
        res = []
        for i in range(len(full_coords[1])):
            if i%3 == pos:
                res.append([full_coords[j][i] for j in range(3)])

        return np.array(res)

    # Helper functions
    def get_dihedral(self, coords1, coords2, coords3, coords4):
        """Returns the dihedral angle in degrees."""

        a1 = coords2 - coords1
        a2 = coords3 - coords2
        a3 = coords4 - coords3

        v1 = np.cross(a1, a2)
        v1 = v1 / (v1 * v1).sum(-1)**0.5
        v2 = np.cross(a2, a3)
        v2 = v2 / (v2 * v2).sum(-1)**0.5
        porm = np.sign((v1 * a3).sum(-1))
        rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
        if not porm == 0:
            rad = rad * porm

        return rad

    def stringify_angle_from_coords(self, vec):
        """ Helper function to save data to .txt file. """
        line = ""
        for v in vec:
            line = line+str(v)+" "
        return line

    def get_ins(self, path = 'default', pssm=None):

        if (path == 'default'):
            path = self.input_AA_path

        """ Gets inputs from both AminoAcids (input_aa) and PSSM (input_pssm)"""
        # handles both files
        if pssm: path = self.input_PSSM_path
        # Open file and read text
        with open(path, "r") as f:
            lines = f.read().split('\n')
        # Extract numeric data from text
        pre_ins = []
        for i,line in enumerate(lines):
            # Read each protein separately
            if line == "NEW":
                prot = []
                raw = lines[i+1:i+(17*2+1)]
                # Read each line as a vector + ensemble one-hot vectors as a matrix
                for r in raw:
                    prot.append(np.array([float(x) for x in r.split(" ") if x != ""]))
                # Add prot to dataset
                pre_ins.append(np.array(prot))  
        
        return np.array(pre_ins)