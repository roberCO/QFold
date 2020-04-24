#####################################################################################################################################
##  Minifold is a project developed by Eric Alcaide (https://github.com/EricAlcaide/MiniFold)                                      ##
##  Thanks to the autor for his awesome project similar to AlphaFold but in miniature                                              ##
##                                                                                                                                 ##
##  Initial project uses Jupyter Notebooks, but in this project everything was translated to plain python for performance reasons  ## 
#####################################################################################################################################

import os
import numpy as np


class MinifoldTrainer():
    
    def __init__(self, inputPath):

        self.inputPath = inputPath

        #HARDCODED paths
        self.extracted_aminoacids_path = '../data/full_under_200.txt'
        self.full_extracted_aminoacids_path = '../data/angles/full_angles_under_200.txt'
        self.output_path = '../data/angles/outputs.txt'
        self.input_AA_path = '../data/angles/input_aa.txt'
        self.input_PSSM_path = '../data/angles/input_pssm.txt'

        if not os.path.isfile(inputPath):
            raise IOError('<!> ERROR: %s does not exist!' %inputPath)
    
    def train(self):

        #Get protein under 200

        #Get angles froom coords

        #Angle data preparation

        #Predicting angles (uses resnet_1d_angles)



    def getProteinsFromRaw(max_aminoacid_length):

        # Scan first n proteins
        names = []
        seqs = []
        coords = []
        pssms = []

        with open(self.inputPath) as f:
            lines = readlines(f)

            for index in range(1, len(lines)):
                
                if len(coords) == 995:
                    break
                
                # Start recording
                if lines[index] == "[ID]":
                    names.append(lines[index+1])
                
                elif lines[index] == "[PRIMARY]":
                    seqs.append(lines[index+1])

                elif lines[index] == "[TERTIARY]":
                    coords.append(coords_split(lines[index+1:index+3], "\t"))

                elif lines[index] == "[EVOLUTIONARY]":
                    pssms.append(coords_split(lines[index+1:i+21], "\t"))

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
                            aad.append(norm([coords[k][1][i],coords[k][2][i],coords[k][3][i]]-[coords[k][1][j],coords[k][2][j],coords[k][3][j]])))
                    
                    dist.append(aad)
            
            dists.append(dist)
            
        # Data is OK. Save it to a file.
        f = open(self.extracted_aminoacids_path, 'w')
            aux = [0]
            for k in under:
                aux.append(aux[length(aux)]+1)
                # ID
                write(f, "\n[ID]\n")
                write(f, names[k])
                # Seq
                write(f, "\n[PRIMARY]\n")
                write(f, seqs[k])
                # PSSMS
                write(f, "\n[EVOLUTIONARY]\n")
                writedlm(f, pssms[k])
                # Coords
                write(f, "\n[TERTIARY]\n")
                writedlm(f, coords[k])
                # Dists
                write(f, "\n[DIST]\n")
                # Check that saved proteins are less than 200 AAs
                if len(dists[aux[length(aux)]][1])>200:
                    print("error when checking protein in dists n: ", aux[length(aux)], " length: ", length(dists[aux[length(aux)]][1]))
                    break
                else
                    writedlm(f, dists[aux[length(aux)]])


    def getAnglesFromCoords():

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
                pssms.append(parse_line(lines[i+1:i+22]))
            elif line == "[TERTIARY]":
                coords.append(parse_line(lines[i+1:i+3+1]))

        # Organize by atom type
        coords_nterm = [separate_coords(full_coords, 0) for full_coords in coords]
        coords_calpha = [separate_coords(full_coords, 1) for full_coords in coords]
        coords_cterm = [separate_coords(full_coords, 2) for full_coords in coords]

        # Compute angles for a protein
        phis, psis = [], [] # phi always starts with a 0 and psi ends with a 0
        ph_angle_dists, ps_angle_dists = [], []
        for k in range(len(coords)):
            phi, psi = [0.0], []
            # Use our own functions inspired from bioPython
            for i in range(len(coords_calpha[k])):
                # Calculate phi, psi
                # CALCULATE PHI - Can't calculate for first residue
                if i>0:
                    phi.append(get_dihedral(coords_cterm[k][i-1], coords_nterm[k][i], coords_calpha[k][i], coords_cterm[k][i])) # my_calc
                    
                # CALCULATE PSI - Can't calculate for last residue
                if i<len(coords_calpha[k])-1: 
                    psi.append(get_dihedral(coords_nterm[k][i], coords_calpha[k][i], coords_cterm[k][i], coords_nterm[k][i+1])) # my_calc
                
            # Add an extra 0 to psi (unable to claculate angle with next aa)
            psi.append(0)
            # Add protein info to register
            phis.append(phi)
            psis.append(psi)

            # Data is OK. Can save it to file.
            with open(self.full_extracted_aminoacids_path, "w") as f:
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
                        f.write(stringify(pssms[k][j])+"\n")
                    # PHI
                    f.write("\n[PHI]\n")
                    f.write(stringify(phis[k]))
                    # PSI
                    f.write("\n[PSI]\n")
                    f.write(stringify(psis[k]))


    def angle_data_preparation():

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
                pssms.append(parse_lines(lines[i+1:i+22]))
            elif lines[i] == "[PHI]":
                phis.append(parse_line(lines[i+1]))
            elif lines[i] == "[PSI]":
                psis.append(parse_line(lines[i+1]))

        input_aa = []
        input_pssm = []
        outputs = []

        long = 0 # Counter to ensure everythings fine

        for i in range(len(seqs)): 
            if len(seqs[i])>17*2:
                long += len(seqs[i])-17*2
                for j in range(17,len(seqs[i])-17):
                # Padd sequence
                    input_aa.append(onehotter_aa(seqs[i], j))
                    input_pssm.append(pssm_cropper(pssms[i], j))
                    outputs.append([phis[i][j], psis[i][j]])
                    # break

        input_aa = np.array(input_aa).reshape(len(input_aa), 17*2, 22)
        input_pssm = np.array(input_pssm).reshape(len(input_pssm), 17*2, 21)

        # Save outputs to txt file
        with open(self.output_path, "w") as f:
            for o in outputs:
                f.write(stringify(o)+"\n")

        # Save AAs & PSSMs data to different files (together makes a 3dims tensor)
        # Will concat later
        with open(self.input_AA_path, "w") as f:
            for aas in input_aa:
                f.write("\nNEW\n")
                for j in range(len(aas)):
                    f.write(stringify(aas[j])+"\n")

        with open(self.input_PSSM_path, "w") as f:
            for k in range(len(input_pssm)):
                f.write("\nNEW\n")
                for j in range(len(input_pssm[k])):
                    f.write(stringify(input_pssm[k][j])+"\n")

    # Helper function to save data to a .txt file
    def stringify(vec):
        return "".join(str(v)+" " for v in vec)

    # Helper functions to extract numeric data from text
    def parse_lines(raw):
        return np.array([[float(x) for x in line.split(" ") if x != ""] for line in raw])

    def parse_line(line):
        return np.array([float(x) for x in line.split(" ") if x != ""])

    def coords_split(lister, splice):
        # Split all passed sequences by "splice" and return an array of them
        # Convert string fragments to float 
        coords = []
        for c in lister
            coords.append([float(a) for a in split(c, splice)])
        
        return coords

    # Length of masking - 17x2 AAs
    def onehotter_aa(seq, pos):
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
    def pssm_cropper(pssm, pos):
        pssm_out = []
        pad = 17
        for i,row in enumerate(pssm):
            pssm_out.append(row[pos-pad:pos+pad])
        # PSSM is Lx21 - solution: transpose
        return np.array(pssm_out)

    # Could use "Using LinearAlgebra + built-in norm()" but gotta learn Julia
    def norm(vector):
        return sqrt(sum([v*v for v in vector]))

    def parse_line(raw):
        return np.array([[float(x) for x in line.split("\t") if x != ""] for line in raw])

    #Get the coordinates for 1 atom type
    def separate_coords(full_coords, pos): # pos can be either 0(n_term), 1(calpha), 2(cterm)
        res = []
        for i in range(len(full_coords[1])):
            if i%3 == pos:
                res.append([full_coords[j][i] for j in range(3)])

        return np.array(res)

    # Helper functions
    def get_dihedral(coords1, coords2, coords3, coords4):
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

    def stringify(vec):
        """ Helper function to save data to .txt file. """
        line = ""
        for v in vec:
            line = line+str(v)+" "
        return line

    

