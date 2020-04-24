#####################################################################################################################################
##  Minifold is a project developed by Eric Alcaide (https://github.com/EricAlcaide/MiniFold)                                      ##
##  Thanks to the autor for his awesome project similar to AlphaFold but in miniature                                              ##
##                                                                                                                                 ##
##  Initial project uses Jupyter Notebooks, but in this project everything was translated to plain python for performance reasons  ## 
#####################################################################################################################################

import os

class MinifoldTrainer():
    
    def __init__(self, inputPath):

        self.inputPath = inputPath

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
        f = open('../data/full_under_200.txt', 'w')
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




    def coords_split(lister, splice):
        # Split all passed sequences by "splice" and return an array of them
        # Convert string fragments to float 
        coords = []
        for c in lister
            coords.append([float(a) for a in split(c, splice)])
        
        return coords

    # Could use "Using LinearAlgebra + built-in norm()" but gotta learn Julia
    def norm(vector):
        return sqrt(sum([v*v for v in vector]))