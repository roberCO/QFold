import atom
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import struct
import json
import psiFour
import copy

class Utils():

    def distance(self, atom, atom2):
        return np.sqrt((atom.x-atom2.x)**2+(atom.y-atom2.y)**2+(atom.z-atom2.z)**2)

    def calculateAtomConnection(self, atoms):

        for at1 in atoms:
            for at2 in atoms:
                if at1 != at2:
                    if at1.element == 'O' and at2.element == 'C' and self.distance(at1,at2)<2:
                        at2.c_type = 'Carboxy'
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to
                        
                    if at1.element == 'N' and at2.element == 'C' and self.distance(at1,at2)<2:
                        if at2.c_type != 'Carboxy':
                            at2.c_type = 'C_alpha'
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to
                        
                    if at1.element == 'H' and self.distance(at1,at2)<1.3: 
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to
                        
                    if at1.element == 'C' and at2.element == 'C' and self.distance(at1,at2)<2  and (at1 not in at2.linked_to) : 
                        at1.linked_to = [at2] + at1.linked_to 
                        at2.linked_to = [at1] + at2.linked_to 

        return atoms

    def findAtom(self, atoms, element, cType, connections):

        for at in atoms:

            #The element for the angle phi is nitrogen (N)
            if ((element != '' and at.element == element) or (cType != '' and at.c_type == cType)):

                for conn in connections:

                    linkedElement = conn[0]
                    numberLinkedElements = conn[1]
                    counterLinkedElements = 0

                    #The valid nitrogen is the connected with two carbons
                    for elementConn in at.linked_to:

                        if(elementConn.element == linkedElement):
                            counterLinkedElements += 1

                    #Atom found and saved in variable
                    if (counterLinkedElements == numberLinkedElements):
                        return at

                    else:
                        raise Exception('Element '+at.element+' not found with the proper connections of '+linkedElement+'! '+str(counterLinkedElements)+' found but there should be ' + str(numberLinkedElements))

    def rotate(self, angle_type, angle, starting_atom):
        if angle_type == 'phi':
            if starting_atom.element != 'N':
                raise Exception('Wrong starting atom for the angle phi:',starting_atom.c_type )
            for atom in starting_atom.linked_to:
                if atom.c_type == 'C_alpha':
                    atom_c_alpha = atom
                    
        elif angle_type == 'psi':
            if starting_atom.c_type != 'Carboxy':
                raise Exception('Wrong starting atom for the angle phi:',starting_atom.c_type )
            for atom in starting_atom.linked_to:
                if atom.c_type == 'C_alpha':
                    atom_c_alpha = atom
        else:
            raise Exception('Angle not recognised!:',angle_type)
        
        # Define the list of atoms to rotate and then rotate them
        
        list_of_atoms_to_rotate = []
        list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type,starting_atom)
        list_of_atoms_to_rotate = self.decorations_to_rotate(list_of_atoms_to_rotate,starting_atom)

        for atom in list_of_atoms_to_rotate:
            atom.rotate(atom_c_alpha, starting_atom, angle, angle_type)  

    def backbone_to_rotate(self, angle_type, starting_atom):
        # Define the list of atoms to rotate and then rotate them
        list_of_atoms_to_rotate = []
        if angle_type == 'phi': # Follows the structure N -> Carboxy -> C_alpha -> N
            for atom in starting_atom.linked_to:
                if starting_atom.element == 'N' and atom.c_type == 'Carboxy':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type, atom)

                elif starting_atom.c_type == 'Carboxy' and atom.c_type == 'C_alpha':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type, atom)

                elif starting_atom.c_type == 'C_alpha' and atom.element == 'N':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type, atom)

                #elif atom.c_type == None:
                #list_of_atoms_to_rotate += additional_atoms()
                ## Leave this until we have backbone

        elif angle_type == 'psi': # Follows the structure N -> C_alpha -> Carboxy -> N
            for atom in starting_atom.linked_to:
                if starting_atom.c_type == 'Carboxy' and atom.element == 'N':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type = angle_type, starting_atom = atom)

                elif starting_atom.c_type == 'C_alpha' and atom.c_type == 'Carboxy':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type = angle_type, starting_atom = atom)

                elif starting_atom.element == 'N' and atom.c_type == 'C_alpha':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += self.backbone_to_rotate(angle_type = angle_type, starting_atom = atom)            
        else: 
            raise Exception('Angle_type should be either phi or psi. However it currently is {}'.format(angle_type))
            
        return list_of_atoms_to_rotate


    #Input: list of backbone atoms. Returns a list of all atoms that must be rotated.
    def decorations_to_rotate(self, backbone_list,starting_atom):

        newly_added = backbone_list
        
        while newly_added != []:
            previously_added = newly_added #Don't want the changes made to newly_added to affect to previously added
            newly_added = []
            for prev_atom in previously_added:
                for atom in prev_atom.linked_to:
                    if atom not in backbone_list and atom not in newly_added and atom != starting_atom:
                        newly_added += [atom]
            backbone_list += newly_added
            previously_added = []
        
        for atom in starting_atom.linked_to:
            if atom.element != 'N' and atom.c_type != 'Carboxy' and atom.c_type != 'C_alpha':
                backbone_list += [atom]
        
        #Perhaps use a method to eliminate duplicate references in the list? I think there shouldn't be any.
        return backbone_list

    def plotting(self, list_of_atoms, title):
    
        #-----
        VecStart_x = []
        VecStart_y = []
        VecStart_z = []
        VecEnd_x = []
        VecEnd_y = []
        VecEnd_z  = []
        
        #Make list of conections
        list_of_connections = []
        for at1 in list_of_atoms:
            for at2 in list_of_atoms:
                if (at1,at2) not in list_of_connections and (at2,at1) not in list_of_connections and at2 in at1.linked_to:
                    list_of_connections += [(at1,at2)]    
                    
        for tupl in list_of_connections:
            VecStart_x += [tupl[0].x]
            VecStart_y += [tupl[0].y]
            VecStart_z += [tupl[0].z]
            VecEnd_x += [tupl[1].x]
            VecEnd_y += [tupl[1].y]
            VecEnd_z  += [tupl[1].z]
            
        fig = plt.figure()
        fig.canvas.set_window_title(title)
        ax = fig.add_subplot(111, projection='3d')
        
        for i in range(len(list_of_connections)):
            ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]],color='grey')
        #-----    
        
        xs = []
        ys = []
        zs = []
        c = []

        planePhiPoints = []
        planePsiPoints = []

        for at in list_of_atoms:
            xs += [at.x]
            ys += [at.y]
            zs += [at.z]
            if at.element == 'N':
                c += ['blue']
            elif at.element == 'C':
                c += ['black']
            elif at.element == 'O':
                c += ['red']
            elif at.element == 'H':
                c += ['green']  

            if at.atomId == 3 or at.atomId == 5 or at.atomId == 6:
                planePhiPoints += [np.array([at.x, at.y, at.z])]
            
            if at.atomId == 4 or at.atomId == 6 or at.atomId == 7:
                planePsiPoints += [np.array([at.x, at.y, at.z])]

        ax.scatter(xs, ys, zs,c=c,depthshade= False)


        XPhi, YPhi, ZPhi = self.calculatePlane(planePhiPoints)
        XPsi, YPsi, ZPsi = self.calculatePlane(planePsiPoints)


        # plot the mesh. Each array is 2D, so we flatten them to 1D arrays
        ax.plot_surface(XPhi, YPhi, ZPhi, color = 'r', alpha = '0.5')
        ax.plot_surface(XPsi, YPsi, ZPsi, color = 'g', alpha = '0.5')


        for i in range(len(xs)): 
            ax.text(xs[i],ys[i],zs[i],  '%s' % (str(i)))
        plt.show()

    def calculateEnergies(self, proteinName, numberBitsRotation):

        energiesJson = {}
        energiesJson['protein'] = proteinName
        energiesJson['numberBitsRotation'] = numberBitsRotation
        energiesJson['energies'] = []
        rotationSteps = pow(2, int(numberBitsRotation))

        #call psi4 to get the atoms of the protein
        psi = psiFour.PsiFour()
        atoms = psi.getAtomsFromProtein(proteinName)

        #Calculate the connection between atoms
        atoms = self.calculateAtomConnection(atoms)

        nitroConnections = [['C', 2]]
        carboxyConnections = [['C', 1], ['O', 2]]

        nitroAtom = self.findAtom(atoms, 'N', '', nitroConnections)
        carboxyAtom = self.findAtom(atoms, '', 'Carboxy', carboxyConnections)

        inputFilenameEnergyPSI4 = 'inputRotations'
        outputFilenameEnergyPSI4 = 'outputRotations'

        anglePhi = 1/rotationSteps
        anglePsi = 1/rotationSteps

        anglesEnergy = []
        #These two nested loops are hardcoded (it could be n nested loops, 1 per AA) because QFold is going to be used just with two and three aminoacids
        #if it is scale to more aminoacids, it should be necessary to implement a recursive function
        for x in range(0, rotationSteps):

            for y in range(0, rotationSteps):

                #Perform the rotations over a copy
                copied_atoms = copy.deepcopy(atoms)
                copied_nitroAtom = self.findAtom(copied_atoms, 'N', '', nitroConnections)
                copied_carboxyAtom = self.findAtom(copied_atoms, '', 'Carboxy', carboxyConnections)

                #Always rotate from state (0,0)
                self.rotate('phi', x * anglePhi, copied_nitroAtom) 

                self.rotate('psi', y * anglePsi, copied_carboxyAtom)
                
                #Write the file with the actual rotations
                psi.writeFileEnergies(copied_atoms, inputFilenameEnergyPSI4)

                #Calculate the energy of the actual rotations using PSI4
                psi.executePsiCommand(inputFilenameEnergyPSI4, outputFilenameEnergyPSI4)

                #Read the PSI4 output file and get the energy
                energy = psi.readEnergyFromFile(outputFilenameEnergyPSI4)
                normalizedEnergy = energy*-1
                normalizedEnergy = "%.6f" % normalizedEnergy

                anglesEnergy.append([anglePhi, anglePsi, normalizedEnergy])
                anglePsi += 1/rotationSteps

                #print ('⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤\n⬤   Phi('+str(x)+'): ' + str(anglePhi) +'\n⬤   Psi('+str(y)+'): '+ str(anglePsi)+ '\n⬤   Energy: ' + str(energy) +'\n⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤⬤\n\n')
                energiesJson['energies'].append({
                    'phi': x,
                    'psi': y,
                    'energy': energy,
                })

                # We eliminate previous copies
                del copied_atoms
                del copied_carboxyAtom
                del copied_nitroAtom

            anglePhi += 1/rotationSteps

        #Create json with calculated energies
        with open('./precalculated_energies/energies_'+proteinName+'_'+str(numberBitsRotation)+'.json', 'w') as outfile:
            json.dump(energiesJson, outfile)

    def readEnergyJson(self, proteinName, numberBitsRotation):

        rotationSteps = pow(2, int(numberBitsRotation))

        with open('./precalculated_energies/energies_'+proteinName+'_'+str(numberBitsRotation)+'.json') as json_file:
            data = json.load(json_file)

            #Create an empty memory structure with the rotation steps dimension
            energyList = [[0 for x in range(rotationSteps)] for y in range(rotationSteps)]

            #Get each entry with energies of the json
            for angle in data['energies']:

                x = angle['phi']
                y = angle['psi']
                energyList[x][y] = angle['energy']

            return energyList

    def calculatePlane(self, planePoints):

        # These two vectors are in the plane
        v1 = planePoints[2] - planePoints[0]
        v2 = planePoints[1] - planePoints[0]

        # the cross product is a vector normal to the plane
        cp = np.cross(v1, v2)
        a, b, c = cp

        # This evaluates a * x3 + b * y3 + c * z3 which equals d
        d = np.dot(cp, planePoints[2])


        #HARCODED values
        X, Y = np.meshgrid(planePoints[0], planePoints[0])

        Z = (d - a * X - b * Y) / c

        return X, Y, Z

    def number2binary(self, number, numberBits):

        if number < 1 and number != 0:

            #Round HARDCODED
            number = round(number, 3)
            while(number < 1000):
                number *= 10

        numberBinary = "{0:b}".format(int(number))
        numberBinary = numberBinary.zfill(numberBits)

        return numberBinary

    def sortAngleMovement(self, value):

        #the number to sort is composed by values of phi, psi, m (angle 0->phi/1->psi), m (rotation 0->-1/1->1)
        composedNumber = str(value[0]) + str(value[1]) + str(value[2]) + str(value[3])
        return int(composedNumber)

    def sortByAngleMovements(self, listToSort):
        listToSort.sort(key=self.sortAngleMovement)