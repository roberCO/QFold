import atom
import numpy as np

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
        list_of_atoms_to_rotate += self.decorations_to_rotate(list_of_atoms_to_rotate)

        for atom in list_of_atoms_to_rotate:
            # The rotation angle is pi/16 and in order to maintain internal coherence 
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


    def decorations_to_rotate(self, backbone_list):
        '''Input: list of backbone atoms. Returns a list of other atoms that must be rotated'''
        list_of_atoms_to_rotate = []
        addition_list = []
        for back_atom in backbone_list:
            for atom in back_atom.linked_to:
                if atom.c_type == '' and atom not in list_of_atoms_to_rotate: #if empty check if not nitrogen
                    addition_list += [atom]
                    list_of_atoms_to_rotate += [atom]

            while addition_list != []:
                old_list = addition_list
                addition_list = []
                for prev_atom in old_list:
                    for atom in prev_atom.linked_to:
                        if atom.c_type == '' and atom not in list_of_atoms_to_rotate: #if empty check if not nitrogen
                            addition_list += [atom]
                            list_of_atoms_to_rotate += [atom]

        return list_of_atoms_to_rotate