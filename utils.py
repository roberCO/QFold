import atom

class Utils():

    def generateInitialConfig(self, aminoacids):
        print("generateInitialConfig")

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
                if starting_atom.c_type == 'N' and atom.c_type == 'Carboxy':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += backbone_to_rotate(angle_type, atom)

                elif starting_atom.c_type == 'Carboxy' and atom.c_type == 'C_alpha':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += backbone_to_rotate(angle_type, atom)

                elif starting_atom.c_type == 'C_alpha' and atom.c_type == 'N':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += backbone_to_rotate(angle_type, atom)

                #elif atom.c_type == None:
                #list_of_atoms_to_rotate += additional_atoms()
                ## Leave this until we have backbone

        elif angle_type == 'psi': # Follows the structure N -> C_alpha -> Carboxy -> N
            for atom in starting_atom.linked_to:
                if starting_atom.c_type == 'Carboxy' and atom.c_type == 'N':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += backbone_to_rotate(angle_type = angle_type, starting_atom = atom)

                elif starting_atom.c_type == 'C_alpha' and atom.c_type == 'Carboxy':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += backbone_to_rotate(angle_type = angle_type, starting_atom = atom)

                elif starting_atom.c_type == 'N' and atom.c_type == 'C_alpha':
                    list_of_atoms_to_rotate += [atom]
                    list_of_atoms_to_rotate += backbone_to_rotate(angle_type = angle_type, starting_atom = atom)            
        else: 
            raise Exception('Angle_type should be either phi or psi. However it currently is {}'.format(angle_type))
            
        return list_of_atoms_to_rotate


    def decorations_to_rotate(self, backbone_list):
        '''Input: list of backbone atoms. Returns a list of other atoms that must be rotated'''
        list_of_atoms_to_rotate = []
        addition_list = []
        for back_atom in backbone_list:
            for atom in back_atom.linked:
                if atom.c_type == None and atom not in list_of_atoms_to_rotate:
                    addition_list += [atom]
                    list_of_atoms_to_rotate += [atom]

            while addition_list != []:
                old_list = addition_list
                addition_list = []
                for prev_atom in old_list:
                    for atom in prev_atom.linked:
                        if atom.c_type == None and atom not in list_of_atoms_to_rotate:
                            addition_list += [atom]
                            list_of_atoms_to_rotate += [atom]

        return list_of_atoms_to_rotate