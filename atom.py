import numpy as np
import warnings
from scipy.spatial.transform import Rotation as R


class Atom():
    def __init__(self,element,x,y,z,mass,aminoacid = "" ,c_type = "",linked_to = [],vanderWals_r = 0):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.aminoacid = aminoacid
        self.c_type = c_type
        self.linked_to = linked_to
        self.vanderWalls_r = 0
        
    def rotate(self,atom1,atom2,angle,angle_type): # Angle should be a fraction since it will be multiplied by 2pi
        
        # Check whether the angle is between the right atoms. Else raise exception
        if angle_type == 'phi':
            if atom1.c_type != 'C_alpha' or atom2.element != 'N' :
                raise Exception('Not the correct angle between N and C_alpha. The atoms are ',atom1.element,' and ',atom2.c_type)
        elif angle_type == 'psi':
            if atom1.c_type != 'C_alpha' or atom2.c_type != 'Carboxy':
                raise Exception('Not the correct angle between C_alpha and Carboxy. The atoms are ',atom1.c_type,' and ',atom2.c_type)
        else: 
            raise Exception('The angle is not phi or psi!')
            
        v1 = np.array([atom1.x,atom1.y,atom1.z])
        v2 = np.array([atom2.x,atom2.y,atom2.z])
        v = v2-v1
        if np.linalg.norm(v) < 0.001:
            warnings.warn('Are you sure atom1 and atom2 are different? They seem to have the same position')
        v /= np.linalg.norm(v)
        
        r = R.from_rotvec(2*np.pi*angle * v)
        M = r.as_matrix()
        
        #Firs thing we set a new coordinate system where v1 = [0,0,0] (translation only) so that we can perform the rotation of M
        p = np.array([self.x,self.y,self.z])
        p -= v1
        
        
        #Next we perform the rotation along the axis 
        p = M@p
        
        #Undo the translation
        p += v1
        self.x = p[0]
        self.y = p[1]
        self.z = p[2]