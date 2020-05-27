import os.path
from keras.models import load_model
import numpy as np
from keras.losses import mean_squared_error, mean_absolute_error

class Minifold:

    def __init__(self, model_path, window_size, max_aa_length):

        self.model_path = model_path+'protein_under_'+str(max_aa_length)+'.h5'
        self.window_size = window_size

        #Hidde info messages from tensorflow
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'

        if not os.path.isfile(self.model_path):

            raise IOError('<!> ERROR: Knowledge model not existing!\nTo generate a model execute: initialAngleTrainer/initialAngleTrainer.py')

    def predictAngles(self, aminoacids):

        print('    ⬤ Loading knowledge model')
        #Load existing model
        model = load_model(self.model_path, custom_objects={'custom_mse_mae': self.custom_mse_mae})

        print('    ⬤ Generating input values')
        input_values = self.generate_input_values(aminoacids)

        #input_aas is a rows x 34 (aas windows size) x 42 (20 aas, Van der Waals distance, Surface, 20 pssms)
        predicted_angles = model.predict(input_values)

        #get angles from sin and cos. Angles are expresed between π and -π
        angles = self.extract_angles(predicted_angles)

        return angles

    # Metric defined in: https://github.com/EricAlcaide/MiniFold/blob/master/models/angles/resnet_1d_angles.py#L14
    def custom_mse_mae(self, y_true, y_pred):
        """ Custom loss function - MSE + MAE """
        return mean_squared_error(y_true, y_pred)+ mean_absolute_error(y_true, y_pred)

    def generate_input_values(self, aminoacids):

        protein_sequence = 'GG'

        input_aas = []

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


        #Method 1

        row_input_aas = []
        protein_sequence_index = 0

        print('<i> Generating input for minifold prediction')
        for index_window in range(0, self.window_size*2):

            column_window = []

            #Add the first aa
            if index_window == self.window_size/2 or index_window == (self.window_size/2)+1:

                #This index stores the position in the key that the aminoacid was found
                key_index = 0

                #in each window is necessary to create an array with 22 positions
                for index_onehotter_vector in range(0, 22):
            
                    #If the aa is the same than in the key, it inserts a 1 (one-hot encoding)
                    if index_onehotter_vector < len(key) and protein_sequence[protein_sequence_index] == key[index_onehotter_vector]:
                        column_window.append(1)
                        #The aminoacid is found in this index value so it is stored in a variable for the VdW and surface
                        key_index = index_onehotter_vector

                    #If the aa is not the same than in the key, it inserts a 0 (one-hot encoding)
                    elif index_onehotter_vector < len(key) and protein_sequence[protein_sequence_index] != key[index_onehotter_vector]:
                        column_window.append(0)

                    #Van der Waals radius
                    elif index_onehotter_vector == 20:
                        column_window.append(vdw_radius[key[key_index]]/max(radius_rel)-basis)
                        
                    #Surface
                    elif index_onehotter_vector == 21:
                        column_window.append(surface[key[key_index]]/max(surface_rel)-surface_basis)

                #The protein was found, so the point is to next protein in the sequence
                protein_sequence_index += 1 

            #Out of the columns 17 and 18, it is necessary to insert 0 padding
            else: 

                for _ in range(0, 22):
                    column_window.append(0)

            row_input_aas.append(column_window)

        input_aas.append(row_input_aas)
        aas = np.array(input_aas)

        return aas

    '''
        #Method 2

        row_input_aas = []
        protein_sequence_index = 0
        for index_window in range(0, 34):

            column_window = []

            #Add the first aa
            if index_window == 17 or index_window == 18:

                #This index stores the position in the key that the aminoacid was found
                key_index = 0

                #in each window is necessary to create an array with 42 positions
                for index_onehotter_vector in range(0, 42):
            
                    #If the aa is the same than in the key, it inserts a 1 (one-hot encoding)
                    if index_onehotter_vector < len(key) and protein_sequence[protein_sequence_index] == key[index_onehotter_vector]:
                        column_window.append(1)
                        #The aminoacid is found in this index value so it is stored in a variable for the VdW and surface
                        key_index = index_onehotter_vector

                    #If the aa is not the same than in the key, it inserts a 0 (one-hot encoding)
                    elif index_onehotter_vector < len(key) and protein_sequence[protein_sequence_index] != key[index_onehotter_vector]:
                        column_window.append(0)

                    #Van der Waals radius
                    elif index_onehotter_vector == 20:
                        column_window.append(vdw_radius[key[key_index]]/max(radius_rel)-basis)
                        
                    #Surface
                    elif index_onehotter_vector == 21:
                        column_window.append(surface[key[key_index]]/max(surface_rel)-surface_basis)
                    
                    #PSSMs that are always with 0 values
                    elif index_onehotter_vector > 21:
                        column_window.append(0)

                #The protein was found, so the point is to next protein in the sequence
                protein_sequence_index += 1 

            #Out of the columns 17 and 18, it is necessary to insert 0 padding
            else: 

                for _ in range(0, 42):
                    column_window.append(0)

            row_input_aas.append(column_window)
            
        input_aas.append(row_input_aas)

        row_input_aas = []
        protein_sequence_index = 0
        for index_window in range(0, 34):

            column_window = []

            #Add the first aa
            if index_window == 18 or index_window == 19:

                #This index stores the position in the key that the aminoacid was found
                key_index = 0

                #in each window is necessary to create an array with 42 positions
                for index_onehotter_vector in range(0, 42):
            
                    #If the aa is the same than in the key, it inserts a 1 (one-hot encoding)
                    if index_onehotter_vector < len(key) and protein_sequence[protein_sequence_index] == key[index_onehotter_vector]:
                        column_window.append(1)
                        #The aminoacid is found in this index value so it is stored in a variable for the VdW and surface
                        key_index = index_onehotter_vector

                    #If the aa is not the same than in the key, it inserts a 0 (one-hot encoding)
                    elif index_onehotter_vector < len(key) and protein_sequence[protein_sequence_index] != key[index_onehotter_vector]:
                        column_window.append(0)

                    #Van der Waals radius
                    elif index_onehotter_vector == 20:
                        column_window.append(vdw_radius[key[key_index]]/max(radius_rel)-basis)
                        
                    #Surface
                    elif index_onehotter_vector == 21:
                        column_window.append(surface[key[key_index]]/max(surface_rel)-surface_basis)
                    
                    #PSSMs that are always with 0 values
                    elif index_onehotter_vector > 21:
                        column_window.append(0)

                #The protein was found, so the point is to next protein in the sequence
                protein_sequence_index += 1 

            #Out of the columns 17 and 18, it is necessary to insert 0 padding
            else: 

                for _ in range(0, 42):
                    column_window.append(0)

            row_input_aas.append(column_window)
            
        input_aas.append(row_input_aas)

        return np.array(input_aas)
    '''
    def extract_angles(self, predicted_angles):

        angles = []

        for prediction in predicted_angles:
            angles_row = []
            phi_sin, phi_cos, psi_sin, psi_cos = prediction[0], prediction[1], prediction[2], prediction[3]

            angles_row.append(np.arctan2(phi_sin, phi_cos))
            angles_row.append(np.arctan2(psi_sin, psi_cos))
            
            angles.append(angles_row)

        return angles