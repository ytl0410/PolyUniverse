import argparse
import pandas as pd

# Create the parser
parser = argparse.ArgumentParser(description='Process some CSV files.')

# Add the arguments
parser.add_argument('input_file', type=str, help='The input CSV file')

# Parse the arguments
args = parser.parse_args()

# Read the input file into a DataFrame
df = pd.read_csv(args.input_file)

import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import DataStructs

import seaborn as sns
from sklearn.model_selection import train_test_split


from rdkit.Chem import rdMolDescriptors as rdmd
from tqdm import tqdm
from functools import wraps

Corr_df = pickle.load(open("Corr_Gas.pickle","rb"))
unique_list = pickle.load(open("unique_list_Gas.pickle","rb"))
Columns = pickle.load(open("Columns_Gas.pickle","rb"))
Substructure_list = pickle.load(open("polymer.keys_Gas.pickle","rb"))

import tensorflow as tf
from tensorflow.python.ops import math_ops
from tensorflow.python.framework import ops
tf.keras.backend.set_floatx('float64')

import numpy as np
import pandas as pd
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.utils import resample
from tensorflow.keras.models import Sequential, save_model, load_model
from tensorflow.keras.layers import Dense, Activation, Dropout
from tensorflow.keras.optimizers import Adam
from tensorflow import keras
import pickle
import argparse
import os
from DNN_functions import nanmean_squared_error, evaluate_model, ensemble_predictions

modelname = 'DNN_BLR_fing'
modeltype = modelname.split('_')[0]
imputation = modelname.split('_')[1]
features = modelname.split('_')[2]

maindirectory = 'pretrained_models/' + modelname

if modeltype == 'DNN':
    folders = os.listdir(maindirectory)
    indices = []
    for name in folders[1:17]:
        #if os.path.isdir(name):
            indices.append(int(name.split('_')[1]))
    max_index = max(indices)

    models = list()
    for i in range(max_index+1):
        directory = maindirectory + '/DNN_' + str(i)
        models.append(tf.keras.models.load_model(directory, custom_objects={'nanmean_squared_error': nanmean_squared_error}))

DatasetA_Smiles_P = pd.read_csv("pretrained_models/datasetA_imputed_all.csv")
DatasetA_grouped = DatasetA_Smiles_P.groupby('Smiles').mean().reset_index()

if imputation == 'BLR':
    Y = DatasetA_grouped.iloc[:,-12:-6]
if imputation == 'ERT':
    Y = DatasetA_grouped.iloc[:,-6:]

Y = np.array(Y)
Yscaler = StandardScaler()
Yscaler.fit(Y)

molecules = df.Smiles.apply(Chem.MolFromSmiles)
fp = molecules.apply(lambda m: AllChem.GetMorganFingerprint(m, radius=3))
fp_n = fp.apply(lambda m: m.GetNonzeroElements())
MY_finger = []
for polymer in fp_n:
    my_finger = [0] * len(unique_list)
    for key in polymer.keys():
        if key in list(Corr_df[0]):
            index = Corr_df[Corr_df[0] == key]['index'].values[0]
            my_finger[index] = polymer[key]         
    MY_finger.append(my_finger)
X_MD = pd.DataFrame(MY_finger)
X_MD = X_MD[Columns]

Y_pred, Y_var = ensemble_predictions(models, X_MD)
Y_pred = Yscaler.inverse_transform(Y_pred)
df[['He', 'H2', 'O2', 'N2', 'CO2', 'CH4']] = Y_pred

df.to_csv('Prediction_results_Gas.csv',index=False)