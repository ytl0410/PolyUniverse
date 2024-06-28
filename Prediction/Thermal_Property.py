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

import tensorflow as tf
model_1 = tf.keras.models.load_model('model_tg_MFF.h5')
model_2 = tf.keras.models.load_model('model_tm_MFF.h5')
model_3 = tf.keras.models.load_model('model_td_MFF.h5')

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

Corr_df = pickle.load(open("Corr_Tg.pickle","rb"))
unique_list = pickle.load(open("unique_list_Tg.pickle","rb"))
Columns = pickle.load(open("Columns_Tg.pickle","rb"))
Substructure_list = pickle.load(open("polymer.keys_Tg.pickle","rb"))

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

prediction = model_1.predict(X_MD)
score =  prediction.mean(axis=1)
df['Tg'] = score


prediction = model_2.predict(X_MD)
score =  prediction.mean(axis=1)
df['Tm'] = score


prediction = model_3.predict(X_MD)
score =  prediction.mean(axis=1)
df['Td'] = score

df.to_csv('Prediction_results.csv',index=False)
