from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd

df_1 = pd.read_csv('epoxy_diN.csv')
df_2 = pd.read_csv('epoxy_diE.csv')

def get_neiid_bysymbol(mol,marker):
    try:
        for atom in mol.GetAtoms():
            if atom.GetSymbol()==marker:
                neighbors=atom.GetNeighbors()
                if len(neighbors)>1:
                    print ('Cannot process more than one neighbor, will only return one of them')
                atom_nb=neighbors[0]
                return atom_nb.GetIdx()
    except Exception as e:
        print (e)
        return None

def get_id_bysymbol(mol,marker):
    for atom in mol.GetAtoms():
        if atom.GetSymbol()==marker:
            return atom.GetIdx()
        
def combine2frags(mol_a,mol_b,maker_b='Cs',maker_a='Rb'):
    merged_mol = Chem.CombineMols(mol_a,mol_b)
    bind_pos_a=get_neiid_bysymbol(merged_mol,maker_a)
    bind_pos_b=get_neiid_bysymbol(merged_mol,maker_b)
    ed_merged_mol= Chem.EditableMol(merged_mol)
    ed_merged_mol.AddBond(bind_pos_a,bind_pos_b,order=Chem.rdchem.BondType.SINGLE)
    marker_a_idx=get_id_bysymbol(merged_mol,maker_a)
    ed_merged_mol.RemoveAtom(marker_a_idx)
    temp_mol = ed_merged_mol.GetMol()
    marker_b_idx=get_id_bysymbol(temp_mol,maker_b)
    ed_merged_mol=Chem.EditableMol(temp_mol)
    ed_merged_mol.RemoveAtom(marker_b_idx)
    final_mol = ed_merged_mol.GetMol()
    final_smiles = Chem.MolToSmiles(final_mol)
    return final_smiles

substructure_smarts = '[NX3;H2].[NX3;H2]'

def contains_substructure(smiles, substructure):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol.HasSubstructMatch(substructure)
    else:
        return False 

tepo_smiles = []

substructure_molecule = Chem.MolFromSmarts(substructure_smarts)

for smiles in df_1['Smiles']:
    if contains_substructure(smiles, substructure_molecule):
        tepo_smiles.append(smiles)

tepo = pd.DataFrame({'Smiles': tepo_smiles})

new_smiles_list = []
for smiles in df_1['Smiles']:
    # 创建一个Molecule对象
    m = Chem.MolFromSmiles(smiles)
    
    patt = Chem.MolFromSmarts('[NX3;H2]')
    repl1 = Chem.MolFromSmiles('[Rb]')
    repl2 = Chem.MolFromSmiles('[Cs]')
    rms = AllChem.ReplaceSubstructs(m, patt, repl1, replacementConnectionPoint=0)
    rms[0]
    rms = AllChem.ReplaceSubstructs(rms[0], patt, repl2, replacementConnectionPoint=0)
    modified_smiles = Chem.MolToSmiles(rms[0])
    modified_smiles
    mol1=Chem.MolFromSmiles('[Sr]N[*]')
    mol2=Chem.MolFromSmiles(modified_smiles)
    modified_smiles = combine2frags(mol1,mol2,'Rb','Sr')
    modified_smiles
    mol1=Chem.MolFromSmiles('[Sr]NC[Rb]')
    mol2=Chem.MolFromSmiles(modified_smiles)
    modified_smiles = combine2frags(mol1,mol2,'Cs','Sr')
    modified_smiles
    
    new_smiles_list.append(modified_smiles)


df_1['New_Smiles'] = new_smiles_list

new_smiles_list = []

for smiles in df_2['Smiles']:
    try:
        m = Chem.MolFromSmiles(smiles)
        patt = Chem.MolFromSmarts('C1OC1')
        repl1 = Chem.MolFromSmiles('[Rb]')
        repl2 = Chem.MolFromSmiles('[Cs]')
        rms = AllChem.ReplaceSubstructs(m, patt, repl1, replacementConnectionPoint=0)
        rms[0]
        rms = AllChem.ReplaceSubstructs(rms[0], patt, repl2, replacementConnectionPoint=0)
        modified_smiles = Chem.MolToSmiles(rms[0])
        modified_smiles
        mol1 = Chem.MolFromSmiles('C(C(O)[Sr])[*]')
        mol2 = Chem.MolFromSmiles(modified_smiles)
        modified_smiles = combine2frags(mol1, mol2, 'Cs', 'Sr')
        modified_smiles
        mol1 = Chem.MolFromSmiles('C(C(O)[Sr])[Cs]')
        mol2 = Chem.MolFromSmiles(modified_smiles)
        modified_smiles = combine2frags(mol1, mol2, 'Rb', 'Sr')
        modified_smiles
        new_smiles_list.append(modified_smiles)
    except Exception as e:
        print(f"Error processing SMILES: {e}")
        new_smiles_list.append('False')
df_2['New_Smiles'] = new_smiles_list

new_df = pd.DataFrame(columns=['Smiles', 'Smiles_Compound_1', 'Smiles_Compound_2'])

for smiles_1 in df_1['New_Smiles']:
    mol1 = Chem.MolFromSmiles(smiles_1)
    for smiles_2 in df_2['New_Smiles']:
        if smiles_2 is not False: 
            mol2 = Chem.MolFromSmiles(smiles_2)
            if mol2 is None:
                    continue
            new_smiles = combine2frags(mol1, mol2)
            
            new_df = new_df.append({'Smiles': new_smiles, 'Smiles_Compound_1': smiles_1, 'Smiles_Compound_2': smiles_2}, ignore_index=True)

for index, row in new_df.iterrows():
    smiles1 = row['Smiles_Compound_1']
    matching_row_df_1 = df_1[df_1['New_Smiles'] == smiles1]
    
    if not matching_row_df_1.empty:
        new_df.at[index, 'Smiles_Compound_1'] = matching_row_df_1.iloc[0]['Smiles']

for index, row in new_df.iterrows():
    smiles2 = row['Smiles_Compound_2']
    matching_row_df_2 = df_2[df_2['New_Smiles'] == smiles2]
    if not matching_row_df_2.empty:
        new_df.at[index, 'Smiles_Compound_2'] = matching_row_df_2.iloc[0]['Smiles']

new_df.to_csv('epoxy.csv',index=False)
