from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd

df = pd.read_csv('diCOOH.csv')

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

new_smiles_list = []


for smiles in df['Smiles']:

    m = Chem.MolFromSmiles(smiles)
    
    patt = Chem.MolFromSmarts('C(=O)O')
    repl1 = Chem.MolFromSmiles('[*]')
    repl2 = Chem.MolFromSmiles('[Cs]')
    rms = AllChem.ReplaceSubstructs(m, patt, repl1, replacementConnectionPoint=0)
    rms[0]
    rms = AllChem.ReplaceSubstructs(rms[0], patt, repl2, replacementConnectionPoint=0)
    modified_smiles = Chem.MolToSmiles(rms[0])
    modified_smiles
    mol1=Chem.MolFromSmiles('C1=CC4=C(C=C1C2=CC3=C(C=C2)N=C([N]3)[Rb])N=C([N]4)[*]')
    mol2=Chem.MolFromSmiles(modified_smiles)
    new_smile = combine2frags(mol1,mol2)
    new_smile
    
    new_smiles_list.append(new_smile)

df['Smiles_Compound_1'] = 'C1=CC(=C(C=C1C2=CC(=C(C=C2)N)N)N)N'

df.rename(columns={'Smiles': 'Smiles_Compound_2'}, inplace=True)

df['Smiles'] = new_smiles_list

df.to_csv('PBI.csv', index=False)
