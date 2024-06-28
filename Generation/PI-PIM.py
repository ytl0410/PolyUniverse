from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import DataStructs
import pandas as pd

df = pd.read_csv('Polyimide.csv')

def calculate_morgan_fingerprint(smiles, radius=3, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
    return fp

def calculate_similarity(fp1, fp2):
    if fp1 is None or fp2 is None:
        return None
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity

def find_n_connected_to_one_c(mol):
    n_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    c_n_atoms = []
    for n_atom in n_atoms:
        neighbors = mol.GetAtomWithIdx(n_atom).GetNeighbors()
        c_count = sum(1 for neighbor in neighbors if neighbor.GetSymbol() == 'C')
        if c_count == 1:
            c_n_atoms.append(n_atom)
    return c_n_atoms

def replace_first_atom_with_Rb(i,mol):
    if mol.GetNumAtoms() > 0:
        new_mol = Chem.RWMol(mol)
        new_mol.ReplaceAtom(i, Chem.Atom(37)) 
        return new_mol.GetMol()
    else:
        return None


def replace_first_atom_with_star(i,mol):
    if mol.GetNumAtoms() > 0:
        new_mol = Chem.RWMol(mol)
        new_mol.ReplaceAtom(i, Chem.Atom('*')) 
        return new_mol.GetMol()
    else:
        return None


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

def replace_star_with_Sr(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSymbol() == '*':
                    neighbor.SetAtomicNum(38) 
                    break
    new_smiles = Chem.MolToSmiles(mol)
    return new_smiles

for index, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['mon2'])
    m = Chem.MolFromSmiles(row['polym'])

    mol1 = replace_star_with_Sr(row['polym'])
    mol1=Chem.MolFromSmiles(mol1)
    c_n_atoms = find_n_connected_to_one_c(mol)
    max_similarity = 0
    best_new_mol = None
    new_mol =[]

    for i in c_n_atoms:
        new_mol = replace_first_atom_with_Rb(i, mol)

        for j in c_n_atoms:
            if j != i:
                new_moll = replace_first_atom_with_star(j, new_mol)
                new_moll = combine2frags(mol1, new_moll, 'Rb', 'Sr')
                similarity = calculate_similarity(calculate_morgan_fingerprint(row['polym']), calculate_morgan_fingerprint(new_moll))
                if similarity is not None and similarity > max_similarity:
                    max_similarity = similarity
                    best_new_mol = new_moll

    df.loc[index, 'New'] = best_new_mol

for index, row in df.iterrows():
    new_smiles = row['New']
    first_star_index = new_smiles.find('*')
    if first_star_index != -1:  
        new_smiles = new_smiles[:first_star_index] + '[Rb]' + new_smiles[first_star_index + 1:]
    df.at[index, 'New'] = new_smiles
new_smiles_list = []

for smiles in df['New']:
    # 创建一个Molecule对象
    modified_smiles = smiles
    mol2=Chem.MolFromSmiles('C1=CC(=CC=C1C=NCCCCCN=CC2=CC=C(C=C2)C=N[Cs])C=N[*]')
    mol1=Chem.MolFromSmiles(modified_smiles)
    new_smile = combine2frags(mol1,mol2)
    new_smile
    
    # 将新的Smiles添加到列表中
    new_smiles_list.append(new_smile)

# 将新的Smiles列表添加到DataFrame中
df['PI-PIM'] = new_smiles_list
df.drop(columns=['New'], inplace=True)

df.to_csv('PI-PIM.csv', index=False)