import json
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from multi_tqdm import run_imap_mp
from tqdm import tqdm

import argparse
from rdkit import RDLogger
from collections import Counter
import os

RDLogger.DisableLog('rdApp.*')

data = []
with open('smarts-frequency.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        sample = {}
        rea = line.split('>>')[0]
        pro = line.split('>>')[1]
        times = int(pro.split(': ')[1])
        pro = pro.split(': ')[0]
        if times > 5:
            sample['reactants'] = rea
            sample['products'] = pro
            sample['times'] = times
            data.append(sample)


print('data length:', len(data))
reactants_lst = []
products_lst = []
times_lst = []

for i in data:
    try:
        x = i['reactants']
        y = i['products']
        z = i['times']
    except:
        continue
    reactants_lst.append(x)
    products_lst.append(y)
    times_lst.append(z)


print('list length:', len(reactants_lst))
print('list length:', len(products_lst))

def smarts_to_smiles(smarts):
    mol = Chem.MolFromSmarts(smarts)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    smiles = Chem.MolToSmiles(mol)
    return smiles

def single_handle(data):
    x = data
    mol = Chem.MolFromSmarts(x)
    try:
        Chem.SanitizeMol(mol)
    except:
        return ''
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    smiles = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smiles)
    try:
        Chem.SanitizeMol(mol)
    except:
        return ''
    return smiles

times = 5
dir = 'test-templates'
pro_file = os.path.join(dir, 'uspto.templates-pro.more{}times.smiles.json'.format(times))
rea_file = os.path.join(dir, 'uspto.templates-rea.more{}times.smiles.json'.format(times))
time_file = os.path.join(dir, 'uspto.templates-times.more{}times.smiles.json'.format(times))
with open(pro_file, 'w') as f:
    json.dump(reactants_lst, f)
with open(rea_file, 'w') as f:
    json.dump(products_lst, f) 
with open(time_file, 'w') as f:
    json.dump(times_lst, f)