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

#with open('uspto.templates.json', 'r') as f:
#    data = json.load(f)
data = []
with open('smarts-frequency.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        sample = {}
        rea = line.split('>>')[0]
        pro = line.split('>>')[1]
        times = int(pro.split(': ')[1])
        pro = pro.split(': ')[0]
        #if times < 5 and times > 1:
        if times > 10:
            #break
           # continue
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


print('list length:',len(reactants_lst))
print('list length:',len(products_lst))

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

# reactants = run_imap_mp(single_handle, reactants_lst, 16, True)
# result = []
# times = 1
# for x in reactants:
#    if x != '':
#        unique_tpls[x] += 1
# for x in unique_tpls:
#     if unique_tpls[x] == times:
#        result.append(x)
# print('available length:',len(result))
times = 10
dir = 'test-templates'
pro_file = os.path.join(dir, 'uspto.templates-pro.more{}times.smiles.json'.format(times))
rea_file = os.path.join(dir, 'uspto.templates-rea.more{}times.smiles.json'.format(times))
time_file = os.path.join(dir, 'uspto.templates-times.more{}times.smiles.json'.format(times))
with open(pro_file, 'w') as f:
    json.dump(reactants_lst, f) #indeed products
with open(rea_file, 'w') as f:
    json.dump(products_lst, f) #indeed reactants
with open(time_file, 'w') as f:
    json.dump(times_lst, f) #indeed reactants
