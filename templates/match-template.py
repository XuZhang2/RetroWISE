import json
from rdkit import Chem
from multi_tqdm import run_imap_mp
from tqdm import tqdm
import argparse
from rdkit import RDLogger
from rdkit.Chem import rdChemReactions
from functools import partial
import multiprocessing

RDLogger.DisableLog('rdApp.*')
parser = argparse.ArgumentParser()
parser.add_argument('-screen_file', type=str)
parser.add_argument('-pro_file', type=str)
parser.add_argument('-out_screen_file', type=str)
parser.add_argument('-out_pro_file', type=str)
parser.add_argument('-out_screen_file_r', type=str)
parser.add_argument('-out_pro_file_r', type=str)
parser.add_argument('-rea_template', type=str)
parser.add_argument('-pro_template', type=str)
opt = parser.parse_args()
screen_file = opt.screen_file
batch_size = int(5e5)

def Smiles2Mol(reactant):
    if reactant == '':
        return None
    smile = reactant.strip().replace(' ', '')
    if len(smile) > 250:
        return None
    try:
        mol = Chem.MolFromSmiles(smile)
        Chem.SanitizeMol(mol)
    except:
        return None
    return mol

def Smarts2Mol(reactant):
    smile = reactant.strip().replace(' ', '')
    try:
        mol = Chem.MolFromSmarts(smile)
        Chem.SanitizeMol(mol)
    except:
        return None
    return mol

def takeSecond(elem):
    return elem[1]

def tplmatch(mol, tpls):
    r_tpl = tpls[0]
    p_tpl = tpls[1]

    mol_reactant = mol[0]
    mol_product = mol[1]

    x = mol_reactant.HasSubstructMatch(r_tpl, useChirality=True)
    y = mol_product.HasSubstructMatch(p_tpl, useChirality=True)

    returned = x and y

    return returned

def structmatch(mol):
    f = partial(tplmatch, mol)
    multiprocessing.current_process().daemon = False
    i_s = run_imap_mp(f, list(zip(rea_templates, pro_templates)), 32, False)
    multiprocessing.current_process().daemon = True
    i_s = [i for i, x in enumerate(i_s) if x]

    for i in i_s:
        tpl = rea_templates_s[i]+'>>'+pro_templates_s[i]
        tpl_m = rdChemReactions.ReactionFromSmarts(tpl)
        reacts = Chem.MolToSmiles(mol[0]).split('.')
        reacts = tuple([Chem.MolFromSmiles(i) for i in reacts])
            #tpl_m.RunReactantInPlace(mol_reactant)
        try:
            pros = tpl_m.RunReactants(reacts)
        except:
            continue
        if pros:
            pros = pros[0]
        if len(pros) != 1:
            continue
        if (Chem.MolToSmiles(mol[1], canonical=True) == Chem.MolToSmiles(pros[0], canonical=True)):
            returned = i+1, True
    return 0, False


def old_structmatch(mol):
    mol_reactant = mol[0]
    mol_product = mol[1]
    for i, (r_tpl, p_tpl) in enumerate(zip(rea_templates, pro_templates)):
        try:
            x = mol_reactant.HasSubstructMatch(r_tpl, useChirality=True)
            y = mol_product.HasSubstructMatch(p_tpl, useChirality=True)
            if x and y:
                tpl = rea_templates_s[i]+'>>'+pro_templates_s[i]
                tpl_m = rdChemReactions.ReactionFromSmarts(tpl)
                reacts = Chem.MolToSmiles(mol_reactant).split('.')
                reacts = tuple([Chem.MolFromSmiles(i) for i in reacts])
                try:
                    pros = tpl_m.RunReactants(reacts)
                except:
                    continue
                if pros:
                    pros = pros[0]
                if len(pros) != 1:
                    continue
                if (Chem.MolToSmiles(mol_product, canonical=True) == Chem.MolToSmiles(pros[0], canonical=True)):
                    returned = i+1, True
                    return returned
        except:
            break
    return 0, False

class AverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

global rea_templates
global pro_templates
global rea_templates_s
global pro_templates_s
with open(opt.rea_template, 'r') as f:
    rea_templates_s = json.load(f)
with open(opt.pro_template, 'r') as f:
    pro_templates_s = json.load(f)

rea_templates = run_imap_mp(Smarts2Mol, rea_templates_s, 16, True)
pro_templates = run_imap_mp(Smarts2Mol, pro_templates_s, 16, True)

none_list = [i for i,_ in enumerate(rea_templates) if _ is None]
for del_index in sorted(none_list, reverse=True):
        del rea_templates[del_index]
        del pro_templates[del_index]
        del rea_templates_s[del_index]
        del pro_templates_s[del_index]
none_list = [i for i,_ in enumerate(pro_templates) if _ is None]
for del_index in sorted(none_list, reverse=True):
        del pro_templates[del_index]
        del rea_templates[del_index]
        del rea_templates_s[del_index]
        del pro_templates_s[del_index]

with open(opt.screen_file, 'r') as f:
    reactants_all = f.readlines()

with open(opt.pro_file, 'r') as f:
    products_all = f.readlines()

count_all = AverageMeter()

with open(opt.out_pro_file, 'w') as fp, open(opt.out_screen_file, 'w') as fs, open(opt.out_pro_file_r, 'w') as fpr, open(opt.out_screen_file_r, 'w') as fsr:    # Use context manager here
    for i in tqdm(range(0, len(reactants_all), batch_size)):
        reactants = reactants_all[i:i + batch_size]
        products = products_all[i:i+batch_size]
        mol_reactants = run_imap_mp(Smiles2Mol, reactants, 16, True)
        mol_products = run_imap_mp(Smiles2Mol, products, 16, True)
        none_list = [i for i,_ in enumerate(mol_reactants) if _ is None]
        none_list2 = [i for i,_ in enumerate(mol_products) if _ is None]
        none_list.extend(none_list2)
        none_list = list(set(none_list))
        for del_index in sorted(none_list, reverse=True):
            del mol_reactants[del_index]
            del mol_products[del_index]
            del reactants[del_index]
            del products[del_index]
        print('#Satinized Molecule:{}'.format(len(mol_reactants)))
        matched_reactions = run_imap_mp(old_structmatch, list(zip(mol_reactants, mol_products)), 62, True)
        #matched_reactions = []
        #for x,y in tqdm(list(zip(mol_reactants, mol_products))):
        #    matched_reaction = old_structmatch((x,y))
        #    matched_reactions.append(matched_reaction)
        count = AverageMeter()
        for c,j in tqdm(enumerate(matched_reactions)):
            if not j[1]:
                print(reactants[c], file=fs, end='')
                print(products[c], file=fp, end='')
            else:
                count.update(j[0])
                count_all.update(j[0])
                print(reactants[c], file=fsr, end='')
                print(products[c], file=fpr, end='')
        print('selected #reactants:{} in batch {}, average template used to match:{}'.format(count.count, i // batch_size, count.avg))

    print('selected #reactants:{} in all batch, average template used to match:{}'.format(count_all.count, count_all.avg))

