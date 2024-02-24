import numpy as np
import pandas as pd
import argparse
import os
import re
import random
import textdistance
import multiprocessing

from rdkit import Chem
from tqdm import tqdm


from rdkit import RDLogger
from multi_tqdm import run_imap_mp
RDLogger.DisableLog('rdApp.*')

def takeSecond(elem):
    return elem[1]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-source_dir', type=str, help='A existing tgt file in PtoR')
    parser.add_argument('-process', type=int, default=8)
    parser.add_argument('-augment', type=int, default=5)

    opt = parser.parse_args()

    dirs = ['train', 'val', 'test']
    for cur_dir in dirs:
        source_file = os.path.join(opt.source_dir, cur_dir, 'tgt-{}.txt'.format(cur_dir))
        target_dir = 'dataset/RtoR/{}/{}-rand_tgt'.format(os.path.basename(opt.source_dir), cur_dir)
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        src_file = os.path.join(target_dir, 'src-{}.txt'.format(cur_dir))
        tgt_file = os.path.join(target_dir, 'tgt-{}.txt'.format(cur_dir))
        with open(source_file, 'r') as file:

            lines = file.readlines()

            largest_reactants = []
            res_reactants = []
            for smi in tqdm(lines):
                sub_smi = smi.split(".")
                sub_mol = [Chem.MolFromSmiles(smiles.strip().replace(' ', ''), isomericSmiles=True) for smiles in sub_smi]
                sub_mol_size = [(sub_smi[i].strip(), len(m.GetAtoms())) for i, m in enumerate(sub_mol) if m is not None]
                sub_mol_size.sort(key=takeSecond, reverse=True) # the first is the largest

                if len(sub_mol_size) <= 1:
                    continue
                cur_res_reactants = []
                for i, smi in enumerate(sub_mol_size):
                    smi = smi[0]
                    if i == 0:
                        cur_largest_reactants = smi
                    else:
                        cur_res_reactants.append(smi)

                cur_largest_reactants += '\n'
                cur_res_reactants = '.'.join(cur_res_reactants)# + '\n'
                cur_res_mols = Chem.MolFromSmiles(cur_res_reactants.strip().replace(' ', ''), isomericSmiles=True)
                cur_res_reactants = Chem.MolToSmiles(cur_res_mols, isomericSmiles=True, rootedAtAtom=(i%opt.augment), canonical=True)

                largest_reactants.append(cur_largest_reactants)
                res_reactants.append(cur_res_reactants)

            with open(src_file, 'w') as f:
                f.writelines(largest_reactants)
            with open(tgt_file, 'w') as f:
                f.writelines(res_reactants)




