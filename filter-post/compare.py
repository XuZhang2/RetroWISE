from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import os
from rdkit import RDLogger
from tqdm import tqdm
from tqdm.contrib import tzip
import argparse

parser = argparse.ArgumentParser('')
parser.add_argument('-thresh', type=float)

args = parser.parse_args()

RDLogger.DisableLog('rdApp.*')

global basedir
basedir = '../dataset/USPTO_extra_RtoP_aug5/train/filter-5/'

def filter_files(file1, file2, threshold, file3):
    filter_dir = os.path.join(basedir, str(threshold))
    if not os.path.isdir(filter_dir):
        os.makedirs(filter_dir)
    with open(file1, 'r') as f1, open(file2, 'r') as f2, open(file3, 'r') as f3:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        lines3 = f3.readlines()
        print('origin length:{}'.format(len(lines3)))
        filtered_lines1 = []
        filtered_lines3 = []
        for line1, line2, line3 in tzip(lines1, lines2, lines3):
            try:
                mol1 = Chem.MolFromSmiles(line1.strip().replace(' ', ''))
                mol2 = Chem.MolFromSmiles(line2.strip().replace(' ', ''))
                fp1 = AllChem.GetMorganFingerprint(mol1, 2)
                fp2 = AllChem.GetMorganFingerprint(mol2, 2)
            except:
                continue

            similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            if similarity >= threshold:
                filtered_lines1.append(line1)
                filtered_lines3.append(line3)
        screen_file = os.path.join(filter_dir, 'tgt-train.txt')
        company_file = os.path.join(filter_dir, 'src-train.txt')
        print('filter length:{}'.format(len(filtered_lines3)))
        with open(screen_file, 'w') as f5:
            f5.writelines(filtered_lines1)
        with open(company_file, 'w') as f4:
            f4.writelines(filtered_lines3)
        print('done')

file1 = os.path.join(basedir, 'tgt-train.txt')
file2 = os.path.join(basedir, 'tgt-train-backT.txt')
file3 = os.path.join(basedir, 'src-train.txt')
threshold = args.thresh

filter_files(file1, file2, threshold, file3)
