import os
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

from sklearn.preprocessing import MaxAbsScaler, MinMaxScaler

from pmapper.customize import load_smarts

from miqssr.conformer_generation.gen_conformers import gen_confs
from miqssr.descriptor_calculation.pmapper_3d import calc_pmapper_descriptors


def rdkit_numpy_convert(fp):
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)


def scale_data(x_train, x_test):
    scaler = MinMaxScaler()
    scaler.fit(np.vstack(x_train))
    x_train_scaled = x_train.copy()
    x_test_scaled = x_test.copy()
    for i, bag in enumerate(x_train):
        x_train_scaled[i] = scaler.transform(bag)
    for i, bag in enumerate(x_test):
        x_test_scaled[i] = scaler.transform(bag)
    return np.array(x_train_scaled), np.array(x_test_scaled)


def load_svm_data(fname):
    def str_to_vec(dsc_str, dsc_num):

        tmp = {}
        for i in dsc_str.split(' '):
            tmp[int(i.split(':')[0])] = int(i.split(':')[1])
        #
        tmp_sorted = {}
        for i in range(dsc_num + 1):
            tmp_sorted[i] = tmp.get(i, 0)
        vec = list(tmp_sorted.values())

        return vec

    #
    with open(fname) as f:
        dsc_tmp = [i.strip() for i in f.readlines()]

    with open(fname.replace('txt', 'rownames')) as f:
        idx_tmp = [int(i.split('_')[0]) for i in f.readlines()]
    #
    dsc_num = max([max([int(j.split(':')[0]) for j in i.strip().split(' ')]) for i in dsc_tmp])
    #
    bags, idx = [], []
    for cat_idx in list(np.unique(idx_tmp)):
        bag, idx_ = [], []
        for dsc_str, i in zip(dsc_tmp, idx_tmp):
            if i == cat_idx:
                bag.append(str_to_vec(dsc_str, dsc_num))
                idx_.append(i)

        bags.append(np.array(bag).astype('uint8'))
        idx.append(idx_[0])
    #
    bags, idx = np.array(bags, dtype=object), np.array(idx, dtype=object)
    bags_dict = {k: v for k, v in zip(idx, bags)}

    return bags_dict


def read_input_data(fname, cols=None):
    data = pd.read_csv(fname)

    reactants = data[cols['REACTANTS']]
    catalysts = data[cols['CATALYST']]

    try:
        activity = np.array(list(data[cols['ACTIVITY']]))
    except:
        activity = [None for i in reactants.values]

    # indexes
    cat_idx = {i: n for n, i in enumerate(catalysts.unique())}

    react_cat_idx = {}
    for n, cat in zip(range(len(reactants)), catalysts):
        react_cat_idx[n] = cat_idx[cat]

    return reactants, catalysts, activity, react_cat_idx


def gen_catalyst_confs(mols, idx, num_confs=[50], energy=50, num_cpu=1, path='.'):
    #
    res = []
    for smi, i in zip(mols, idx):
        res.append({'SMILES': smi, 'MOL_ID': idx[i], 'ACT': None})
    #
    cat_file = os.path.join(path, 'catalyst_data.smi')
    res = pd.DataFrame(res).drop_duplicates()
    res.to_csv(cat_file, index=False, header=False)
    #
    conf_files = gen_confs(cat_file,
                           nconfs_list=[num_confs],
                           stereo=False,
                           energy=energy,
                           ncpu=num_cpu,
                           path=path,
                           verbose=False)

    return conf_files


def calc_reaction_descr(reacts, path='.'):

    # train descr
    descr = []
    for react in reacts.values:
        x = np.concatenate(rdkit_numpy_convert([AllChem.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(mol), 2) for mol in react]))
        descr.append(x)
    descr = pd.DataFrame(np.array(descr))
    #
    descr.to_csv(os.path.join(path, 'reaction_descr.csv'), index=False)

    return descr


def calc_catalyst_descr(conf_file=None, smarts_features=None, num_descr=[3], ncpu=5, path='.'):
    smarts_features = load_smarts(smarts_features)

    conf_file = os.path.join(conf_file, 'catalyst_conformers.pkl')
    out_fname = os.path.join(path, 'catalyst_descriptors.txt')

    calc_pmapper_descriptors(inp_fname=conf_file,
                             out_fname=out_fname,
                             smarts_features=smarts_features,
                             descr_num=num_descr,
                             remove=0.05,
                             keep_temp=False,
                             ncpu=ncpu,
                             verbose=False)

    # read catalyst descr file
    bags_dict = load_svm_data(out_fname)

    return bags_dict


def concat_react_cat_descr(reacts_descr, bags_dict, react_cat_idx):
    bags = []
    for react_i, cat_i in react_cat_idx.items():
        cat_bag = bags_dict[cat_i]
        react_vec = reacts_descr.iloc[[react_i]].values
        react_bag = np.repeat(react_vec, len(cat_bag), axis=0)
        #
        bags.append(np.concatenate((react_bag, cat_bag), axis=1))

    return np.array(bags, dtype=object)

