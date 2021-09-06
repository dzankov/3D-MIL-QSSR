import sys
sys.path.append('/home/zankov/dev/miqsar')

import os
import pickle
import joblib
import pkg_resources
import numpy as np
import pandas as pd
from itertools import groupby
from sklearn.pipeline import Pipeline
from CGRtools import RDFRead, RDFWrite
from CIMtools.preprocessing import Fragmentor, CGR, EquationTransformer, SolventVectorizer
from miqssr.descriptor_calculation.pmapper_3d import calc_pmapper_descriptors
from miqssr.conformer_generation.gen_conformers import gen_confs

fragmentor_path = pkg_resources.resource_filename(__name__, '.')
os.environ['PATH'] += ':{}'.format(fragmentor_path)

def read_pkl(fname):
    with open(fname, 'rb') as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break

def calc_3d_pmapper(conf_files, path='.', ncpu=10):
    for conf in conf_files:
        dsc_file = calc_pmapper_descriptors(conf, path=path, ncpu=ncpu, col_clean=None, del_undef=True)

        with open(dsc_file, 'rb') as inp:
            data = joblib.load(inp)

        if 'mol_title' not in data.columns:
            data = data.reset_index()
        data['mol_id'] = data['mol_id'].str.lower()

        out_fname = dsc_file.replace('_proc.pkl', '.csv')
        data.to_csv(out_fname, index=False)

    return out_fname

def calc_2d_isida(fname, path='.'):
    reacts = RDFRead(fname, remap=False).read()
    for reaction in reacts:
        reaction.standardize()
        reaction.kekule()
        reaction.implicify_hydrogens()
        reaction.thiele()
        reaction.clean2d()

    frag = Pipeline(
        [('CGR', CGR()), ('frg', Fragmentor(fragment_type=9, max_length=5, useformalcharge=True, version='2017.x'))])
    res = frag.fit_transform(reacts)
    res['react_id'] = [i.meta['ID'] for i in reacts]
    #
    out_fname = os.path.join(path, '2DDescrISIDA_cgr-data_0.csv')
    res.to_csv(out_fname, index=False)

    import shutil
    del frag
    frg_files = [i for i in os.listdir() if i.startswith('frg')]
    for file in frg_files:
        if os.path.isfile(file):
            os.remove(file)
        else:
            shutil.rmtree(file)

    return out_fname

def create_catalyst_input_file(input_fname=None):
    data = RDFRead(input_fname, remap=False).read()

    groups = []
    for k, g in groupby(data, lambda x: x.meta['CATALYST_SMILES']):
        groups.append(list(g))
    #
    smiles = [i[0].meta['CATALYST_SMILES'] for i in groups]
    act = [np.mean([float(i.meta['SELECTIVITY']) for i in x]) for x in groups]
    #
    res = []
    ids = {}
    for i in range(len(smiles)):
        res.append({'SMILES': smiles[i], 'MOL_ID': i, 'ACT': act[i]})
        ids[i] = [x.meta['ID'] for x in groups[i]]
    #
    out_fname = 'catalyst_data.smi'
    res = pd.DataFrame(res)
    res.to_csv(out_fname, index=False, header=False)

    return out_fname, ids


def calc_descriptors(input_fname=None, nconfs=5, energy=10, ncpu=5, path='.'):
    cat_data_file, cat_ids = create_catalyst_input_file(input_fname)

    conf_files = gen_confs(cat_data_file, ncpu=ncpu, nconfs_list=[nconfs], stereo=False, energy=energy, path=path)
    os.remove(cat_data_file)

    react_out_fname = calc_2d_isida(input_fname, path=path)
    cat_out_fname = calc_3d_pmapper(conf_files, ncpu=ncpu, path=path)

    #
    reacts = pd.read_csv(react_out_fname, index_col='react_id')
    catalysts = pd.read_csv(cat_out_fname, index_col='mol_id').sort_index()
    #
    res = []
    for i in catalysts.index.unique():
        for j in cat_ids[i]:
            cat = catalysts.loc[i:i]
            react = pd.concat([reacts.loc[j:j]] * len(cat))
            react_cat = pd.concat([react, cat.set_index(react.index)], axis=1)

            res.append(react_cat)

    out_fname = os.path.join(path, 'PhFprPmapper_concat-data_{}.csv'.format(nconfs))
    pd.concat(res).to_csv(out_fname)

    return out_fname