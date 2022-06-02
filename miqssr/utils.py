import os
import pickle
import joblib
import pkg_resources
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import groupby
from sklearn.pipeline import Pipeline
from CGRtools import RDFRead, RDFWrite
from CIMtools.preprocessing import Fragmentor, CGR, EquationTransformer, SolventVectorizer
from pmapper.customize import load_smarts, load_factory
from .conformer_generation.gen_conformers import gen_confs
from .descriptor_calculation.pmapper_3d import calc_pmapper_descriptors
from .descriptor_calculation.rdkit_morgan import calc_morgan_descriptors

fragmentor_path = pkg_resources.resource_filename(__name__, './')
os.environ['PATH'] += ':{}'.format(fragmentor_path)


def read_pkl(fname):
    with open(fname, 'rb') as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break

                
def calc_3d_rdkit(conf_file, path='.', ncpu=10):
    out_fname = calc_rdkit_descriptors(conf_file, path=path, ncpu=ncpu)
    return out_fname


def calc_3d_pmapper(conf_file, path='.', ncpu=10, verbose=True):
    
    out_fname = os.path.join(path, 'PhFprPmapper_{}.txt'.format(os.path.basename(conf_file).split('.')[0]))
    

    smarts_features = load_smarts('./miqssr/smarts_features.txt')
    factory = load_factory('./miqssr/smarts_features.fdef')
    
    calc_pmapper_descriptors(inp_fname=conf_file, out_fname=out_fname, 
                             smarts_features=smarts_features, factory=factory,
                             descr_num=[3], remove=0.05, keep_temp=False, ncpu=ncpu, verbose=verbose)

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
    res['act'] = [i.meta['SELECTIVITY'] for i in reacts]
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


def calc_2d_descr(input_fname=None, ncpu=10, path='.'):
    cat_data_file, cat_ids = create_catalyst_input_file(input_fname)
    react_out_fname = calc_2d_isida(input_fname, path=path)
    
    out_fnames = []
    for dsc_calc_func, dsc_name in zip([calc_2d_descriptors, calc_morgan_descriptors], ['2DDescrRDKit', 'MorganFprRDKit']):

        cat_out_fname = dsc_calc_func(cat_data_file, ncpu=ncpu, path=path)
        #
        reacts = pd.read_csv(react_out_fname, index_col='react_id')
        catalysts = pd.read_csv(cat_out_fname, index_col='mol_id').sort_index().drop(['act'], axis=1)
        #
        res = []
        for i in catalysts.index.unique():
            for j in cat_ids[i]:
                cat = catalysts.loc[i:i]
                react = pd.concat([reacts.loc[j:j]] * len(cat))
                react_cat = pd.concat([react, cat.set_index(react.index)], axis=1)

                res.append(react_cat)

        out_fname = os.path.join(path, f'{dsc_name}_concat-data_0.csv')
        pd.concat(res).to_csv(out_fname)
        out_fnames.append(out_fname)
    # 
    os.remove(cat_data_file)

    return out_fnames


def calc_3d_descr(input_fname=None, nconfs=[1, 50], energy=10, ncpu=5, path='.', verbose=True):
    cat_data_file, cat_ids = create_catalyst_input_file(input_fname)
    
    conf_files = gen_confs(cat_data_file, ncpu=ncpu, nconfs_list=nconfs, stereo=False, energy=energy, path=path, verbose=verbose)

    os.remove(cat_data_file)
    
    react_out_fname = calc_2d_isida(input_fname, path=path)
    
    out_fnames = []
    for nconf, conf_file in zip(nconfs, conf_files):

        for dsc_calc_func, dsc_name in zip([calc_3d_pmapper], ['PhFprPmapper']):

            cat_out_fname = dsc_calc_func(conf_file, ncpu=ncpu, path=path, verbose=verbose)
            raws_out_fname = cat_out_fname.replace('txt', 'rownames')

            with open(cat_out_fname) as f:
                cat_dsc = f.readlines()

            with open(raws_out_fname) as f:
                cat_names = [i.strip() for i in f.readlines()]

            cats_dict = defaultdict(list)
            for i, j in zip(cat_names, cat_dsc):
                cats_dict[int(i.split('_')[0])].append(j)
            #
            reacts = pd.read_csv(react_out_fname, index_col='react_id')
            labels = reacts['act']

            reacts_dict = {}
            for i in reacts.drop(['act'], axis=1).index:

                tmp = []
                for n, d in enumerate(reacts.loc[i]):
                    if d != 0:
                        tmp.append(f'{n}:{d:.0f}')
                reacts_dict[i] = tmp
        #
        shift = len(reacts.columns) + 1
        react_rownames = []

        out_fname_dsc = os.path.join(path, f'{dsc_name}_concat_data_{nconf}.txt')
        with open(out_fname_dsc, 'w') as f:
            for cat_name, cat_dsc in cats_dict.items():
                for i in cat_ids[cat_name]:
                    react = reacts_dict[i]
                    for cat_conf in cats_dict[cat_name]:
                        cat_conf = cat_conf.split(' ')
                        cat_conf = [':'.join([str(int(i.split(':')[0]) + shift), i.split(':')[1].strip()]) for i in  cat_conf]

                        react_cat_str = ' '.join(react + cat_conf) + '\n'

                        f.write(react_cat_str)

                        react_rownames.append(f'{i}:{labels[i]}' + '\n')
        #
        out_fname_rows = os.path.join(path, f'{dsc_name}_concat_data_{nconf}.rownames')
        with open(out_fname_rows, 'w') as f:
            f.write(''.join(react_rownames))

    return out_fnames