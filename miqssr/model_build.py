import os
import yaml
import logging
import pandas as pd

from miqssr.estimators.wrappers import BagWrapperMLPRegressor, InstanceWrapperMLPRegressor
from miqssr.estimators.attention_nets import AttentionNetRegressor
from miqssr.estimators.mi_nets import BagNetRegressor, InstanceNetRegressor

from miqssr.utils import (read_input_data, gen_catalyst_confs, calc_catalyst_descr,
                          calc_reaction_descr, concat_react_cat_descr, scale_data)


logging.basicConfig(format='%(message)s', level=logging.NOTSET)

def build_model(config_path):

    # read config file
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    # create results folder
    results_path = config['General']['results_path']
    if not os.path.exists(results_path):
        os.mkdir(results_path)
        os.mkdir(os.path.join(results_path, 'train'))
        os.mkdir(os.path.join(results_path, 'test'))

    # read input data
    cols = {'REACTANTS': config['Data']['reactants'],
            'CATALYST': config['Data']['catalyst'][0],
            'ACTIVITY': config['Data']['activity'][0]}

    reacts_train, cat_train, activity_train, react_cat_idx_train = read_input_data(
        config['Data']['training_dataset_path'], cols=cols)

    reacts_test, cat_test, activity_test, react_cat_idx_test = read_input_data(
        config['Data']['test_dataset_path'], cols=cols)
    #
    logging.info(f'Read input data: {len(reacts_train)} training reactions, {len(reacts_test)} test reactions')

    # generate catalyst conformers
    confs_train_path = os.path.join(results_path, 'train', 'conformers')
    confs_test_path = os.path.join(results_path, 'test', 'conformers')

    if not os.path.exists(confs_train_path):
        os.mkdir(confs_train_path)

        logging.info(f'Generate conformers for training catalysts ...')

        gen_catalyst_confs(cat_train, react_cat_idx_train,
                           num_confs=config['Conformers']['max_num_conf'],
                           energy=config['Conformers']['energy_window'],
                           num_cpu=config['General']['num_cpu'],
                           path=confs_train_path)

    if not os.path.exists(confs_test_path):
        os.mkdir(confs_test_path)

        logging.info(f'Generate conformers for test catalysts ...')

        gen_catalyst_confs(cat_test, react_cat_idx_test,
                           num_confs=config['Conformers']['max_num_conf'],
                           energy=config['Conformers']['energy_window'],
                           num_cpu=config['General']['num_cpu'],
                           path=confs_test_path)

    # calculate reaction descriptors
    reacts_train_path = os.path.join(results_path, 'train', 'descriptors')
    reacts_test_path = os.path.join(results_path, 'test', 'descriptors')

    if not os.path.exists(reacts_train_path):
        os.mkdir(reacts_train_path)

    logging.info(f'Calculate descriptors for training reactions ...')

    reacts_descr_train = calc_reaction_descr(reacts_train, path=reacts_train_path)

    if not os.path.exists(reacts_test_path):
        os.mkdir(reacts_test_path)

    logging.info(f'Calculate descriptors for test reactions ...')

    reacts_descr_test = calc_reaction_descr(reacts_test, path=reacts_test_path)

    # calculate catalyst descriptors
    logging.info(f'Calculate descriptors for training catalysts ...')

    bags_dict_train = calc_catalyst_descr(conf_file=confs_train_path,
                                          smarts_features=config['Descriptors']['smarts_path'],
                                          ncpu=config['General']['num_cpu'],
                                          num_descr=[config['Descriptors']['num_descr']],
                                          path=reacts_train_path)

    logging.info(f'Calculate descriptors for test catalysts ...')

    bags_dict_test = calc_catalyst_descr(conf_file=confs_test_path,
                                         smarts_features=config['Descriptors']['smarts_path'],
                                         ncpu=config['General']['num_cpu'],
                                         num_descr=[config['Descriptors']['num_descr']],
                                         path=reacts_test_path)

    # concatenate reaction and catalyst descriptors
    train_cols = pd.read_csv(os.path.join(reacts_train_path, 'catalyst_descriptors.colnames'), header=None)
    test_cols = pd.read_csv(os.path.join(reacts_test_path, 'catalyst_descriptors.colnames'), header=None)

    common_train_cols = []
    common_test_cols = []
    for n, dsc in enumerate(train_cols[0]):
        if dsc in list(test_cols[0]):
            common_train_cols.append(n)
            common_test_cols.append((test_cols == dsc)[0].argmax())
            #
    for k, v in bags_dict_train.items():
        bags_dict_train[k] = v[:, common_train_cols]
    for k, v in bags_dict_test.items():
        bags_dict_test[k] = v[:, common_test_cols]
    #
    x_train = concat_react_cat_descr(reacts_descr_train, bags_dict_train, react_cat_idx_train)
    x_test = concat_react_cat_descr(reacts_descr_test, bags_dict_test, react_cat_idx_test)
    x_train_scaled, x_test_scaled = scale_data(x_train, x_test)

    # train model
    logging.info(f'Model training ...')

    algorithms = {'BagWrapperMLPRegressor': BagWrapperMLPRegressor,
                  'InstanceWrapperMLPRegressor': InstanceWrapperMLPRegressor,
                  'BagNetRegressor': BagNetRegressor,
                  'InstanceNetRegressor': InstanceNetRegressor,
                  'AttentionNetRegressor': AttentionNetRegressor}
    #
    init_cuda = True if config['General']['num_gpu'] > 0 else False
    #
    n_dim = [x_test[0].shape[-1]] + [256, 128, 64]
    net = algorithms[config['Algorithm']['type']](ndim=n_dim,
                                                  pool=config['Algorithm']['pooling'],
                                                  init_cuda=init_cuda)
    #
    net.fit(x_train_scaled, activity_train,
            n_epoch=config['Algorithm']['n_epoch'],
            lr=config['Algorithm']['learning_rate'],
            weight_decay=config['Algorithm']['weight_decay'],
            batch_size=9999999,
            verbose=False)

    # save predictions
    logging.info(f'Save test set predictions')

    data_test = pd.read_csv(config['Data']['test_dataset_path'])
    data_test['ACTIVITY_PREDICTION'] = net.predict(x_test_scaled)
    #
    data_test.to_csv(os.path.join(config['General']['results_path'], 'test_set_predictions.csv'), index=False)

    logging.info(f'Finished!')

