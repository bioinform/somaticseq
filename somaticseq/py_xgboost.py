#!/usr/bin/env python3

import argparse
import xgboost as xgb
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import logging
from copy import copy
from somaticseq._version import  __version__


DEFAULT_PARAM = {'max_depth': 12, 'nthread': 1, 'seed': 0, 'objective': 'binary:logistic'}
NON_FEATURE   = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Strelka_QSS', 'Strelka_TQSS', 'if_COSMIC', 'COSMIC_CNT', 'TrueVariant_or_False']




def builder(input_tsv, param=DEFAULT_PARAM, non_feature=NON_FEATURE, num_rounds=200, model=None):

    if not model:
        model = input_tsv + '.xgb.v{}.model'.format( __version__ )
    
    data        = pd.read_csv(input_tsv, sep='\t')
    
    train_data  = data.drop(non_feature, axis=1)
    train_label = data['TrueVariant_or_False']
    
    dtrain      = xgb.DMatrix(train_data, label=train_label)
    bst         = xgb.train(param, dtrain, num_boost_round=num_rounds)
    
    bst.save_model(model)
    bst.dump_model(model+'.txt')

    return model




def predictor(model, input_tsv, output_tsv, non_feature=NON_FEATURE):
    
    input_data = pd.read_csv(input_tsv, sep='\t')
    test_data  = input_data.drop(non_feature, axis=1)
    dtest      = xgb.DMatrix(test_data)
    
    xgb_model = xgb.Booster()
    xgb_model.load_model(model)
    
    scores = model.predict(dtest)
    input_data.assign(SCORE = scores)
    
    input_data.to_csv(output_tsv, sep='\t')
    
    return output_tsv






################################################################################################
# Execute:
if __name__ == '__main__':
    
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger = logging.getLogger('pyXGBOOST')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(ch)


    parser = argparse.ArgumentParser(description="Run XGBoost", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    sample_parsers = parser.add_subparsers(title="mode")

    # TRAINING mode
    parser_paired = sample_parsers.add_parser('train')
    parser_paired.add_argument('-tsv',     '--tsv-in',           type=str, help='labeled tsv file',  required=True)
    parser_paired.add_argument('-out',     '--model-out',        type=str, help='output model file name')
    parser_paired.add_argument('-iter',    '--num-boost-rounds', type=int, default=200)
    parser_paired.add_argument('-threads', '--num-threads',      type=int, help='num threads')
    parser_paired.add_argument('-depth',   '--max-depth',        type=int, help='tree max depth')
    parser_paired.add_argument('-seed',    '--seed',             type=int, help='random seed')
    parser_paired.set_defaults(which='train')


    # PREDICTION mode
    parser_single = sample_parsers.add_parser('predict')
    parser_single.add_argument('-model', '--model',         type=str, help='xgboost model',  required=True)
    parser_single.add_argument('-tsv',   '--tsv-in',        type=str, help='tsv file in',    required=True)
    parser_single.add_argument('-out',   '--predicted-tsv', type=str, help='tsv file out',   required=True)
    parser_single.set_defaults(which='predict')

    args = parser.parse_args()



    if args.which == 'train':
        
        PARAM = DEFAULT_PARAM
        
        if args.num_threads:
            PARAM['nthread'] = args.num_threads
            
        if args.max_depth:
            PARAM['max_depth'] = args.max_depth
            
        if args.seed:
            PARAM['seed'] = args.seed
        
        builder(args.tsv_in, param=PARAM, non_feature=NON_FEATURE, num_rounds=args.num_boost_rounds, model=args.model_out):



    elif args.which == 'predict':
        predictor(args.model, args.tsv_in, args.predicted_tsv, non_feature=NON_FEATURE)
