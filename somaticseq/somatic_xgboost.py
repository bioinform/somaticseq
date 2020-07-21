#!/usr/bin/env python3

import argparse
import xgboost as xgb
import pandas as pd
import numpy as np
import re
import logging
import somaticseq.ntchange_type as ntchange
from copy import copy
from somaticseq._version import  __version__


DEFAULT_PARAM = {'max_depth': 12, 'nthread': 1, 'seed': 0, 'objective': 'binary:logistic'}
NON_FEATURE   = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Strelka_QSS', 'Strelka_TQSS', 'if_COSMIC', 'COSMIC_CNT', 'TrueVariant_or_False']




def builder(input_tsvs, param=DEFAULT_PARAM, non_feature=NON_FEATURE, num_rounds=200, model=None):

    logger = logging.getLogger( 'XGBOOST_' + builder.__name__)
    logger.info('TRAINING {} for XGBOOST'.format( ','.join(input_tsvs)) )
    logger.info('Columns removed before training: {}'.format( ', '.join(non_feature)) )
    logger.info('Number of boosting rounds = {}'.format(num_rounds) )
    logger.info( 'PARAMETER: ' + ', '.join( [ '{}: {}'.format(i, param[i]) for i in param ] ) )


    if not model:
        model = input_tsvs[0] + '.xgb.v{}.classifier'.format( __version__ )
    
    input_data = pd.concat( [pd.read_csv(input_tsv_i, sep='\t', low_memory=False) for input_tsv_i in input_tsvs] )
    
    train_data = ntchange.ntchange(input_data)
    for non_feature_i in non_feature:
        if non_feature_i in train_data:
            train_data.drop([non_feature_i,], axis=1, inplace=True)
    
    train_label   = input_data['TrueVariant_or_False']
    
    dtrain      = xgb.DMatrix(train_data, label=train_label)
    bst         = xgb.train(param, dtrain, num_boost_round=num_rounds)
    
    bst.save_model(model)
    bst.dump_model(model+'.txt')

    return model




def predictor(model, input_tsv, output_tsv, non_feature=NON_FEATURE):

    logger = logging.getLogger( 'XGBOOST_' + predictor.__name__)
    logger.info('Columns removed for prediction: {}'.format( ','.join(non_feature)) )

    xgb_model = xgb.Booster()
    xgb_model.load_model(model)

    # Read in chunks, start with write/Header, and append/NoHeader afterwards.
    chunksize = 10000
    writeMode, writeHeader = 'w', True
    
    for input_data in pd.read_csv(input_tsv, sep='\t', chunksize=chunksize, low_memory=False):
    
        test_data = ntchange.ntchange(input_data)
        for non_feature_i in non_feature:
            if non_feature_i in test_data:
                test_data.drop(non_feature_i, axis=1, inplace=True)
    
        dtest     = xgb.DMatrix(test_data)
        scores    = xgb_model.predict(dtest)
        predicted = input_data.assign(SCORE = scores)
    
        predicted.to_csv(output_tsv, sep='\t', index=False, mode=writeMode, header=writeHeader, na_rep='nan')
        
        writeMode, writeHeader = 'a', False

    return output_tsv






################################################################################################
# Execute:
if __name__ == '__main__':
    
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger = logging.getLogger('SomaticSeq_XGBOOST')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(ch)


    parser = argparse.ArgumentParser(description="Run XGBoost", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    sample_parsers = parser.add_subparsers(title="mode")

    # TRAINING mode
    parser_paired = sample_parsers.add_parser('train')
    parser_paired.add_argument('-tsvs',    '--tsvs-in',          type=str, nargs='+', help='labeled tsv file(s)',  required=True)
    parser_paired.add_argument('-out',     '--model-out',        type=str, help='output model file name')
    parser_paired.add_argument('-iter',    '--num-boost-rounds', type=int, default=200)
    parser_paired.add_argument('-threads', '--num-threads',      type=int, help='num threads')
    parser_paired.add_argument('-depth',   '--max-depth',        type=int, help='tree max depth')
    parser_paired.add_argument('-seed',    '--seed',             type=int, help='random seed')
    parser_paired.add_argument('--features-excluded',            type=str, nargs='*', help='features to exclude for xgboost training. Must be same for train/predict.', default=[] )
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
            
        for feature_i in args.features_excluded:
            NON_FEATURE.append( feature_i )
        
        builder(args.tsvs_in, param=PARAM, non_feature=NON_FEATURE, num_rounds=args.num_boost_rounds, model=args.model_out)



    elif args.which == 'predict':
        predictor(args.model, args.tsv_in, args.predicted_tsv, non_feature=NON_FEATURE)
