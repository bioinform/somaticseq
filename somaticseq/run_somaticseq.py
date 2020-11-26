#!/usr/bin/env python3

import argparse, os, subprocess, logging
import somaticseq.combine_callers as combineCallers
from somaticseq._version import  __version__

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'

#ch = logging.StreamHandler()
#ch.setLevel(logging.DEBUG)
#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
#ch.setFormatter(formatter)

logger = logging.getLogger('SomaticSeq')
logger.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.INFO, format=FORMAT)

#logger.addHandler(ch)

DEFAULT_XGB_BOOST_ROUNDS  = 500
DEFAULT_NUM_TREES_PREDICT = 100



def modelTrainer(input_file, algo, threads=1, seed=0, max_depth=12, iterations=200, features_to_exclude=[]):
    
    logger = logging.getLogger(modelTrainer.__name__)
    
    if algo == 'ada' or algo == 'ada.R':
        command_item = ('ada_model_builder_ntChange.R', input_file)
        logger.info(' '.join(command_item))
        exit_code = subprocess.call( command_item )
        assert exit_code == 0

        return input_file + '.ada.Classifier.RData'

    if algo == 'xgboost':
        
        import somaticseq.somatic_xgboost as somatic_xgboost
        
        xgb_param = somatic_xgboost.DEFAULT_PARAM
        xgb_param['nthread']   = threads
        xgb_param['max_depth'] = max_depth
        xgb_param['seed']      = seed
        
        non_features = somatic_xgboost.NON_FEATURE
        for feature_i in features_to_exclude:
            non_features.append( feature_i )
        
        logger.info( 'PARAMETER: ' + ', '.join( [ '{}: {}'.format(i, xgb_param[i]) for i in xgb_param ] ) )
        
        xgb_model = somatic_xgboost.builder( [input_file,], param=xgb_param, non_feature=non_features, num_rounds=iterations )
        
        return xgb_model


    
def modelPredictor(input_file, output_file, algo, classifier, iterations=100, features_to_exclude=[]):
    
    logger = logging.getLogger(modelPredictor.__name__)
    
    if algo == 'ada' or algo == 'ada.R':
        command_item = ('ada_model_predictor.R', classifier, input_file, output_file)
        logger.info(' '.join(command_item))
        exit_code = subprocess.call( command_item )
        assert exit_code == 0

        return output_file


    if algo == 'xgboost':
        import somaticseq.somatic_xgboost as somatic_xgboost
        
        non_features = somatic_xgboost.NON_FEATURE
        for feature_i in features_to_exclude:
            non_features.append( feature_i )

        somatic_xgboost.predictor(classifier, input_file, output_file, non_features, iterations)
        
        return output_file



def runPaired(outdir, ref, tbam, nbam, tumor_name='TUMOR', normal_name='NORMAL', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, inclusion=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, platypus=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', algo='ada', keep_intermediates=False, train_seed=0, tree_depth=12, iterations=None, features_excluded=[]):

    logger = logging.getLogger(runPaired.__name__)

    import somaticseq.somatic_vcf2tsv as somatic_vcf2tsv
    import somaticseq.SSeq_tsv2vcf as tsv2vcf

    files_to_delete = set()

    snvCallers = []
    if mutect or mutect2: snvCallers.append('MuTect')
    if varscan_snv:       snvCallers.append('VarScan2')
    if jsm:               snvCallers.append('JointSNVMix2')
    if sniper:            snvCallers.append('SomaticSniper')
    if vardict:           snvCallers.append('VarDict')
    if muse:              snvCallers.append('MuSE')
    if lofreq_snv:        snvCallers.append('LoFreq')
    if strelka_snv:       snvCallers.append('Strelka')
    if tnscope:           snvCallers.append('TNscope')
    if platypus:          snvCallers.append('Platypus')


    indelCallers = []
    if indelocator or mutect2: indelCallers.append('MuTect')
    if varscan_indel:          indelCallers.append('VarScan2')
    if vardict:                indelCallers.append('VarDict')
    if lofreq_indel:           indelCallers.append('LoFreq')
    if scalpel:                indelCallers.append('Scalpel')
    if strelka_indel:          indelCallers.append('Strelka')
    if tnscope:                indelCallers.append('TNscope')
    if platypus:               indelCallers.append('Platypus')

    # Function to combine individual VCFs into a simple VCF list of variants:
    outSnv, outIndel, intermediateVcfs, tempFiles = combineCallers.combinePaired(outdir=outdir, ref=ref, tbam=tbam, nbam=nbam, inclusion=inclusion, exclusion=exclusion, mutect=mutect, indelocator=indelocator, mutect2=mutect2, varscan_snv=varscan_snv, varscan_indel=varscan_indel, jsm=jsm, sniper=sniper, vardict=vardict, muse=muse, lofreq_snv=lofreq_snv, lofreq_indel=lofreq_indel, scalpel=scalpel, strelka_snv=strelka_snv, strelka_indel=strelka_indel, tnscope=tnscope, platypus=platypus, keep_intermediates=True)

    files_to_delete.add(outSnv)
    files_to_delete.add(outIndel)
    [ files_to_delete.add(i) for i in tempFiles ]


    ensembleSnv   = os.sep.join(( outdir, ensembleOutPrefix + 'sSNV.tsv' ))
    ensembleIndel = os.sep.join(( outdir, ensembleOutPrefix + 'sINDEL.tsv' ))


    ######################  SNV  ######################
    mutect_infile = intermediateVcfs['MuTect2']['snv'] if intermediateVcfs['MuTect2']['snv'] else mutect

    somatic_vcf2tsv.vcf2tsv(is_vcf=outSnv, nbam_fn=nbam, tbam_fn=tbam, truth=truth_snv, cosmic=cosmic, dbsnp=dbsnp, mutect=mutect_infile, varscan=varscan_snv, jsm=jsm, sniper=sniper, vardict=intermediateVcfs['VarDict']['snv'], muse=muse, lofreq=lofreq_snv, scalpel=None, strelka=strelka_snv, tnscope=intermediateVcfs['TNscope']['snv'], platypus=intermediateVcfs['Platypus']['snv'], dedup=True, min_mq=min_mq, min_bq=min_bq, min_caller=min_caller, ref_fa=ref, p_scale=None, outfile=ensembleSnv)


    # Classify SNV calls
    if classifier_snv:
        classifiedSnvTsv = os.sep.join(( outdir, classifiedOutPrefix + 'sSNV.tsv' ))
        classifiedSnvVcf = os.sep.join(( outdir, classifiedOutPrefix + 'sSNV.vcf' ))

        iterations = iterations if iterations else DEFAULT_NUM_TREES_PREDICT
        modelPredictor(ensembleSnv, classifiedSnvTsv, algo, classifier_snv, iterations=iterations, features_to_exclude=features_excluded)

        extra_header = ['##SomaticSeqClassifier={}'.format(classifier_snv), ]

        tsv2vcf.tsv2vcf(classifiedSnvTsv, classifiedSnvVcf, snvCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name=normal_name, tumor_sample_name=tumor_name, print_reject=True, phred_scaled=True, extra_headers=extra_header)


    else:
        # Train SNV classifier:
        if somaticseq_train and truth_snv:
            
            iterations = iterations if iterations else DEFAULT_XGB_BOOST_ROUNDS
            modelTrainer(ensembleSnv, algo, threads=1, seed=train_seed, max_depth=tree_depth, iterations=iterations, features_to_exclude=features_excluded)
        
        consensusSnvVcf = os.sep.join(( outdir, consensusOutPrefix + 'sSNV.vcf' ))
        tsv2vcf.tsv2vcf(ensembleSnv, consensusSnvVcf, snvCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name=normal_name, tumor_sample_name=tumor_name, print_reject=True)



    ###################### INDEL ######################
    mutect_infile = intermediateVcfs['MuTect2']['indel'] if intermediateVcfs['MuTect2']['indel'] else indelocator

    somatic_vcf2tsv.vcf2tsv(is_vcf=outIndel, nbam_fn=nbam, tbam_fn=tbam, truth=truth_indel, cosmic=cosmic, dbsnp=dbsnp, mutect=mutect_infile, varscan=varscan_indel, vardict=intermediateVcfs['VarDict']['indel'], lofreq=lofreq_indel, scalpel=scalpel, strelka=strelka_indel, tnscope=intermediateVcfs['TNscope']['indel'], platypus=intermediateVcfs['Platypus']['indel'], dedup=True, min_mq=min_mq, min_bq=min_bq, min_caller=min_caller, ref_fa=ref, p_scale=None, outfile=ensembleIndel)


    # Classify INDEL calls
    if classifier_indel:
        classifiedIndelTsv = os.sep.join(( outdir, classifiedOutPrefix + 'sINDEL.tsv' ))
        classifiedIndelVcf = os.sep.join(( outdir, classifiedOutPrefix + 'sINDEL.vcf' ))

        iterations = iterations if iterations else DEFAULT_NUM_TREES_PREDICT
        modelPredictor(ensembleIndel, classifiedIndelTsv, algo, classifier_indel, iterations=iterations, features_to_exclude=features_excluded)
        
        extra_header = ['##SomaticSeqClassifier={}'.format(classifier_indel), ]
        
        tsv2vcf.tsv2vcf(classifiedIndelTsv, classifiedIndelVcf, indelCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name=normal_name, tumor_sample_name=tumor_name, print_reject=True, phred_scaled=True, extra_headers=extra_header)

    else:
        # Train INDEL classifier:
        if somaticseq_train and truth_indel:
            
            iterations = iterations if iterations else DEFAULT_XGB_BOOST_ROUNDS
            modelTrainer(ensembleIndel, algo, threads=1, seed=train_seed, max_depth=tree_depth, iterations=iterations, features_to_exclude=features_excluded)

        consensusIndelVcf = os.sep.join(( outdir, consensusOutPrefix + 'sINDEL.vcf' ))
        tsv2vcf.tsv2vcf(ensembleIndel, consensusIndelVcf, indelCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name=normal_name, tumor_sample_name=tumor_name, print_reject=True)


    ## Clean up after yourself ##
    if not keep_intermediates:
        for file_i in files_to_delete:
            os.remove(file_i)
            logger.info('Removed {}'.format( file_i ) )






def runSingle(outdir, ref, bam, sample_name='TUMOR', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, inclusion=None, exclusion=None, mutect=None, mutect2=None, varscan=None, vardict=None, lofreq=None, scalpel=None, strelka=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', algo='ada', keep_intermediates=False, train_seed=0, tree_depth=12, iterations=None, features_excluded=[]):

    logger = logging.getLogger(runSingle.__name__)

    import somaticseq.single_sample_vcf2tsv as single_sample_vcf2tsv
    import somaticseq.SSeq_tsv2vcf as tsv2vcf

    files_to_delete = set()

    snvCallers = []
    if mutect or mutect2: snvCallers.append('MuTect')
    if varscan:           snvCallers.append('VarScan2')
    if vardict:           snvCallers.append('VarDict')
    if lofreq:            snvCallers.append('LoFreq')
    if strelka:           snvCallers.append('Strelka')


    indelCallers = []
    if mutect2: indelCallers.append('MuTect2')
    if varscan: indelCallers.append('VarScan2')
    if vardict: indelCallers.append('VarDict')
    if lofreq:  indelCallers.append('LoFreq')
    if scalpel: indelCallers.append('Scalpel')
    if strelka: indelCallers.append('Strelka')


    # Function to combine individual VCFs into a simple VCF list of variants:
    outSnv, outIndel, intermediateVcfs, tempFiles = combineCallers.combineSingle(outdir=outdir, ref=ref, bam=bam, inclusion=inclusion, exclusion=exclusion, mutect=mutect, mutect2=mutect2, varscan=varscan, vardict=vardict, lofreq=lofreq, scalpel=scalpel, strelka=strelka, keep_intermediates=True)

    files_to_delete.add(outSnv)
    files_to_delete.add(outIndel)
    [ files_to_delete.add(i) for i in tempFiles ]

    ensembleSnv   = os.sep.join(( outdir, ensembleOutPrefix + 'sSNV.tsv' ))
    ensembleIndel = os.sep.join(( outdir, ensembleOutPrefix + 'sINDEL.tsv' ))


    ######################  SNV  ######################
    mutect_infile = intermediateVcfs['MuTect2']['snv'] if intermediateVcfs['MuTect2']['snv'] else mutect

    single_sample_vcf2tsv.vcf2tsv(is_vcf=outSnv, bam_fn=bam, truth=truth_snv, cosmic=cosmic, dbsnp=dbsnp, mutect=mutect_infile, varscan=intermediateVcfs['VarScan2']['snv'], vardict=intermediateVcfs['VarDict']['snv'], lofreq=intermediateVcfs['LoFreq']['snv'], scalpel=None, strelka=intermediateVcfs['Strelka']['snv'], dedup=True, min_mq=min_mq, min_bq=min_bq, min_caller=min_caller, ref_fa=ref, p_scale=None, outfile=ensembleSnv)


    # Classify SNV calls
    if classifier_snv:
        classifiedSnvTsv = os.sep.join(( outdir, classifiedOutPrefix + 'sSNV.tsv' ))
        classifiedSnvVcf = os.sep.join(( outdir, classifiedOutPrefix + 'sSNV.vcf' ))

        iterations = iterations if iterations else DEFAULT_NUM_TREES_PREDICT
        modelPredictor(ensembleSnv, classifiedSnvTsv, algo, classifier_snv, iterations=iterations, features_to_exclude=features_excluded)
        
        extra_header = ['##SomaticSeqClassifier={}'.format(classifier_snv), ]
        
        tsv2vcf.tsv2vcf(classifiedSnvTsv, classifiedSnvVcf, snvCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=True, paired_mode=False, tumor_sample_name=sample_name, print_reject=True, phred_scaled=True, extra_headers=extra_header)


    else:
        # Train SNV classifier:
        if somaticseq_train and truth_snv:
            
            iterations = iterations if iterations else DEFAULT_XGB_BOOST_ROUNDS
            modelTrainer(ensembleSnv, algo, threads=1, seed=train_seed, max_depth=tree_depth, iterations=iterations, features_to_exclude=features_excluded)

        consensusSnvVcf = os.sep.join(( outdir, consensusOutPrefix + 'sSNV.vcf' ))
        tsv2vcf.tsv2vcf(ensembleSnv, consensusSnvVcf, snvCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=True, paired_mode=False, tumor_sample_name=sample_name, print_reject=True)



    ###################### INDEL ######################
    single_sample_vcf2tsv.vcf2tsv(is_vcf=outIndel, bam_fn=bam, truth=truth_indel, cosmic=cosmic, dbsnp=dbsnp, mutect=intermediateVcfs['MuTect2']['indel'], varscan=intermediateVcfs['VarScan2']['indel'], vardict=intermediateVcfs['VarDict']['indel'], lofreq=intermediateVcfs['LoFreq']['indel'], scalpel=scalpel, strelka=intermediateVcfs['Strelka']['indel'], dedup=True, min_mq=min_mq, min_bq=min_bq, min_caller=min_caller, ref_fa=ref, p_scale=None, outfile=ensembleIndel)


    # Classify INDEL calls
    if classifier_indel:
        classifiedIndelTsv = os.sep.join(( outdir, classifiedOutPrefix + 'sINDEL.tsv' ))
        classifiedIndelVcf = os.sep.join(( outdir, classifiedOutPrefix + 'sINDEL.vcf' ))

        iterations = iterations if iterations else DEFAULT_NUM_TREES_PREDICT
        modelPredictor(ensembleIndel, classifiedIndelTsv, algo, classifier_indel, iterations=iterations, features_to_exclude=features_excluded)
        
        extra_header = ['##SomaticSeqClassifier={}'.format(classifier_indel), ]
        
        tsv2vcf.tsv2vcf(classifiedIndelTsv, classifiedIndelVcf, indelCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=True, paired_mode=False, tumor_sample_name=sample_name, print_reject=True, phred_scaled=True, extra_headers=extra_header)

    else:
        # Train INDEL classifier:
        if somaticseq_train and truth_indel:
            
            iterations = iterations if iterations else DEFAULT_XGB_BOOST_ROUNDS
            modelTrainer(ensembleIndel, algo, threads=1, seed=train_seed, max_depth=tree_depth, iterations=iterations, features_to_exclude=features_excluded)

        consensusIndelVcf = os.sep.join(( outdir, consensusOutPrefix + 'sINDEL.vcf' ))
        tsv2vcf.tsv2vcf(ensembleIndel, consensusIndelVcf, indelCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=True, paired_mode=False, tumor_sample_name=sample_name, print_reject=True)


    ## Clean up after yourself ##
    if not keep_intermediates:
        for file_i in files_to_delete:
            os.remove(file_i)
            logger.info('Removed {}'.format( file_i ) )






################################################
def run():

    parser = argparse.ArgumentParser(description='Somaticseq v{}: a method to combine results from multiple somatic mutation callers, extract genomic and sequencing features for each variant call from the BAM files, and then use machine learning to score the variants.'.format(__version__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-outdir',      '--output-directory',   type=str, help='output directory', default='.')
    parser.add_argument('-ref',         '--genome-reference',   type=str, help='.fasta.fai file to get the contigs', required=True)

    parser.add_argument('--truth-snv',         type=str, help='VCF of true hits')
    parser.add_argument('--truth-indel',       type=str, help='VCF of true hits')
    parser.add_argument('--classifier-snv',    type=str, help='RData for SNV')
    parser.add_argument('--classifier-indel',  type=str, help='RData for INDEL')
    parser.add_argument('--pass-threshold',    type=float, help='SCORE for PASS', default=0.5)
    parser.add_argument('--lowqual-threshold', type=float, help='SCORE for LowQual', default=0.1)

    parser.add_argument('-algo', '--algorithm', type=str, help='ada or xgboost', default='xgboost', choices=('ada', 'xgboost', 'ada.R'))
    parser.add_argument('-hom',  '--homozygous-threshold', type=float, help='VAF for homozygous', default=0.85)
    parser.add_argument('-het',  '--heterozygous-threshold', type=float, help='VAF for heterozygous', default=0.01)

    parser.add_argument('-minMQ',     '--minimum-mapping-quality',type=float, help='Minimum mapping quality below which is considered poor', default=1)
    parser.add_argument('-minBQ',     '--minimum-base-quality',   type=float, help='Minimum base quality below which is considered poor', default=5)
    parser.add_argument('-mincaller', '--minimum-num-callers',    type=float, help='Minimum number of tools to be considered', default=0.5)

    parser.add_argument('-dbsnp',  '--dbsnp-vcf',          type=str,   help='dbSNP VCF',)
    parser.add_argument('-cosmic', '--cosmic-vcf',         type=str,   help='COSMIC VCF')

    parser.add_argument('-include',  '--inclusion-region', type=str,   help='inclusion bed')
    parser.add_argument('-exclude',  '--exclusion-region', type=str,   help='exclusion bed')

    parser.add_argument('-nt', '--threads',  type=int, help='number of threads', default=1)

    parser.add_argument('-train',  '--somaticseq-train', action='store_true', help='Invoke training mode with ground truths', default=False)
    parser.add_argument('-seed',   '--seed',        type=int, help='seed for xgboost training', default=0)
    parser.add_argument('-tdepth', '--tree-depth',  type=int, help='max tree depth for xgboost training', default=12)
    parser.add_argument('-iters',  '--iterations',  type=int, help='num boosting rounds for xgboost: default is 500 for training and 100 for predicting, i.e., by default, 500 trees are built for classifier, but only the first 100 trees are used.')
    parser.add_argument('--features-excluded',      type=str, nargs='*', help='features to exclude for xgboost training. Must be same for train/predict.', default=[] )

    parser.add_argument('--keep-intermediates', action='store_true', help='Keep intermediate files', default=False)

    # Modes:
    sample_parsers = parser.add_subparsers(title="sample_mode")

    # Paired Sample mode
    parser_paired = sample_parsers.add_parser('paired')
    parser_paired.add_argument('-tbam',          '--tumor-bam-file',    type=str,   help='Tumor BAM File',  required=True)
    parser_paired.add_argument('-nbam',          '--normal-bam-file',   type=str,   help='Normal BAM File', required=True)

    parser_paired.add_argument('-tumorSM',       '--tumor-sample',      type=str,   help='Tumor Name',  default='TUMOR')
    parser_paired.add_argument('-normalSM',      '--normal-sample',     type=str,   help='Normal Name', default='NORMAL')

    parser_paired.add_argument('-mutect',        '--mutect-vcf',        type=str,   help='MuTect VCF',        )
    parser_paired.add_argument('-indelocator',   '--indelocator-vcf',   type=str,   help='Indelocator VCF',   )
    parser_paired.add_argument('-mutect2',       '--mutect2-vcf',       type=str,   help='MuTect2 VCF',       )
    parser_paired.add_argument('-varscansnv',    '--varscan-snv',       type=str,   help='VarScan2 VCF',      )
    parser_paired.add_argument('-varscanindel',  '--varscan-indel',     type=str,   help='VarScan2 VCF',      )
    parser_paired.add_argument('-jsm',           '--jsm-vcf',           type=str,   help='JointSNVMix2 VCF',  )
    parser_paired.add_argument('-sniper',        '--somaticsniper-vcf', type=str,   help='SomaticSniper VCF', )
    parser_paired.add_argument('-vardict',       '--vardict-vcf',       type=str,   help='VarDict VCF',       )
    parser_paired.add_argument('-muse',          '--muse-vcf',          type=str,   help='MuSE VCF',          )
    parser_paired.add_argument('-lofreqsnv',     '--lofreq-snv',        type=str,   help='LoFreq VCF',        )
    parser_paired.add_argument('-lofreqindel',   '--lofreq-indel',      type=str,   help='LoFreq VCF',        )
    parser_paired.add_argument('-scalpel',       '--scalpel-vcf',       type=str,   help='Scalpel VCF',       )
    parser_paired.add_argument('-strelkasnv',    '--strelka-snv',       type=str,   help='Strelka VCF',       )
    parser_paired.add_argument('-strelkaindel',  '--strelka-indel',     type=str,   help='Strelka VCF',       )
    parser_paired.add_argument('-tnscope',       '--tnscope-vcf',       type=str,   help='TNscope VCF',       )
    parser_paired.add_argument('-platypus',      '--platypus-vcf',      type=str,   help='Platypus VCF',      )

    parser_paired.set_defaults(which='paired')


    # Single Sample mode
    parser_single = sample_parsers.add_parser('single')
    parser_single.add_argument('-bam',     '--bam-file',    type=str, help='BAM File',     required=True)
    parser_single.add_argument('-SM',      '--sample-name', type=str, help='Sample Name',  default='TUMOR')
    parser_single.add_argument('-mutect',  '--mutect-vcf',  type=str, help='MuTect VCF',   )
    parser_single.add_argument('-mutect2', '--mutect2-vcf', type=str, help='MuTect2 VCF',  )
    parser_single.add_argument('-varscan', '--varscan-vcf', type=str, help='VarScan2 VCF', )
    parser_single.add_argument('-vardict', '--vardict-vcf', type=str, help='VarDict VCF',  )
    parser_single.add_argument('-lofreq',  '--lofreq-vcf',  type=str, help='LoFreq VCF',   )
    parser_single.add_argument('-scalpel', '--scalpel-vcf', type=str, help='Scalpel VCF',  )
    parser_single.add_argument('-strelka', '--strelka-vcf', type=str, help='Strelka VCF',  )

    parser_single.set_defaults(which='single')

    args = parser.parse_args()
    #inputParameters = vars(args)

    logger.info( 'SomaticSeq Input Arguments: ' + ', '.join( [ '{}={}'.format(i, vars(args)[i])  for i in vars(args) ] ) )

    return args




################################################################################################
# Execute:
if __name__ == '__main__':
    
    args = run()

    os.makedirs(args.output_directory, exist_ok=True)

    if args.which == 'paired':

        runPaired( outdir             = args.output_directory, \
                   ref                = args.genome_reference, \
                   tbam               = args.tumor_bam_file, \
                   nbam               = args.normal_bam_file, \
                   tumor_name         = args.tumor_sample, \
                   normal_name        = args.normal_sample, \
                   truth_snv          = args.truth_snv, \
                   truth_indel        = args.truth_indel, \
                   classifier_snv     = args.classifier_snv, \
                   classifier_indel   = args.classifier_indel, \
                   pass_threshold     = args.pass_threshold, \
                   lowqual_threshold  = args.lowqual_threshold, \
                   hom_threshold      = args.homozygous_threshold, \
                   het_threshold      = args.heterozygous_threshold, \
                   min_mq             = args.minimum_mapping_quality, \
                   min_bq             = args.minimum_base_quality, \
                   min_caller         = args.minimum_num_callers, \
                   dbsnp              = args.dbsnp_vcf, \
                   cosmic             = args.cosmic_vcf, \
                   inclusion          = args.inclusion_region, \
                   exclusion          = args.exclusion_region, \
                   mutect             = args.mutect_vcf, \
                   indelocator        = args.indelocator_vcf, \
                   mutect2            = args.mutect2_vcf, \
                   varscan_snv        = args.varscan_snv, \
                   varscan_indel      = args.varscan_indel, \
                   jsm                = args.jsm_vcf, \
                   sniper             = args.somaticsniper_vcf, \
                   vardict            = args.vardict_vcf, \
                   muse               = args.muse_vcf, \
                   lofreq_snv         = args.lofreq_snv, \
                   lofreq_indel       = args.lofreq_indel, \
                   scalpel            = args.scalpel_vcf, \
                   strelka_snv        = args.strelka_snv, \
                   strelka_indel      = args.strelka_indel, \
                   tnscope            = args.tnscope_vcf, \
                   platypus           = args.platypus_vcf, \
                   algo               = args.algorithm, \
                   somaticseq_train   = args.somaticseq_train, \
                   train_seed         = args.seed, \
                   tree_depth         = args.tree_depth, \
                   iterations         = args.iterations, \
                   features_excluded  = args.features_excluded, \
                   keep_intermediates = args.keep_intermediates, \
                   )

    elif args.which == 'single':

        runSingle( outdir             = args.output_directory, \
                   ref                = args.genome_reference, \
                   bam                = args.bam_file, \
                   sample_name        = args.sample_name, \
                   truth_snv          = args.truth_snv, \
                   truth_indel        = args.truth_indel, \
                   classifier_snv     = args.classifier_snv, \
                   classifier_indel   = args.classifier_indel, \
                   pass_threshold     = args.pass_threshold, \
                   lowqual_threshold  = args.lowqual_threshold, \
                   hom_threshold      = args.homozygous_threshold, \
                   het_threshold      = args.heterozygous_threshold, \
                   min_mq             = args.minimum_mapping_quality, \
                   min_bq             = args.minimum_base_quality, \
                   min_caller         = args.minimum_num_callers, \
                   dbsnp              = args.dbsnp_vcf, \
                   cosmic             = args.cosmic_vcf, \
                   inclusion          = args.inclusion_region, \
                   exclusion          = args.exclusion_region, \
                   mutect             = args.mutect_vcf, \
                   mutect2            = args.mutect2_vcf, \
                   varscan            = args.varscan_vcf, \
                   vardict            = args.vardict_vcf, \
                   lofreq             = args.lofreq_vcf, \
                   scalpel            = args.scalpel_vcf, \
                   strelka            = args.strelka_vcf, \
                   algo               = args.algorithm, \
                   somaticseq_train   = args.somaticseq_train, \
                   train_seed         = args.seed, \
                   tree_depth         = args.tree_depth, \
                   iterations         = args.iterations, \
                   features_excluded  = args.features_excluded, \
                   keep_intermediates = args.keep_intermediates, \
                   )
