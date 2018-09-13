#!/usr/bin/env python3

import sys, argparse, gzip, os, re, subprocess, logging

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomicFileHandler.genomic_file_handlers as genome
import vcfModifier.copy_TextFile as copy_TextFile
import somaticseq.combine_callers as combineCallers


ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger = logging.getLogger('SomaticSeq')
logger.setLevel(logging.DEBUG)
logger.addHandler(ch)


adaTrainer   = os.sep.join( (PRE_DIR, 'r_scripts', 'ada_model_builder_ntChange.R') )
adaPredictor = os.sep.join( (PRE_DIR, 'r_scripts', 'ada_model_predictor.R') )


def runPaired(outdir, ref, tbam, nbam, tumor_name='TUMOR', normal_name='NORMAL', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, inclusion=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, platypus=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', keep_intermediates=False):

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
    outSnv, outIndel, intermediateVcfs, tempFiles = combineCallers.combinePaired(outdir=outdir, ref=ref, tbam=tbam, nbam=nbam, inclusion=inclusion, exclusion=exclusion, mutect=mutect, indelocator=indelocator, mutect2=mutect2, varscan_snv=varscan_snv, varscan_indel=varscan_indel, jsm=jsm, sniper=sniper, vardict=vardict, muse=muse, lofreq_snv=lofreq_snv, lofreq_indel=lofreq_indel, scalpel=scalpel, strelka_snv=strelka_snv, strelka_indel=strelka_indel, tnscope=tnscope, keep_intermediates=True)

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

        subprocess.call( (adaPredictor, classifier_snv, ensembleSnv, classifiedSnvTsv) )

        tsv2vcf.tsv2vcf(classifiedSnvTsv, classifiedSnvVcf, snvCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name=normal_name, tumor_sample_name=tumor_name, print_reject=True, phred_scaled=True)


    else:
        # Train SNV classifier:
        if somaticseq_train and truth_snv:
            subprocess.call( (adaTrainer, ensembleSnv, 'Consistent_Mates', 'Inconsistent_Mates') )

        consensusSnvVcf = os.sep.join(( outdir, consensusOutPrefix + 'sSNV.vcf' ))
        tsv2vcf.tsv2vcf(ensembleSnv, consensusSnvVcf, snvCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name=normal_name, tumor_sample_name=tumor_name, print_reject=True)



    ###################### INDEL ######################
    mutect_infile = intermediateVcfs['MuTect2']['indel'] if intermediateVcfs['MuTect2']['indel'] else indelocator

    somatic_vcf2tsv.vcf2tsv(is_vcf=outIndel, nbam_fn=nbam, tbam_fn=tbam, truth=truth_indel, cosmic=cosmic, dbsnp=dbsnp, mutect=mutect_infile, varscan=varscan_indel, vardict=intermediateVcfs['VarDict']['indel'], lofreq=lofreq_indel, scalpel=scalpel, strelka=strelka_indel, tnscope=intermediateVcfs['TNscope']['indel'], platypus=intermediateVcfs['Platypus']['indel'], dedup=True, min_mq=min_mq, min_bq=min_bq, min_caller=min_caller, ref_fa=ref, p_scale=None, outfile=ensembleIndel)


    # Classify INDEL calls
    if classifier_indel:
        classifiedIndelTsv = os.sep.join(( outdir, classifiedOutPrefix + 'sINDEL.tsv' ))
        classifiedIndelVcf = os.sep.join(( outdir, classifiedOutPrefix + 'sINDEL.vcf' ))

        subprocess.call( (adaPredictor, classifier_indel, ensembleIndel, classifiedIndelTsv) )

        tsv2vcf.tsv2vcf(classifiedIndelTsv, classifiedIndelVcf, indelCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name=normal_name, tumor_sample_name=tumor_name, print_reject=True, phred_scaled=True)

    else:
        # Train INDEL classifier:
        if somaticseq_train and truth_indel:
            subprocess.call( (adaTrainer, ensembleIndel, 'Strelka_QSS', 'Strelka_TQSS', 'Consistent_Mates', 'Inconsistent_Mates') )

        consensusIndelVcf = os.sep.join(( outdir, consensusOutPrefix + 'sINDEL.vcf' ))
        tsv2vcf.tsv2vcf(ensembleIndel, consensusIndelVcf, indelCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name=normal_name, tumor_sample_name=tumor_name, print_reject=True)


    ## Clean up after yourself ##
    if not keep_intermediates:
        for file_i in files_to_delete:
            os.remove(file_i)
            logger.info('Removed {}'.format( file_i ) )






def runSingle(outdir, ref, bam, sample_name='TUMOR', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, inclusion=None, exclusion=None, mutect=None, mutect2=None, varscan=None, vardict=None, lofreq=None, scalpel=None, strelka=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', keep_intermediates=False):

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

        subprocess.call( (adaPredictor, classifier_snv, ensembleSnv, classifiedSnvTsv) )

        tsv2vcf.tsv2vcf(classifiedSnvTsv, classifiedSnvVcf, snvCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=True, paired_mode=False, tumor_sample_name=sample_name, print_reject=True, phred_scaled=True)


    else:
        # Train SNV classifier:
        if somaticseq_train and truth_snv:
            subprocess.call( (adaTrainer, ensembleSnv, 'Consistent_Mates', 'Inconsistent_Mates') )

        consensusSnvVcf = os.sep.join(( outdir, consensusOutPrefix + 'sSNV.vcf' ))
        tsv2vcf.tsv2vcf(ensembleSnv, consensusSnvVcf, snvCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=True, paired_mode=False, tumor_sample_name=sample_name, print_reject=True)



    ###################### INDEL ######################
    single_sample_vcf2tsv.vcf2tsv(is_vcf=outIndel, bam_fn=bam, truth=truth_indel, cosmic=cosmic, dbsnp=dbsnp, mutect=intermediateVcfs['MuTect2']['indel'], varscan=intermediateVcfs['VarScan2']['indel'], vardict=intermediateVcfs['VarDict']['indel'], lofreq=intermediateVcfs['LoFreq']['indel'], scalpel=scalpel, strelka=intermediateVcfs['Strelka']['indel'], dedup=True, min_mq=min_mq, min_bq=min_bq, min_caller=min_caller, ref_fa=ref, p_scale=None, outfile=ensembleIndel)


    # Classify INDEL calls
    if classifier_indel:
        classifiedIndelTsv = os.sep.join(( outdir, classifiedOutPrefix + 'sINDEL.tsv' ))
        classifiedIndelVcf = os.sep.join(( outdir, classifiedOutPrefix + 'sINDEL.vcf' ))

        subprocess.call( (adaPredictor, classifier_indel, ensembleIndel, classifiedIndelTsv) )

        tsv2vcf.tsv2vcf(classifiedIndelTsv, classifiedIndelVcf, indelCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=True, paired_mode=False, tumor_sample_name=sample_name, print_reject=True, phred_scaled=True)

    else:
        # Train INDEL classifier:
        if somaticseq_train and truth_indel:
            subprocess.call( (adaTrainer, ensembleIndel, 'Strelka_QSS', 'Strelka_TQSS', 'Consistent_Mates', 'Inconsistent_Mates') )

        consensusIndelVcf = os.sep.join(( outdir, consensusOutPrefix + 'sINDEL.vcf' ))
        tsv2vcf.tsv2vcf(ensembleIndel, consensusIndelVcf, indelCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=True, paired_mode=False, tumor_sample_name=sample_name, print_reject=True)


    ## Clean up after yourself ##
    if not keep_intermediates:
        for file_i in files_to_delete:
            os.remove(file_i)
            logger.info('Removed {}'.format( file_i ) )






################################################
def run():

    inputParameters = {}

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-outdir',      '--output-directory',   type=str, help='output directory', default='.')
    parser.add_argument('-ref',         '--genome-reference',   type=str, help='.fasta.fai file to get the contigs', required=True)

    parser.add_argument('--truth-snv',         type=str, help='VCF of true hits')
    parser.add_argument('--truth-indel',       type=str, help='VCF of true hits')
    parser.add_argument('--classifier-snv',    type=str, help='RData for SNV')
    parser.add_argument('--classifier-indel',  type=str, help='RData for INDEL')
    parser.add_argument('--pass-threshold',    type=float, help='SCORE for PASS', default=0.5)
    parser.add_argument('--lowqual-threshold', type=float, help='SCORE for LowQual', default=0.1)
    parser.add_argument('-hom', '--homozygous-threshold', type=float, help='VAF for homozygous', default=0.85)
    parser.add_argument('-het', '--heterozygous-threshold', type=float, help='VAF for heterozygous', default=0.01)

    parser.add_argument('-minMQ',     '--minimum-mapping-quality',type=float, help='Minimum mapping quality below which is considered poor', default=1)
    parser.add_argument('-minBQ',     '--minimum-base-quality',   type=float, help='Minimum base quality below which is considered poor', default=5)
    parser.add_argument('-mincaller', '--minimum-num-callers',    type=float, help='Minimum number of tools to be considered', default=0.5)

    parser.add_argument('-dbsnp',  '--dbsnp-vcf',          type=str,   help='dbSNP VCF',)
    parser.add_argument('-cosmic', '--cosmic-vcf',         type=str,   help='COSMIC VCF')

    parser.add_argument('-include',  '--inclusion-region', type=str,   help='inclusion bed')
    parser.add_argument('-exclude',  '--exclusion-region', type=str,   help='exclusion bed')

    parser.add_argument('-nt', '--threads',  type=int, help='number of threads', default=1)

    parser.add_argument('--keep-intermediates',         action='store_true', help='Keep intermediate files', default=False)
    parser.add_argument('-train', '--somaticseq-train', action='store_true', help='Invoke training mode with ground truths', default=False)


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
    inputParameters = vars(args)

    logger.info( 'SomaticSeq Input Arguments: ' + ', '.join( [ '{}={}'.format(i, inputParameters[i])  for i in inputParameters ] ) )

    return inputParameters




################################################################################################
# Execute:
if __name__ == '__main__':
    runParameters = run()

    os.makedirs(runParameters['output_directory'], exist_ok=True)

    if runParameters['which'] == 'paired':

        runPaired( outdir             = runParameters['output_directory'], \
                   ref                = runParameters['genome_reference'], \
                   tbam               = runParameters['tumor_bam_file'], \
                   nbam               = runParameters['normal_bam_file'], \
                   tumor_name         = runParameters['tumor_sample'], \
                   normal_name        = runParameters['normal_sample'], \
                   truth_snv          = runParameters['truth_snv'], \
                   truth_indel        = runParameters['truth_indel'], \
                   classifier_snv     = runParameters['classifier_snv'], \
                   classifier_indel   = runParameters['classifier_indel'], \
                   pass_threshold     = runParameters['pass_threshold'], \
                   lowqual_threshold  = runParameters['lowqual_threshold'], \
                   hom_threshold      = runParameters['homozygous_threshold'], \
                   het_threshold      = runParameters['heterozygous_threshold'], \
                   min_mq             = runParameters['minimum_mapping_quality'], \
                   min_bq             = runParameters['minimum_base_quality'], \
                   min_caller         = runParameters['minimum_num_callers'], \
                   dbsnp              = runParameters['dbsnp_vcf'], \
                   cosmic             = runParameters['cosmic_vcf'], \
                   inclusion          = runParameters['inclusion_region'], \
                   exclusion          = runParameters['exclusion_region'], \
                   mutect             = runParameters['mutect_vcf'], \
                   indelocator        = runParameters['indelocator_vcf'], \
                   mutect2            = runParameters['mutect2_vcf'], \
                   varscan_snv        = runParameters['varscan_snv'], \
                   varscan_indel      = runParameters['varscan_indel'], \
                   jsm                = runParameters['jsm_vcf'], \
                   sniper             = runParameters['somaticsniper_vcf'], \
                   vardict            = runParameters['vardict_vcf'], \
                   muse               = runParameters['muse_vcf'], \
                   lofreq_snv         = runParameters['lofreq_snv'], \
                   lofreq_indel       = runParameters['lofreq_indel'], \
                   scalpel            = runParameters['scalpel_vcf'], \
                   strelka_snv        = runParameters['strelka_snv'], \
                   strelka_indel      = runParameters['strelka_indel'], \
                   tnscope            = runParameters['tnscope_vcf'], \
                   platypus           = runParameters['platypus_vcf'], \
                   somaticseq_train   = runParameters['somaticseq_train'], \
                   keep_intermediates = runParameters['keep_intermediates'] )

    elif runParameters['which'] == 'single':

        runSingle( outdir             = runParameters['output_directory'], \
                   ref                = runParameters['genome_reference'], \
                   bam                = runParameters['bam_file'], \
                   sample_name        = runParameters['sample_name'], \
                   truth_snv          = runParameters['truth_snv'], \
                   truth_indel        = runParameters['truth_indel'], \
                   classifier_snv     = runParameters['classifier_snv'], \
                   classifier_indel   = runParameters['classifier_indel'], \
                   pass_threshold     = runParameters['pass_threshold'], \
                   lowqual_threshold  = runParameters['lowqual_threshold'], \
                   hom_threshold      = runParameters['homozygous_threshold'], \
                   het_threshold      = runParameters['heterozygous_threshold'], \
                   min_mq             = runParameters['minimum_mapping_quality'], \
                   min_bq             = runParameters['minimum_base_quality'], \
                   min_caller         = runParameters['minimum_num_callers'], \
                   dbsnp              = runParameters['dbsnp_vcf'], \
                   cosmic             = runParameters['cosmic_vcf'], \
                   inclusion          = runParameters['inclusion_region'], \
                   exclusion          = runParameters['exclusion_region'], \
                   mutect             = runParameters['mutect_vcf'], \
                   mutect2            = runParameters['mutect2_vcf'], \
                   varscan            = runParameters['varscan_vcf'], \
                   vardict            = runParameters['vardict_vcf'], \
                   lofreq             = runParameters['lofreq_vcf'], \
                   scalpel            = runParameters['scalpel_vcf'], \
                   strelka            = runParameters['strelka_vcf'], \
                   somaticseq_train   = runParameters['somaticseq_train'], \
                   keep_intermediates = runParameters['keep_intermediates'] )
