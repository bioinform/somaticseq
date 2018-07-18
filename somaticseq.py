#!/usr/bin/env python3

import sys, argparse, gzip, os, re, subprocess
import genomicFileHandler.genomic_file_handlers as genome
import vcfModifier.copy_TextFile as copy_TextFile
import somaticseq.combine_callers as combineCallers


MY_DIR = os.path.dirname(os.path.realpath(__file__))

adaTrainer   = MY_DIR +  os.sep + 'r_scripts' + os.sep + 'ada_model_builder_ntChange.R'
adaPredictor = MY_DIR +  os.sep + 'r_scripts' + os.sep + 'ada_model_predictor.R'



def runPaired(outdir, ref, tbam, nbam, tumor_name='TUMOR', normal_name='NORMAL', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, inclusion=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, keep_intermediates=False):
    
    import somaticseq.somatic_vcf2tsv as somatic_vcf2tsv
    import somaticseq.SSeq_tsv2vcf as tsv2vcf
    
    files_to_delete = set()
    
    snvCallers = []
    if mutect or mutect2:
        snvCallers.append('MuTect')
    if varscan_snv:
        snvCallers.append('VarScan2')
    if jsm:
        snvCallers.append('JointSNVMix2')
    if sniper:
        snvCallers.append('SomaticSniper')
    if vardict:
        snvCallers.append('VarDict')
    if muse:
        snvCallers.append('MuSE')
    if lofreq_snv:
        snvCallers.append('LoFreq')
    if strelka_snv:
        snvCallers.append('Strelka')
    if tnscope:
        snvCallers.append('TNscope')
    
    
    indelCallers = []
    if indelocator or mutect2:
        indelCallers.append('MuTect')
    if varscan_indel:
        indelCallers.append('VarScan2')
    if vardict:
        indelCallers.append('VarDict')
    if lofreq_indel:
        indelCallers.append('LoFreq')
    if scalpel:
        indelCallers.append('Scalpel')
    if strelka_indel:
        indelCallers.append('Strelka')
    if tnscope:
        indelCallers.append('TNscope')
    
    
    
    # Function to combine individual VCFs into a simple VCF list of variants:
    outSnv, outIndel, intermediateVcfs, tempFiles = combineCallers.combinePaired(outdir, ref, tbam, nbam, inclusion, exclusion, mutect, indelocator, mutect2, varscan_snv, varscan_indel, jsm, sniper, vardict, muse, lofreq_snv, lofreq_indel, scalpel, strelka_snv, strelka_indel, tnscope, keep_intermediates=True)
    
    files_to_delete.add(outSnv)
    files_to_delete.add(outIndel)
    [ files_to_delete.add(i) for i in tempFiles ]
    
    ensembleSnv   = outdir + os.sep + 'Ensemble.sSNV.tsv'
    ensembleIndel = outdir + os.sep + 'Ensemble.sINDEL.tsv'
    
    
    
    ######################  SNV  ######################
    somatic_vcf2tsv.vcf2tsv(is_vcf=outSnv, nbam_fn=nbam, tbam_fn=tbam, truth=truth_snv, cosmic=cosmic, dbsnp=dbsnp, mutect=intermediateVcfs['MuTect2']['snv'], varscan=varscan_snv, jsm=jsm, sniper=sniper, vardict=intermediateVcfs['VarDict']['snv'], muse=muse, lofreq=lofreq_snv, scalpel=None, strelka=strelka_snv, tnscope=intermediateVcfs['TNscope']['snv'], ref=ref, dedup=True, min_mq=min_mq, min_bq=min_bq, min_caller=min_caller, ref_fa=ref, p_scale=None, outfile=ensembleSnv)
    
        
    # Classify SNV calls
    if classifier_snv:
        classifiedSnvTsv = outdir + os.sep + 'SSeq.Classified.sSNV.tsv'
        classifiedSnvVcf = outdir + os.sep + 'SSeq.Classified.sSNV.vcf'
        
        subprocess.call( (adaPredictor, classifier_snv, ensembleSnv, classifiedSnvTsv) )
        
        tsv2vcf.tsv2vcf(classifiedSnvTsv, classifiedSnvVcf, snvCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name='NORMAL', tumor_sample_name='TUMOR', print_reject=True, phred_scaled=True)
    
    
    else:
        # Train SNV classifier:
        if somaticseq_train and truth_snv:
            subprocess.call( (adaTrainer, ensembleSnv, 'Consistent_Mates' 'Inconsistent_Mates') )
        
        
        consensusSnvVcf = outdir + os.sep + 'Consensus.sSNV.vcf'
        tsv2vcf.tsv2vcf(ensembleSnv, consensusSnvVcf, snvCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name='NORMAL', tumor_sample_name='TUMOR', print_reject=True)
    
    
    
    ###################### INDEL ######################
    somatic_vcf2tsv.vcf2tsv(is_vcf=outIndel, nbam_fn=nbam, tbam_fn=tbam, truth=truth_indel, cosmic=cosmic, dbsnp=dbsnp, mutect=intermediateVcfs['MuTect2']['indel'], varscan=varscan_indel, vardict=intermediateVcfs['VarDict']['indel'], lofreq=lofreq_indel, scalpel=scalpel, strelka=strelka_indel, tnscope=intermediateVcfs['TNscope']['indel'], ref=ref, dedup=True, min_mq=min_mq, min_bq=min_bq, min_caller=min_caller, ref_fa=ref, p_scale=None, outfile=ensembleIndel)


    # Classify INDEL calls
    if classifier_indel:
        classifiedIndelTsv = outdir + os.sep + 'SSeq.Classified.sINDEL.tsv'
        classifiedIndelVcf = outdir + os.sep + 'SSeq.Classified.sINDEL.vcf'
        
        subprocess.call( (adaPredictor, classifier_indel, ensembleIndel, classifiedIndelTsv) )
        
        tsv2vcf.tsv2vcf(classifiedIndelTsv, classifiedIndelVcf, indelCallers, pass_score=pass_threshold, lowqual_score=lowqual_threshold, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name='NORMAL', tumor_sample_name='TUMOR', print_reject=True, phred_scaled=True)

    else:
        # Train INDEL classifier:
        if somaticseq_train and truth_indel:
            subprocess.call( (adaTrainer, ensembleIndel, 'Strelka_QSS', 'Strelka_TQSS', 'Consistent_Mates', 'Inconsistent_Mates') )
        
        
        consensusIndelVcf = outdir + os.sep + 'Consensus.sINDEL.vcf'
        tsv2vcf.tsv2vcf(ensembleIndel, consensusIndelVcf, indelCallers, hom_threshold=hom_threshold, het_threshold=het_threshold, single_mode=False, paired_mode=True, normal_sample_name='NORMAL', tumor_sample_name='TUMOR', print_reject=True)



    ## Clean up after yourself ##
    if not keep_intermediates:
        for file_i in files_to_delete:
            subprocess.call( ('rm', '-v', file_i ) )




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
    
    parser.add_argument('--keep-intermediates', action='store_true', help='Keep intermediate files', default=False)
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
    parser_paired.add_argument('-strelka-snv',   '--strelka-snv',       type=str,   help='Strelka VCF',       )
    parser_paired.add_argument('-strelka-indel', '--strelka-indel',       type=str,   help='Strelka VCF',     )
    parser_paired.add_argument('-tnscope',       '--tnscope-vcf',       type=str,   help='TNscope VCF',       )
    
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
    
    ##
    inputParameters['outdir']            = args.output_directory
    inputParameters['ref']               = args.genome_reference
    inputParameters['truth_snv']         = args.truth_snv
    inputParameters['truth_indel']       = args.truth_indel
    inputParameters['classifier_snv']    = args.classifier_snv
    inputParameters['classifier_indel']  = args.classifier_indel
    inputParameters['pass_threshold']    = args.pass_threshold
    inputParameters['lowqual_threshold'] = args.lowqual_threshold
    inputParameters['hom_threshold']     = args.homozygous_threshold
    inputParameters['het_threshold']     = args.heterozygous_threshold
    inputParameters['minMQ']             = args.minimum_mapping_quality
    inputParameters['minBQ']             = args.minimum_base_quality
    inputParameters['mincaller']         = args.minimum_num_callers
    inputParameters['dbsnp']             = args.dbsnp_vcf
    inputParameters['cosmic']            = args.cosmic_vcf
    inputParameters['inclusion']         = args.inclusion_region
    inputParameters['exclusion']         = args.exclusion_region
    
    inputParameters['somaticseq_train']   = args.somaticseq_train
    inputParameters['keep_intermediates'] = args.keep_intermediates
    
    ################## paired sample ################## 
    if parser.parse_args().which == 'paired':
        inputParameters['tbam']          = args.tumor_bam_file
        inputParameters['nbam']          = args.normal_bam_file
        inputParameters['tumor_name']    = args.tumor_sample
        inputParameters['normal_name']   = args.normal_sample
        inputParameters['mutect']        = args.mutect_vcf
        inputParameters['indelocator']   = args.indelocator_vcf
        inputParameters['mutect2']       = args.mutect2_vcf
        inputParameters['varscan_snv']   = args.varscan_snv
        inputParameters['varscan_indel'] = args.varscan_indel
        inputParameters['jsm']           = args.jsm_vcf
        inputParameters['sniper']        = args.somaticsniper_vcf
        inputParameters['vardict']       = args.vardict_vcf
        inputParameters['muse']          = args.muse_vcf
        inputParameters['lofreq_snv']    = args.lofreq_snv
        inputParameters['lofreq_indel']  = args.lofreq_indel
        inputParameters['scalpel']       = args.scalpel_vcf
        inputParameters['strelka_snv']   = args.strelka_snv
        inputParameters['strelka_indel'] = args.strelka_indel
        inputParameters['tnscope']       = args.tnscope_vcf
        inputParameters['mode']          = 'paired'
        
    ################## single sample ##################
    elif parser.parse_args().which == 'single':
        inputParameters['tbam']        = args.bam_file
        inputParameters['sample_name'] = args.sample_name
        inputParameters['mutect']      = args.mutect_vcf
        inputParameters['mutect2']     = args.mutect2_vcf
        inputParameters['varscan']     = args.varscan_vcf
        inputParameters['vardict']     = args.vardict_vcf
        inputParameters['lofreq']      = args.lofreq_vcf
        inputParameters['scalpel']     = args.scalpel_vcf
        inputParameters['strelka']     = args.strelka_vcf

    return inputParameters




################################################################################################
# Execute:
if __name__ == '__main__':
    runParameters = run()
    if runParameters['mode'] == 'paired':
        
        runPaired( outdir             = runParameters['outdir'], \
                   ref                = runParameters['ref'], \
                   tbam               = runParameters['tbam'], \
                   nbam               = runParameters['nbam'], \
                   tumor_name         = runParameters['tumor_name'], \
                   normal_name        = runParameters['normal_name'], \
                   truth_snv          = runParameters['truth_snv'], \
                   truth_indel        = runParameters['truth_indel'], \
                   classifier_snv     = runParameters['classifier_snv'], \
                   classifier_indel   = runParameters['classifier_indel'], \
                   pass_threshold     = runParameters['pass_threshold'], \
                   lowqual_threshold  = runParameters['lowqual_threshold'], \
                   hom_threshold      = runParameters['hom_threshold'], \
                   het_threshold      = runParameters['het_threshold'], \
                   min_mq             = runParameters['minMQ'], \
                   min_bq             = runParameters['minBQ'], \
                   min_caller         = runParameters['mincaller'], \
                   dbsnp              = runParameters['dbsnp'], \
                   cosmic             = runParameters['cosmic'], \
                   inclusion          = runParameters['inclusion'], \
                   exclusion          = runParameters['exclusion'], \
                   mutect             = runParameters['mutect'], \
                   indelocator        = runParameters['indelocator'], \
                   mutect2            = runParameters['mutect2'], \
                   varscan_snv        = runParameters['varscan_snv'], \
                   varscan_indel      = runParameters['varscan_indel'], \
                   jsm                = runParameters['jsm'], \
                   sniper             = runParameters['sniper'], \
                   vardict            = runParameters['vardict'], \
                   muse               = runParameters['muse'], \
                   lofreq_snv         = runParameters['lofreq_snv'], \
                   lofreq_indel       = runParameters['lofreq_indel'], \
                   scalpel            = runParameters['scalpel'], \
                   strelka_snv        = runParameters['strelka_snv'], \
                   strelka_indel      = runParameters['strelka_indel'], \
                   tnscope            = runParameters['tnscope'], \
                   somaticseq_train   = runParameters['somaticseq_train'], \
                   keep_intermediates = runParameters['keep_intermediates'] )
    
    elif runParameters['mode'] == 'single':
        
        runSingle( outdir             = runParameters['outdir'], \
                   ref                = runParameters['ref'], \
                   bam                = runParameters['tbam'], \
                   sample_name        = runParameters['sample_name'], \
                   truth_snv          = runParameters['truth_snv'], \
                   truth_indel        = runParameters['truth_indel'], \
                   classifier_snv     = runParameters['classifier_snv'], \
                   classifier_indel   = runParameters['classifier_indel'], \
                   pass_threshold     = runParameters['pass_threshold'], \
                   lowqual_threshold  = runParameters['lowqual_threshold'], \
                   hom_threshold      = runParameters['hom_threshold'], \
                   het_threshold      = runParameters['het_threshold'], \                   
                   min_mq             = runParameters['minMQ'], \
                   min_bq             = runParameters['minBQ'], \
                   min_caller         = runParameters['mincaller'], \
                   dbsnp              = runParameters['dbsnp'], \
                   cosmic             = runParameters['cosmic'], \
                   inclusion          = runParameters['inclusion'], \
                   exclusion          = runParameters['exclusion'], \
                   mutect             = runParameters['mutect'], \
                   mutect2            = runParameters['mutect2'], \
                   varscan            = runParameters['varscan'], \
                   vardict            = runParameters['vardict'], \
                   lofreq             = runParameters['lofreq'], \
                   scalpel            = runParameters['scalpel'], \
                   strelka            = runParameters['strelka'], \
                   somaticseq_train   = runParameters['somaticseq_train'], \
                   keep_intermediates = runParameters['keep_intermediates'] )
