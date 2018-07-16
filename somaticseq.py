#!/usr/bin/env python3

import sys, argparse, gzip, os, re
import genomicFileHandler.genomic_file_handlers as genome
import vcfModifier.copy_TextFile as copy_TextFile

MY_DIR = os.path.dirname(os.path.realpath(__file__))


def runSingle(outdir, ref, bam, sample_name='TUMOR', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, dbsnp=None, cosmic=None, inclusion=None, exclusion=None, mutect=None, mutect2=None, varscan=None, vardict=None, lofreq=None, scalpel=None, strelka=None, keep_intermediates=False):
    pass


def runPaired(outdir, ref, tbam, nbam, tumor_name='TUMOR', normal_name='NORMAL', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, dbsnp=None, cosmic=None, inclusion=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, keep_intermediates=False):
    
    intermediate_files  = []
    snv_intermediates   = []
    indel_intermediates = []
    
    hg_dict = re.sub(r'\.fa(sta)?$', '.dict', ref)
    
    # Modify direct VCF outputs for merging:
    if mutect2:
        import vcfModifier.modify_MuTect2 as mod_mutect2
        
        snv_mutect_out = outdir + os.sep + 'snv.mutect.vcf'
        indel_mutect_out = outdir + os.sep + 'indel.mutect.vcf'
        mod_mutect2.convert(mutect2, snv_mutect_out, indel_mutect_out, False)
        
        intermediate_files.extend( [snv_mutect_out, indel_mutect_out] )
    
    if varscan_snv or varscan_indel:
        import vcfModifier.modify_VarScan2 as mod_varscan2
        
        if varscan_snv:
            snv_varscan_out = outdir + os.sep + 'snv.varscan.vcf'
            mod_varscan2.convert(varscan_snv, snv_varscan_out)
            
        if varscan_indel:
            indel_varscan_out = outdir + os.sep + 'indel.varscan.vcf'
            mod_varscan2.convert(varscan_indel, indel_varscan_out)
            
        intermediate_files.extend( [snv_varscan_out, indel_varscan_out] )
    
    if jsm:
        import vcfModifier.modify_JointSNVMix2 as mod_jsm
        
        jsm_out = outdir + os.sep + 'snv.jsm.vcf'
        mod_jsm.convert(jsm, jsm_out)
        
        intermediate_files.append(jsm_out)
        
    if sniper:
        import vcfModifier.modify_SomaticSniper as mod_sniper
        
        sniper_out = outdir + os.sep + 'snv.somaticsniper.vcf'
        mod_sniper.convert(sniper, sniper_out)
        
        intermediate_files.append(sniper_out)

    if vardict:
        import vcfModifier.modify_VarDict as mod_vardict
        
        snv_vardict_out = outdir + os.sep + 'snv.vardict.vcf'
        indel_vardict_out = outdir + os.sep + 'indel.vardict.vcf'
        mod_vardict.convert(vardict, snv_vardict_out, indel_vardict_out)
        
        sorted_snv_vardict_out = outdir + os.sep + 'snv.sort.vardict.vcf'
        sorted_indel_vardict_out = outdir + os.sep + 'indel.sort.vardict.vcf'
        os.system('{}/utilities/vcfsorter.pl {} {} > {}'.format(MY_DIR, hg_dict, snv_vardict_out,   sorted_snv_vardict_out) )
        os.system('{}/utilities/vcfsorter.pl {} {} > {}'.format(MY_DIR, hg_dict, indel_vardict_out, sorted_indel_vardict_out) )
        
        intermediate_files.extend([snv_vardict_out, indel_vardict_out, sorted_snv_vardict_out, indel_vardict_out])

    if muse:
        muse_out = outdir + os.sep + 'snv.muse.vcf'
        copy_TextFile.copy(muse, muse_out)
        intermediate_files.append( muse_out )
        
    if lofreq_snv:
        snv_lofreq_out = outdir + os.sep + 'snv.lofreq.vcf'
        copy_TextFile.copy(lofreq_snv, snv_lofreq_out)
        intermediate_files.append( snv_lofreq_out )

    if lofreq_indel:
        indel_lofreq_out = outdir + os.sep + 'indel.lofreq.vcf'
        copy_TextFile.copy(lofreq_indel, indel_lofreq_out)
        intermediate_files.append( indel_lofreq_out )
        
    if scalpel:
        scalpel_out = outdir + os.sep + 'indel.scalpel.vcf'
        copy_TextFile.copy(scalpel, scalpel_out)
        intermediate_files.append( scalpel_out )
        





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
    
    parser.add_argument('-dbsnp',  '--dbsnp-vcf',          type=str,   help='dbSNP VCF',)
    parser.add_argument('-cosmic', '--cosmic-vcf',         type=str,   help='COSMIC VCF')
    
    parser.add_argument('-include',  '--inclusion-region', type=str,   help='inclusion bed')
    parser.add_argument('-exclude',  '--exclusion-region', type=str,   help='exclusion bed')
    
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
    inputParameters['dbsnp']             = args.dbsnp_vcf
    inputParameters['cosmic']            = args.cosmic_vcf
    inputParameters['inclusion']         = args.inclusion_region
    inputParameters['exclusion']         = args.exclusion_region
    
    inputParameters['keep_intermediates'] = args.keep_intermediates
    
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
                   keep_intermediates = runParameters['keep_intermediates'] )
