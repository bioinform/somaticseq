#!/usr/bin/env python3

import sys, os, argparse, shutil, math, re, subprocess
from multiprocessing import Pool
from functools import partial
from shutil import rmtree

import somaticseq.run_somaticseq as run_somaticseq
import utilities.split_Bed_into_equal_regions as split_bed
import genomicFileHandler.concat as concat

def splitRegions(nthreads, outfiles, bed=None, fai=None):

    assert bed or fai
    if fai and not bed:
        bed = split_bed.fai2bed(fai, outfiles)

    writtenBeds = split_bed.split(bed, outfiles, nthreads)

    return writtenBeds




def runPaired_by_region(inclusion, outdir=None, ref=None, tbam=None, nbam=None, tumor_name='TUMOR', normal_name='NORMAL', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, platypus=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', keep_intermediates=False):

    basename   = inclusion.split(os.sep)[-1].split('.')[0]
    outdir_i   = outdir + os.sep + basename
    os.makedirs(outdir_i, exist_ok=True)

    run_somaticseq.runPaired(outdir_i, ref, tbam, nbam, tumor_name, normal_name, truth_snv, truth_indel, classifier_snv, classifier_indel, pass_threshold, lowqual_threshold, hom_threshold, het_threshold, dbsnp, cosmic, inclusion, exclusion, mutect, indelocator, mutect2, varscan_snv, varscan_indel, jsm, sniper, vardict, muse, lofreq_snv, lofreq_indel, scalpel, strelka_snv, strelka_indel, tnscope, platypus, min_mq, min_bq, min_caller, somaticseq_train, ensembleOutPrefix, consensusOutPrefix, classifiedOutPrefix, keep_intermediates)

    return outdir_i



def runSingle_by_region(inclusion, outdir, ref, bam, sample_name='TUMOR', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, exclusion=None, mutect=None, mutect2=None, varscan=None, vardict=None, lofreq=None, scalpel=None, strelka=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', keep_intermediates=False):

    basename   = inclusion.split(os.sep)[-1].split('.')[0]
    outdir_i   = outdir + os.sep + basename
    os.makedirs(outdir_i, exist_ok=True)

    run_somaticseq.runSingle(outdir_i, ref, bam, sample_name, truth_snv, truth_indel, classifier_snv, classifier_indel, pass_threshold, lowqual_threshold, hom_threshold, het_threshold, dbsnp, cosmic, inclusion, exclusion, mutect, mutect2, varscan, vardict, lofreq, scalpel, strelka, min_mq, min_bq, min_caller, somaticseq_train, ensembleOutPrefix, consensusOutPrefix, classifiedOutPrefix, keep_intermediates)

    return outdir_i



def mergeSubdirTsv(dirList, filename, outdir=os.curdir):
    fileList = ['{}/{}'.format(dir_i, filename) for dir_i in dirList]
    concat.tsv(fileList, outdir + os.sep + filename)

def mergeSubdirVcf(dirList, filename, outdir=os.curdir):
    fileList = ['{}/{}'.format(dir_i, filename) for dir_i in dirList]
    concat.vcf(fileList, outdir + os.sep + filename)



if __name__ == '__main__':

    runParameters = run_somaticseq.run()

    os.makedirs(runParameters['output_directory'], exist_ok=True)

    bed_splitted = splitRegions(runParameters['threads'], runParameters['output_directory']+os.sep+'th.input.bed', runParameters['inclusion_region'], runParameters['genome_reference']+'.fai')

    pool = Pool(processes = runParameters['threads'])

    if runParameters['which'] == 'paired':

        runPaired_by_region_i = partial(runPaired_by_region, \
                   outdir             = runParameters['output_directory'], \
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
                   somaticseq_train   = False, \
                   keep_intermediates = runParameters['keep_intermediates'] )

        subdirs = pool.map(runPaired_by_region_i, bed_splitted)

    elif runParameters['which'] == 'single':

        runSingle_by_region_i = partial(runSingle_by_region, \
                   outdir             = runParameters['output_directory'], \
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
                   exclusion          = runParameters['exclusion_region'], \
                   mutect             = runParameters['mutect_vcf'], \
                   mutect2            = runParameters['mutect2_vcf'], \
                   varscan            = runParameters['varscan_vcf'], \
                   vardict            = runParameters['vardict_vcf'], \
                   lofreq             = runParameters['lofreq_vcf'], \
                   scalpel            = runParameters['scalpel_vcf'], \
                   strelka            = runParameters['strelka_vcf'], \
                   somaticseq_train   = False, \
                   keep_intermediates = runParameters['keep_intermediates'] )

        subdirs = pool.map(runSingle_by_region_i, bed_splitted)

    run_somaticseq.logger.info('Sub-directories created: {}'.format(', '.join(subdirs)) )

    # Merge sub-results
    mergeSubdirTsv(subdirs, 'Ensemble.sSNV.tsv', runParameters['output_directory'])
    mergeSubdirTsv(subdirs, 'Ensemble.sINDEL.tsv', runParameters['output_directory'])

    if runParameters['classifier_snv']:
        mergeSubdirTsv(subdirs, 'SSeq.Classified.sSNV.tsv', runParameters['output_directory'])
        mergeSubdirVcf(subdirs, 'SSeq.Classified.sSNV.vcf', runParameters['output_directory'])
    else:
        mergeSubdirVcf(subdirs, 'Consensus.sSNV.vcf', runParameters['output_directory'])

    if runParameters['classifier_indel']:
        mergeSubdirTsv(subdirs, 'SSeq.Classified.sINDEL.tsv', runParameters['output_directory'])
        mergeSubdirVcf(subdirs, 'SSeq.Classified.sINDEL.vcf', runParameters['output_directory'])
    else:
        mergeSubdirVcf(subdirs, 'Consensus.sINDEL.vcf', runParameters['output_directory'])

    if runParameters['somaticseq_train']:
        subprocess.call( (run_somaticseq.adaTrainer, runParameters['output_directory'] + os.sep + 'Ensemble.sSNV.tsv', 'Consistent_Mates', 'Inconsistent_Mates') )
        subprocess.call( (run_somaticseq.adaTrainer, runParameters['output_directory'] + os.sep + 'Ensemble.sINDEL.tsv', 'Strelka_QSS', 'Strelka_TQSS','Consistent_Mates', 'Inconsistent_Mates') )

    # Clean up after yourself
    if not runParameters['keep_intermediates']:
        for bed_i in bed_splitted:
            os.remove( bed_i )
            run_somaticseq.logger.info('Removed: {}'.format( bed_i ) )

        for dir_i in subdirs:
            rmtree( dir_i )
            run_somaticseq.logger.info('Removed sub-directory: {}'.format( dir_i ) )
