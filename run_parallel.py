#!/usr/bin/env python3

import sys, os, argparse, shutil, math, re, subprocess
from multiprocessing import Pool
from functools import partial

import somaticseq.run_somaticseq as run_somaticseq
import utilities.split_Bed_into_equal_regions as split_bed
import genomicFileHandler.concat as concat



def splitRegions(nthreads, outfiles, bed=None, fai=None):

    if bed:
        writtenBeds = split_bed.split(bed, outfiles, nthreads)
    elif fai:
        bed         = split_bed.fai2bed(fai, outfiles)
        writtenBeds = split_bed.split(bed, outfiles, nthreads)

    return writtenBeds




def runPaired_by_region(inclusion, outdir=None, ref=None, tbam=None, nbam=None, tumor_name='TUMOR', normal_name='NORMAL', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', keep_intermediates=False):

    basename   = inclusion.split(os.sep)[-1].split('.')[0]
    outdir_i   = outdir + os.sep + basename
    os.mkdir(outdir_i)

    run_somaticseq.runPaired(outdir_i, ref, tbam, nbam, tumor_name, normal_name, truth_snv, truth_indel, classifier_snv, classifier_indel, pass_threshold, lowqual_threshold, hom_threshold, het_threshold, dbsnp, cosmic, inclusion, exclusion, mutect, indelocator, mutect2, varscan_snv, varscan_indel, jsm, sniper, vardict, muse, lofreq_snv, lofreq_indel, scalpel, strelka_snv, strelka_indel, tnscope, min_mq, min_bq, min_caller, somaticseq_train, ensembleOutPrefix, consensusOutPrefix, classifiedOutPrefix, keep_intermediates)

    return outdir_i



def runSingle_by_region(inclusion, outdir, ref, bam, sample_name='TUMOR', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, exclusion=None, mutect=None, mutect2=None, varscan=None, vardict=None, lofreq=None, scalpel=None, strelka=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', keep_intermediates=False):

    basename   = inclusion.split(os.sep)[-1].split('.')[0]
    outdir_i   = outdir + os.sep + basename
    os.mkdir(outdir_i)

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

    bed_splitted = splitRegions(runParameters['threads'], runParameters['outdir']+os.sep+'th.input.bed', runParameters['inclusion'], runParameters['ref']+'.fai')

    pool = Pool(processes = runParameters['threads'])

    if runParameters['mode'] == 'paired':

        runPaired_by_region_i = partial(runPaired_by_region, \
                   outdir             = runParameters['outdir'], \
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
                   somaticseq_train   = False, \
                   keep_intermediates = runParameters['keep_intermediates'] )

        subdirs = pool.map(runPaired_by_region_i, bed_splitted)

    elif runParameters['mode'] == 'single':

        runSingle_by_region_i = partial(runSingle_by_region, \
                   outdir             = runParameters['outdir'], \
                   ref                = runParameters['ref'], \
                   bam                = runParameters['bam'], \
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
                   exclusion          = runParameters['exclusion'], \
                   mutect             = runParameters['mutect'], \
                   mutect2            = runParameters['mutect2'], \
                   varscan            = runParameters['varscan'], \
                   vardict            = runParameters['vardict'], \
                   lofreq             = runParameters['lofreq'], \
                   scalpel            = runParameters['scalpel'], \
                   strelka            = runParameters['strelka'], \
                   somaticseq_train   = False, \
                   keep_intermediates = runParameters['keep_intermediates'] )

        subdirs = pool.map(runSingle_by_region_i, bed_splitted)

    print('Sub-directories created: {}'.format(subdirs), file=sys.stderr)

    # Merge sub-results
    mergeSubdirTsv(subdirs, 'Ensemble.sSNV.tsv', runParameters['outdir'])
    mergeSubdirTsv(subdirs, 'Ensemble.sINDEL.tsv', runParameters['outdir'])

    if runParameters['classifier_snv']:
        mergeSubdirTsv(subdirs, 'SSeq.Classified.sSNV.tsv', runParameters['outdir'])
        mergeSubdirVcf(subdirs, 'SSeq.Classified.sSNV.vcf', runParameters['outdir'])
    else:
        mergeSubdirVcf(subdirs, 'Consensus.sSNV.vcf', runParameters['outdir'])

    if runParameters['classifier_indel']:
        mergeSubdirTsv(subdirs, 'SSeq.Classified.sINDEL.tsv', runParameters['outdir'])
        mergeSubdirVcf(subdirs, 'SSeq.Classified.sINDEL.vcf', runParameters['outdir'])
    else:
        mergeSubdirVcf(subdirs, 'Consensus.sINDEL.vcf', runParameters['outdir'])

    if runParameters['somaticseq_train']:
        subprocess.call( (run_somaticseq.adaTrainer, runParameters['outdir'] + os.sep + 'Ensemble.sSNV.tsv', 'Consistent_Mates' 'Inconsistent_Mates') )
        subprocess.call( (run_somaticseq.adaTrainer, runParameters['outdir'] + os.sep + 'Ensemble.sINDEL.tsv', 'Strelka_QSS', 'Strelka_TQSS','Consistent_Mates' 'Inconsistent_Mates') )
