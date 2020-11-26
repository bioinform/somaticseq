#!/usr/bin/env python3

import os, logging
from multiprocessing import Pool
from functools import partial
from shutil import rmtree

import somaticseq.run_somaticseq as run_somaticseq
import somaticseq.utilities.split_Bed_into_equal_regions as split_bed
import somaticseq.genomicFileHandler.concat as concat



def splitRegions(nthreads, outfiles, bed=None, fai=None):

    assert bed or fai
    if fai and not bed:
        bed = split_bed.fai2bed(fai, outfiles)

    writtenBeds = split_bed.split(bed, outfiles, nthreads)

    return writtenBeds




def runPaired_by_region(inclusion, outdir=None, ref=None, tbam=None, nbam=None, tumor_name='TUMOR', normal_name='NORMAL', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, platypus=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', algo='ada', keep_intermediates=False, train_seed=0, tree_depth=12, iterations=200, features_excluded=[]):

    logger = logging.getLogger(runPaired_by_region.__name__)

    basename   = inclusion.split(os.sep)[-1].split('.')[0]
    outdir_i   = outdir + os.sep + basename
    os.makedirs(outdir_i, exist_ok=True)

    run_somaticseq.runPaired(outdir_i, ref, tbam, nbam, tumor_name, normal_name, truth_snv, truth_indel, classifier_snv, classifier_indel, pass_threshold, lowqual_threshold, hom_threshold, het_threshold, dbsnp, cosmic, inclusion, exclusion, mutect, indelocator, mutect2, varscan_snv, varscan_indel, jsm, sniper, vardict, muse, lofreq_snv, lofreq_indel, scalpel, strelka_snv, strelka_indel, tnscope, platypus, min_mq, min_bq, min_caller, somaticseq_train, ensembleOutPrefix, consensusOutPrefix, classifiedOutPrefix, algo, keep_intermediates, train_seed, tree_depth, iterations, features_excluded)

    return outdir_i



def runSingle_by_region(inclusion, outdir, ref, bam, sample_name='TUMOR', truth_snv=None, truth_indel=None, classifier_snv=None, classifier_indel=None, pass_threshold=0.5, lowqual_threshold=0.1, hom_threshold=0.85, het_threshold=0.01, dbsnp=None, cosmic=None, exclusion=None, mutect=None, mutect2=None, varscan=None, vardict=None, lofreq=None, scalpel=None, strelka=None, min_mq=1, min_bq=5, min_caller=0.5, somaticseq_train=False, ensembleOutPrefix='Ensemble.', consensusOutPrefix='Consensus.', classifiedOutPrefix='SSeq.Classified.', algo='ada', keep_intermediates=False, train_seed=0, tree_depth=12, iterations=200, features_excluded=[]):

    logger = logging.getLogger(runSingle_by_region.__name__)
    
    basename   = inclusion.split(os.sep)[-1].split('.')[0]
    outdir_i   = outdir + os.sep + basename
    os.makedirs(outdir_i, exist_ok=True)

    run_somaticseq.runSingle(outdir_i, ref, bam, sample_name, truth_snv, truth_indel, classifier_snv, classifier_indel, pass_threshold, lowqual_threshold, hom_threshold, het_threshold, dbsnp, cosmic, inclusion, exclusion, mutect, mutect2, varscan, vardict, lofreq, scalpel, strelka, min_mq, min_bq, min_caller, somaticseq_train, ensembleOutPrefix, consensusOutPrefix, classifiedOutPrefix, algo, keep_intermediates, train_seed, tree_depth, iterations, features_excluded)

    return outdir_i



def mergeSubdirTsv(dirList, filename, outdir=os.curdir):
    fileList = ['{}/{}'.format(dir_i, filename) for dir_i in dirList]
    concat.tsv(fileList, outdir + os.sep + filename)

def mergeSubdirVcf(dirList, filename, outdir=os.curdir):
    fileList = ['{}/{}'.format(dir_i, filename) for dir_i in dirList]
    concat.vcf(fileList, outdir + os.sep + filename)



if __name__ == '__main__':

    args = run_somaticseq.run()

    os.makedirs(args.output_directory, exist_ok=True)

    bed_splitted = splitRegions(args.threads, args.output_directory+os.sep+'th.input.bed', args.inclusion_region, args.genome_reference+'.fai')

    pool = Pool(processes = args.threads)

    if args.which == 'paired':

        runPaired_by_region_i = partial(runPaired_by_region, \
                   outdir             = args.output_directory, \
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
                   somaticseq_train   = False, \
                   train_seed         = args.seed, \
                   tree_depth         = args.tree_depth, \
                   iterations         = args.iterations, \
                   features_excluded  = args.features_excluded, \
                   keep_intermediates = args.keep_intermediates, \
                   )

        subdirs = pool.map(runPaired_by_region_i, bed_splitted)
        pool.close()

    elif args.which == 'single':

        runSingle_by_region_i = partial(runSingle_by_region, \
                   outdir             = args.output_directory, \
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
                   exclusion          = args.exclusion_region, \
                   mutect             = args.mutect_vcf, \
                   mutect2            = args.mutect2_vcf, \
                   varscan            = args.varscan_vcf, \
                   vardict            = args.vardict_vcf, \
                   lofreq             = args.lofreq_vcf, \
                   scalpel            = args.scalpel_vcf, \
                   strelka            = args.strelka_vcf, \
                   algo               = args.algorithm, \
                   somaticseq_train   = False, \
                   train_seed         = args.seed, \
                   tree_depth         = args.tree_depth, \
                   iterations         = args.iterations, \
                   features_excluded  = args.features_excluded, \
                   keep_intermediates = args.keep_intermediates, \
                   )

        subdirs = pool.map(runSingle_by_region_i, bed_splitted)
        pool.close()


    run_somaticseq.logger.info('Sub-directories created: {}'.format(', '.join(subdirs)) )

    # Merge sub-results
    mergeSubdirTsv(subdirs, 'Ensemble.sSNV.tsv', args.output_directory)
    mergeSubdirTsv(subdirs, 'Ensemble.sINDEL.tsv', args.output_directory)

    if args.classifier_snv:
        mergeSubdirTsv(subdirs, 'SSeq.Classified.sSNV.tsv', args.output_directory)
        mergeSubdirVcf(subdirs, 'SSeq.Classified.sSNV.vcf', args.output_directory)
    else:
        mergeSubdirVcf(subdirs, 'Consensus.sSNV.vcf', args.output_directory)

    if args.classifier_indel:
        mergeSubdirTsv(subdirs, 'SSeq.Classified.sINDEL.tsv', args.output_directory)
        mergeSubdirVcf(subdirs, 'SSeq.Classified.sINDEL.vcf', args.output_directory)
    else:
        mergeSubdirVcf(subdirs, 'Consensus.sINDEL.vcf', args.output_directory)

    # If there is training, it should be done after merging the results
    if args.somaticseq_train:

        snv_training_file   = args.output_directory + os.sep + 'Ensemble.sSNV.tsv'
        indel_training_file = args.output_directory + os.sep + 'Ensemble.sINDEL.tsv'
        
        num_iterations = args.iterations if args.iterations else run_somaticseq.DEFAULT_XGB_BOOST_ROUNDS
            
        run_somaticseq.modelTrainer(snv_training_file,   args.algorithm, threads=args.threads, seed=args.seed, max_depth=args.tree_depth, iterations=num_iterations, features_to_exclude=args.features_excluded)
        run_somaticseq.modelTrainer(indel_training_file, args.algorithm, threads=args.threads, seed=args.seed, max_depth=args.tree_depth, iterations=num_iterations, features_to_exclude=args.features_excluded)
        

    # Clean up after yourself
    if not args.keep_intermediates:
        for bed_i in bed_splitted:
            os.remove( bed_i )
            run_somaticseq.logger.info('Removed: {}'.format( bed_i ) )

        for dir_i in subdirs:
            rmtree( dir_i )
            run_somaticseq.logger.info('Removed sub-directory: {}'.format( dir_i ) )


    run_somaticseq.logger.info('Done.')
