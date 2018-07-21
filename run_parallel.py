#!/usr/bin/env python3

import sys, os, argparse, shutil, math, re
from multiprocessing import Pool
import somaticseq.run_somaticseq as run_somaticseq
import utilities.split_Bed_into_equal_regions as split_bed




def splitRegions(nthreads, outfiles, bed=None, fai=None):
    
    if bed:
        writtenBeds = split_bed.split(bed, outfiles, nthreads)
    elif fai:
        bed         = split_bed.fai2bed(fai, outfiles)
        writtenBeds = split_bed.split(bed, outfiles, nthreads)
        
    return writtenBeds



if __name__ == '__main__':

    runParameters = run_somaticseq.run()
    bed_splitted = splitRegions(runParameters['threads'], runParameters['outdir']+os.sep+'th.input.bed', runParameters['inclusion'], runParameters['ref']+'.fai')

    print( bed_splitted )
