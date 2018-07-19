#!/usr/bin/env python3

import sys, os, argparse, gzip, re, subprocess

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomicFileHandler.genomic_file_handlers as genome



def bed_include(infile, inclusion_region, outfile):
    
    assert infile != outfile
    
    if inclusion_region:
    
        exit_code = os.system( 'intersectBed -header -a {} -b {} | uniq > {}'.format(infile, inclusion_region, outfile) )
        assert exit_code == 0
    else:
    
        outfile = infile
    
    return outfile



def bed_exclude(infile, exclusion_region, outfile):
    
    assert infile != outfile
    
    if exclusion_region:
    
        exit_code = os.system( 'intersectBed -header -a {} -b {} -v | uniq > {}'.format(infile, exclusion_region, outfile) )
        assert exit_code == 0
    else:
    
        outfile = infile
        
    return outfile


# Use utilities/vcfsorter.pl fa.dict unsorted.vcf > sorted.vcf
def vcfsorter(hg_dict, vcfin, vcfout):
    
    vcfsort = '{}/utilities/vcfsorter.pl'.format(PRE_DIR)
    os.system( '{} {} {} > {}'.format(vcfsort, hg_dict, vcfin, vcfout ) )
    
