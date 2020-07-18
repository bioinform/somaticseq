#!/usr/bin/env python3

import sys, os, argparse, gzip, re, subprocess, uuid

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomicFileHandler.genomic_file_handlers as genome



def bed_include(infile, inclusion_region, outfile):
    
    assert infile != outfile
    
    if inclusion_region:
    
        cmd_line = 'intersectBed -header -a {} -b {} | uniq > {}'.format(infile, inclusion_region, outfile)
        subprocess.check_call(cmd_line, shell=True)
        
    else:
        outfile = None
    
    return outfile



def bed_exclude(infile, exclusion_region, outfile):
    
    assert infile != outfile
    
    if exclusion_region:
        cmd_line = 'intersectBed -header -a {} -b {} -v | uniq > {}'.format(infile, exclusion_region, outfile) 
        subprocess.check_call(cmd_line, shell=True)
        
    else:
        outfile = None
        
    return outfile




def bed_intersector(infile, outfile, inclusion_region=None, exclusion_region=None):
    
    assert infile != outfile
    from shutil import copyfile
    
    # Get the input file name minus the extention, and also get the extension
    infile_noext = re.sub(r'\.\w+$', '', infile)
    file_ext     = re.search(r'\.\w+$',  infile).group()
    
    temp_files   = []
    
    if inclusion_region:
        
        included_temp_file = infile_noext + uuid.uuid4().hex + file_ext

        cmd_line = 'intersectBed -header -a {} -b {} | uniq > {}'.format(infile, inclusion_region, included_temp_file)
        subprocess.check_call(cmd_line, shell=True)

        infile  = included_temp_file
        temp_files.append( included_temp_file )
    
    
    if exclusion_region:
        
        cmd_line = 'intersectBed -header -a {} -b {} -v | uniq > {}'.format(infile, exclusion_region, outfile)
        subprocess.check_call(cmd_line, shell=True)
    
    if inclusion_region and not exclusion_region:
        copyfile(included_temp_file, outfile)
    
    elif not (inclusion_region or exclusion_region):
        if infile.endswith('.gz'):
        
            cmd_line = 'gunzip -c {} > {}'.format(infile, outfile)
            subprocess.check_call(cmd_line, shell=True)
        
        else:
            copyfile(infile, outfile)
    
    
    for file_i in temp_files:
        os.remove( file_i )
        
    return outfile




# Use utilities/vcfsorter.pl fa.dict unsorted.vcf > sorted.vcf
def vcfsorter(ref, vcfin, vcfout):
    
    #vcfsort = '{}/utilities/vcfsorter.pl'.format(PRE_DIR)
    #os.system( '{} {} {} > {}'.format(vcfsort, hg_dict, vcfin, vcfout ) )
    fai = ref + '.fai'
    
    cmd_line = 'bedtools sort -faidx {} -header -i {} > {}'.format(fai, vcfin, vcfout)
    subprocess.check_call(cmd_line, shell=True)
