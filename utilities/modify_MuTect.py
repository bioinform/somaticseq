#!/usr/bin/env python3

# For single-sample mode, remove the "none" column from the VCF file to be consistent with VCF files from other single-sample tools. 
# Keep Broad's convention for AD, i.e., allelic depths for the ref and alt alleles in the order listed
# 4/12/2015


import sys, os, argparse, gzip
import regex as re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

# argparse Stuff
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Variant Call Type, i.e., snp or indel
parser.add_argument('-type',     '--variant-type',         type=str, help='Either SNP or INDEL. Required', required=True, default='snp')
parser.add_argument('-infile',   '--input-vcf', nargs='*', type=str, help='Input VCF file', required=False, default=None)
parser.add_argument('-outfile',  '--output-vcf',           type=str, help='Output VCF file', required=False, default=sys.stdout)
parser.add_argument('-tag',      '--somatic-tag',          type=str, required=False, help='Change the SOMATIC tag to something else', default='CGA')
parser.add_argument('-tbam',     '--tumor-bam',            type=str, required=False, help='A tumor bam file for sample name identification.' )
parser.add_argument('-nbam',     '--normal-bam',           type=str, required=False, help='A normal bam file for sample name identification.' )
parser.add_argument('-tsm',      '--tumor-input-name',     type=str, required=False, help='Tumor sample name in the MuTect vcf file.' )
parser.add_argument('-nsm',      '--normal-input-name',    type=str, required=False, help='Normal sample name in the MuTect vcf file.' )
parser.add_argument('-N',        '--normal-sample-name',   type=str, help='1st Sample Name, default=NORMAL',  required=False, default='NORMAL')
parser.add_argument('-T',        '--tumor-sample-name',    type=str, help='2nd Sample Name, default=TUMOR',   required=False, default='TUMOR')
parser.add_argument('-samtools', '--samtools-command',     type=str, required=False, help='Deprecated', default='samtools' )
parser.add_argument('-gz',       '--gz-compressed', action='store_true', help='If the input files are 1.vcf.gz, 2.vcf.gz, ...', required=False)


# Parse the arguments:
args = parser.parse_args()

samtools = args.samtools_command
out_vcf = args.output_vcf

# the -gz flag appropriate IF AND ONLY IF I did not directly supply a list of vcf file(s):
gz = '.gz' if (not args.input_vcf) and args.gz_compressed else ''


# VCF File's index locations
idx_chrom,idx_pos,idx_id,idx_ref,idx_alt,idx_qual,idx_filter,idx_info,idx_format = 0,1,2,3,4,5,6,7,8
idx_SM1, idx_SM2 = 9,10

# The regular expression pattern for "chrXX 1234567" in both VarScan2 Output and VCF files:
pattern_chr_position = re.compile(r'^(?:chr)?(?:[0-9]+|[XYM]|MT)\t[0-9]+\b')


tbam = args.tumor_bam if args.tumor_bam else None
nbam = args.normal_bam if args.normal_bam else None


if tbam:

    paired_mode = True if nbam else False
    
    # Get tumor and normal sample names from the bam files:
    nbam_header = genome.pysam_header(nbam) if nbam else None
    tbam_header = genome.pysam_header(tbam)
    
    # When MuTect is run in a "single sample mode," the "normal" will be named "none."
    n_samplename = nbam_header.SM() if nbam else ['none']
    t_samplename = tbam_header.SM()
    
    if not ( len(n_samplename)==1 and len(t_samplename)==1 ):
        sys.stderr.write('There are multiple Sample Names present in the BAM file!')
    
    n_samplename = n_samplename[0]
    t_samplename = t_samplename[0]
    
elif args.tumor_input_name:
    
    paired_mode = True if args.normal_input_name else False
    n_samplename = args.normal_input_name if args.normal_input_name else 'none'
    t_samplename = args.tumor_input_name

else:
    raise Exception('Cannot tell which column is Tumor and which column is normal.')



def open_my_vcf(file_name):
    
    # See if the input file is a .gz file:
    if file_name.lower().endswith('.gz'):
        return gzip.open(file_name, 'rt')
        
    else:
        return open(file_name)



if not args.input_vcf:

    # Chromosome namings:
    hg19_chrom_sequence = ['chr'+str(i) for i in range(1,23)]
    hg19_chrom_sequence.append('chrX')
    hg19_chrom_sequence.append('chrY')
    hg19_chrom_sequence.append('chrM')
    b37_chrom_sequence = [ i.replace('chr','').replace('M', 'MT') for i in hg19_chrom_sequence ]

    # Basically, decide to read the files as (chr1, chr2, ..., chrX, chrY, and chrM), or (1, 2, ..., X, Y, MT). 
    local_files = os.listdir()
    hg19_file_check = [ os.path.exists(i+'.vcf'+gz) for i in hg19_chrom_sequence ]
    hg19_bool = [True for i in hg19_chrom_sequence]
    got_hg19 = hg19_file_check==hg19_bool
    
    b37_file_check  = [ os.path.exists(i+'.vcf'+gz) for i in b37_chrom_sequence ]
    b37_bool = [True for i in b37_file_check]
    got_b37 = b37_file_check==b37_bool
    
    if not (got_hg19) and not (got_b37):
        raise Exception('chrX.vcf or X.vcf does not exist.')
    elif got_hg19 and got_b37:
        raise Exception('Both chrX.vcf and X.vcf exist. I cannot decide.')
    elif got_hg19:
        chosen_chrom_sequence = hg19_chrom_sequence
    elif got_b37:
        chosen_chrom_sequence = b37_chrom_sequence
        
    right_files = [i+'.vcf'+gz for i in chosen_chrom_sequence]

else:
    right_files = (args.input_vcf)



# Now, finally start merging. 
# First, get the header:
vcf_header = []
with open_my_vcf(right_files[0]) as first_vcf:
    
    line_i = first_vcf.readline()
    
    while line_i.startswith('#'):
        
        if line_i.startswith('##INFO=<ID=SOMATIC,'):
            
            replaced_line = line_i.replace('SOMATIC', args.somatic_tag).replace('Somatic event', 'CGA called somatic event')
            vcf_header.append( replaced_line )
            
        elif line_i.startswith('##FORMAT=<ID=BQ,') or line_i.startswith('##FORMAT=<ID=FA,'):
            
            replaced_line = line_i.replace('Number=A', 'Number=.')
            vcf_header.append( replaced_line )
        
        elif line_i.startswith('##contig=') or line_i.startswith('##GATKCommandLine=') or line_i.startswith('##SID_bam_file'):
            pass
        
        elif line_i.startswith('##'):
            vcf_header.append(line_i)
            
        elif line_i.startswith('#CHROM'):
            header_items = line_i.rstrip().split('\t')
            
            idxN = header_items.index(n_samplename)
            idxT = header_items.index(t_samplename)
                        
            if paired_mode:
                header_items[idx_SM1] = args.normal_sample_name
                header_items[idx_SM2] = args.tumor_sample_name
                
            else:
                
                # Keep up to the first sample column, then make sure it's labeled the TUMOR sample name
                header_items = header_items[:idx_SM1+1]
                header_items[idx_SM1] = args.tumor_sample_name
            
            replaced_header = '\t'.join(header_items) + '\n'
            vcf_header.append(replaced_header)
        
        line_i = first_vcf.readline()


# Open file and start doing work:
vcfout = open(out_vcf, 'w')

# Start printing the headers first:
for line_i in vcf_header:
    vcfout.write(line_i)



for chr_i_vcf in right_files:
    
    with open_my_vcf(chr_i_vcf) as vcf:
        
        line_i = vcf.readline().rstrip()
        
        # Skip headers from now on:
        while line_i.startswith('#'):
            line_i = vcf.readline().rstrip()
            
        
        while line_i:
            
            items_i = line_i.split('\t')
            
            if 'SOMATIC' in items_i[idx_info]:
                items_i[idx_info] = items_i[idx_info].replace('SOMATIC', args.somatic_tag)
            
            if paired_mode:
                items_i[idx_SM1], items_i[idx_SM2] = items_i[idxN], items_i[idxT]
            
            else:
                items_i = items_i[:idx_SM1] + [items_i[idxT]]
                        
            # Print the new stuff:
            new_line = '\t'.join( items_i )
            
            # Have to get rid of "N" in REF, because after snpSift annotation, it changes the ALT and vcf-validator will complain.
            if not ( 'N' in items_i[idx_ref] ):
                vcfout.write( new_line + '\n' )
            
            
            line_i = vcf.readline().rstrip()


vcfout.close()
