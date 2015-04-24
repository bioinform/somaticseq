#!/usr/bin/env python3
# Update for VarDict. Since VarDict may be in flux, this script may not be stable
# Usage: python3 merge.modify_vcfs_for_gatk_combinevariants_r20140808.py -method VarDict -infile sorted.vardict.vcf.gz -outfile vardict.mod.vcf

# For VarDict, consider a p-value < 0.1 as PASS. 

import sys, os
import argparse
import gzip
import regex as re

sys.path.append('/home/ltfang/programming/NGS')
sys.path.append('/home/lethalfang/programming/NGS')
import genomic_file_handlers as genome


# argparse Stuff
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Variant Call Type, i.e., snp or indel
parser.add_argument('-infile', '--input-vcf', nargs='*', type=str, help='Input VCF file', required=False, default=None)
parser.add_argument('-outfile', '--output-vcf', type=str, help='Output VCF file', required=False, default=sys.stdout)
parser.add_argument('-method', '--call-method', type=str, help='VarScan2, JointSNVMix2, SomaticSniper, or VarDict', required=True)

parser.add_argument('-filter', '--filter-version', type=str, help='strict, medium, or lenient for VarDict', required=False, default='v0')


parser.add_argument('-gz', '--gz-compressed', type=int, required=False, help='1 if compressed. 0 if not (vcf directly). This is only necessary when you want the program to look for file names without providing it. At 1, it will look for 1.vcf.gz, 2.vcf.gz, ... and without it, it will look for 1.vcf, 2.vcf, ...', default=0)

# Parse the arguments:
args = parser.parse_args()

# MuTect has it own script:
assert args.call_method in ('VarScan2', 'JointSNVMix2', 'SomaticSniper', 'VarDict', 'VarScan2_single')


out_vcf = args.output_vcf

# If it's VarDict, split the output into snp and indel separately:
if args.call_method == 'VarDict':
    
    out_snp_file       = out_vcf.split(os.sep)
    out_snp_file[-1]   = 'snp.' + out_snp_file[-1]
    
    out_indel_file     = out_vcf.split(os.sep)
    out_indel_file[-1] = 'indel.' + out_indel_file[-1]
        
    out_snp            = os.sep.join(out_snp_file)
    out_indel          = os.sep.join(out_indel_file)




# VCF File's index locations
idx_chrom,idx_pos,idx_id,idx_ref,idx_alt,idx_qual,idx_filter,idx_info,idx_format,idxN,idxT = 0,1,2,3,4,5,6,7,8,9,10


# the -gz flag only appropriate if I did not directly supply a list of vcf files:
if (not args.input_vcf) and args.gz_compressed:
    gz = '.gz'
else:
    gz = ''


# The regular expression pattern for "chrXX 1234567" in both VarScan2 Output and VCF files:
pattern_chr_position = re.compile(r'^(?:chr)?(?:[0-9]+|[XYM]|MT)\t[0-9]+\b')


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

    if not (got_hg19) and not(got_b37):
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



# Open files:
if args.call_method == 'VarDict':
    snpout   = open(out_snp,   'w')
    indelout = open(out_indel, 'w')
else:
    vcfout   = open(out_vcf,   'w')



### First, get the header. IF there are multiple VCF files (e.g., chromosome by chromosome), use the first file for this purpose:
vcf_header = []
with open_my_vcf(right_files[0]) as vcf:
    
    line_i = vcf.readline().rstrip()
    
    # Save the headers and then sort them:
    vcfheader_filter_info_filter = []
    vcfheader_filter_info_filter.append('##INFO=<ID={0},Number=0,Type=Flag,Description="Indicates if record is a {0} called somatic mutation">'.format(args.call_method))
    
    vcfheader_misc = []
    
    while line_i.startswith('#'):
        
        if re.match(r'##fileformat=', line_i):
            vcffileformat = line_i
        
        elif re.match(r'^##FORMAT=<ID=DP4,', line_i):
            line_i = '##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">'
            vcfheader_filter_info_filter.append( line_i )
            
        elif re.match(r'^##FORMAT=<ID=RD,', line_i):
            line_i = '##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">'
            vcfheader_filter_info_filter.append( line_i )
            
        elif re.match(r'^##FORMAT=<ID=AD,', line_i):
            line_i = '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">'
            vcfheader_filter_info_filter.append( line_i )
        
        elif re.match(r'^##FORMAT=<ID=(BQ|FA),', line_i):
            line_i = line_i.replace('Number=A', 'Number=.')
            vcfheader_filter_info_filter.append( line_i )
            
            
        #### Fix VarDict header ####
        elif re.match(r'^##INFO=<ID=(LSEQ|RSEQ),', line_i):
            line_i = line_i.replace('Number=G', 'Number=1')
            vcfheader_filter_info_filter.append( line_i )
        
        elif line_i.startswith('##FORMAT=<ID=BIAS,'):
            line_i = line_i.replace('Number=1', 'Number=.')
            vcfheader_filter_info_filter.append( line_i )
            
        elif line_i.startswith('##FORMAT=<ID=PSTD,') or \
        line_i.startswith('##FORMAT=<ID=QSTD,') or \
        line_i.startswith('##INFO=<ID=SOR,'):
            line_i = line_i.replace('Type=Float', 'Type=String')
            vcfheader_filter_info_filter.append( line_i )

        
        elif line_i.startswith('##INFO=<ID=SOMATIC,'):
            pass
        
        elif re.match(r'^##(INFO|FORMAT|FILTER)', line_i):
            vcfheader_filter_info_filter.append( line_i )
                    
        elif re.match(r'^##', line_i):
            vcfheader_misc.append( line_i )
            
        ##### Done fixing VarDict header #####
            
        


        elif re.match(r'^#CHROM', line_i):
            header_main_item = line_i.split('\t')
            
            if args.call_method == 'VarScan2_single':
                header_main_item[idxN]='TUMOR'
            
            else:
                header_main_item[idxN]='NORMAL'
                header_main_item[idxT]='TUMOR'
            
            vcfheader_main = '\t'.join(header_main_item)

        # Continue:
        line_i = vcf.readline().rstrip()
    
    
    # Add those things for VarDict stuff:
    if args.call_method == 'VarDict':
        vcfheader_filter_info_filter.append('##INFO=<ID=Germline,Number=0,Type=Flag,Description="VarDict Germline">')
        vcfheader_filter_info_filter.append('##INFO=<ID=StrongSomatic,Number=0,Type=Flag,Description="VarDict Strong Somatic">')
        vcfheader_filter_info_filter.append('##INFO=<ID=LikelySomatic,Number=0,Type=Flag,Description="VarDict Likely Somatic">')
        vcfheader_filter_info_filter.append('##INFO=<ID=LikelyLOH,Number=0,Type=Flag,Description="VarDict Likely LOH">')
        vcfheader_filter_info_filter.append('##INFO=<ID=StrongLOH,Number=0,Type=Flag,Description="VarDict Strong LOH">')
        vcfheader_filter_info_filter.append('##INFO=<ID=AFDiff,Number=0,Type=Flag,Description="VarDict AF Diff">')
        vcfheader_filter_info_filter.append('##INFO=<ID=Deletion,Number=0,Type=Flag,Description="VarDict Deletion">')
        vcfheader_filter_info_filter.append('##INFO=<ID=SampleSpecific,Number=0,Type=Flag,Description="VarDict SampleSpecific">')
        #vcfheader_filter_info_filter.append('##INFO=<ID=SOR,Number=1,Type=String,Description="VarDict stuff, really a float, but could be Inf">')
        #vcfheader_filter_info_filter.append('##INFO=<ID=SSF,Number=1,Type=Float,Description="VarDict stuff">')

        #vcfheader_filter_info_filter.append('##FORMAT=<ID=QSTD,Number=1,Type=String,Description="Quality score STD in reads. Use string because it can be NA.">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Oddratio">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=PMEAN,Number=1,Type=Float,Description="Mean position in reads">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=BIAS,Number=.,Type=Integer,Description="Strand Bias Info">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=PSTD,Number=1,Type=String,Description="Position STD in reads. Use string because can be NA.">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=SN,Number=1,Type=Float,Description="Signal to noise">')
        #vcfheader_filter_info_filter.append('##FORMAT=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">')
        vcfheader_filter_info_filter.append('##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">')
        

    # Print headers:
    vcfheader_filter_info_filter.sort()

    # Write the headers:
    if args.call_method == 'VarDict':
        
        snpout.write( vcffileformat+'\n' )
        [ snpout.write(i+'\n') for i in vcfheader_misc ]
        [ snpout.write(i+'\n') for i in vcfheader_filter_info_filter ]
        snpout.write( vcfheader_main+'\n' )
        
        indelout.write( vcffileformat+'\n' )
        [ indelout.write(i+'\n') for i in vcfheader_misc ]
        [ indelout.write(i+'\n') for i in vcfheader_filter_info_filter ]
        indelout.write( vcfheader_main+'\n' )
        
    else:
        vcfout.write( vcffileformat+'\n' )
        [ vcfout.write(i+'\n') for i in vcfheader_misc ]
        [ vcfout.write(i+'\n') for i in vcfheader_filter_info_filter ]
        vcfout.write( vcfheader_main+'\n' )

    # Continue:
    line_i = vcf.readline().rstrip()




### It will be a single file there is a single file, otherwise it goes chr1, chr2, ...., or in the list of files input as ordered:
for chr_i_vcf in right_files:
    
    with open_my_vcf(chr_i_vcf) as vcf:
        
        line_i = vcf.readline().rstrip()
        
        # Skip headers from now on:
        while line_i.startswith('#'):
            line_i = vcf.readline().rstrip()
            
        
        # Doing the work here:
        while line_i:
            
            
            # -------------------- VarScan2 -------------------- #
            if args.call_method == 'VarScan2':
                
                item = line_i.split('\t')
                
                # Replace the wrong "G/A" with the correct "G,A" in ALT column:
                item[idx_alt] = item[idx_alt].replace('/', ',')
                
                # Replace the "SOMATIC" tag in INFO with "VarScan2" tag:
                item[idx_info] = item[idx_info].replace('SOMATIC', args.call_method)
                
                # vcf-validator is not going to accept multiple sequences in the REF, as is the case in VarScan2's indel output:
                item[idx_ref] = re.sub(r'[^\w].*$', '', item[idx_ref])
                
                # Get rid of non-compliant characters in the ALT column:
                item[idx_alt] = re.sub(r'[^\w,.]', '', item[idx_alt])
                
                # Eliminate dupliate entries in ALT:
                item[idx_alt] = re.sub(r'(\w+),\1', r'\1', item[idx_alt])
                
                # Eliminate ALT entries when it matches with the REF column, to address vcf-validator complaints:
                if ',' in item[idx_alt]:
                    alt_item = item[idx_alt].split(',')
                    
                    if item[idx_ref] in alt_item:
                        
                        bad_idx = alt_item.index(item[idx_ref])
                        
                        alt_item.pop(bad_idx)
                        
                        item[idx_alt] = ','.join(alt_item)
                    
                    # To fix this vcf-validator complaints:
                    # Could not parse the allele(s) [GTC], first base does not match the reference
                    for n1,alt_i in enumerate(alt_item[1::]):
                        if not alt_i.startswith(item[idx_ref]):
                            
                            alt_item.pop(n1+1)
                            item[idx_alt] = ','.join(alt_item)
                
                
                # Reform the line:
                line_i = '\t'.join(item)
                
                # VarScan2 output a line with REF allele as "M". GATK CombineVariants complain about that.
                if not re.search(r'[^GCTAU]', item[idx_ref], re.I):
                    vcfout.write(line_i+'\n')
            
            
            
            # -------------------- SomaticSniper -------------------- #
            elif args.call_method == 'SomaticSniper':
                
                # Print "SomaticSniper" into the INFO field if it is called so, otherwise never mind.
                item = line_i.split('\t')
                
                tumor_item = item[idxT].split(':')
                field_item = item[idx_format].split(':')
                
                try:
                    ss_idx = field_item.index('SS')
                    ss_status = tumor_item[ss_idx]
                except ValueError:
                    ss_status = 0
                
                if ss_status == '2':
                    item = line_i.split('\t')
                    item[idx_info] = args.call_method
                    line_i = '\t'.join(item)
                
                vcfout.write(line_i+'\n')
            
            
            
            # -------------------- JointSNVMix2 --------------- #
            elif args.call_method == 'JointSNVMix2':
                
                # Print "JointSNVMix2" into the INFO field.
                item = line_i.split('\t')
                item[idx_info] = ';'.join( (args.call_method, item[idx_info]) )
                line_i = '\t'.join(item)
                
                vcfout.write(line_i+'\n')
                
                
            
            # -------------------- VarDict -------------------- #
            elif args.call_method == 'VarDict':
                
                vcfcall = genome.Vcf_line( line_i )
                item = line_i.split('\t')
                
                if vcfcall.refbase != vcfcall.altbase:
                    
                    if args.filter_version == 'v1':
                        
                        second_vardict_filter = False
                    
                    elif args.filter_version == 'v2':
                        second_vardict_filter = vcfcall.filters == 'P0.05' and float(vcfcall.get_info_value('SSF') ) < 0.1 
                        
                    elif args.filter_version == 'somatic':
                        
                        vardict_filters = vcfcall.filters.split(';')
                        
                        bad_filters = ('d7' in vardict_filters) or ('DIFF0.2' in vardict_filters) or ('LongAT' in vardict_filters) or ('MAF0.05' in vardict_filters) or ('MSI6' in vardict_filters) or ('NM4' in vardict_filters) or ('pSTD' in vardict_filters) or ('SN1.5' in vardict_filters) or ('P0.05' in vardict_filters and float(vcfcall.get_info_value('SSF') ) >= 0.15 ) or ('v4' in vardict_filters and int(vcfcall.get_sample_value('VD', 0))<3)
                        
                        no_bad_filter = not bad_filters
                        
                        filter_fail_times = len(vardict_filters)
                        
                        # No more than 2 fails for the filter:
                        second_vardict_filter = no_bad_filter and filter_fail_times<=2
                        
                    else:
                        raise Exception( '{} is a wrong VarDict Filter Argument'.format(args.variant_type) )
                        
                        
                    # Add the "VarDict" tag if the site is labeled (StrongSomatic or LikelySomatic) and PASS:
                    if ('StrongSomatic' in item[idx_info] or 'LikelySomatic' in item[idx_info] ) and \
                    (vcfcall.filters == 'PASS' or second_vardict_filter):
                    
                        item[idx_info] = ';'.join(( args.call_method, item[idx_info] ))
    
                    # Fix VarDict Typo, should be fixed by VarDict soon
                    elif re.match(r'Germlin;', item[idx_info]):
                        item[idx_info] = item[idx_info].replace('Germlin;', 'Germline;')
                        
                    
                    # In the REF field, non-GCTA characters should be changed to N to fit the VCF standard:
                    item[idx_ref] = re.sub( r'[^GCTA]', 'N', item[idx_ref], flags=re.I )
                    
                    
                    ## To be consistent with other tools, Combine AD:RD into DP4.
                    # VarDict puts Tumor first and Normal next
                    format_field = item[idx_format].split(':')
                    tumor_sample  = item[9].split(':')
                    normal_sample = item[10].split(':')
    
                    idx_rd = format_field.index('RD')
                    
                    normal_dp4 = normal_sample.pop(idx_rd)
                    tumor_dp4  = tumor_sample.pop(idx_rd)
                    format_field.pop(idx_rd)
                    
                    
                    ## Be careful not to find the index and then pop them..... index changes each time you pop something.
                    idx_ad = format_field.index('AD')
                    
                    normal_dp4 = normal_dp4 + ',' + normal_sample.pop(idx_ad)
                    tumor_dp4  = tumor_dp4  + ',' + tumor_sample.pop(idx_ad)
                    format_field.pop(idx_ad)
    
    
                    # Re-format the strings:
                    format_field.append('DP4')
                    normal_sample.append(normal_dp4)
                    tumor_sample.append(tumor_dp4)
                    
                    normal_sample = ':'.join(normal_sample)
                    tumor_sample  = ':'.join(tumor_sample)
                    item[idx_format] = ':'.join(format_field)
    
                    line_i = '\t'.join(( '\t'.join(item[0:idx_format+1]), normal_sample, tumor_sample ))
                    
                    
                    # Write to snp and indel into different files:
                    if 'TYPE=SNV' in item[idx_info]:
                        snpout.write(line_i+'\n')
                        
                    elif 'TYPE=Deletion' in item[idx_info] or 'TYPE=Insertion' in item[idx_info]:
                        indelout.write(line_i+'\n')
            


            # -------------------- VarScan2_single -------------------- #
            elif args.call_method == 'VarScan2_single':
                
                item = line_i.split('\t')
                
                # Replace the wrong "G/A" with the correct "G,A" in ALT column:
                item[idx_alt] = item[idx_alt].replace('/', ',')
                
                # Replace the "SOMATIC" tag in INFO with "VarScan2" tag:
                item[idx_info] = ';'.join( (args.call_method, item[idx_info]) )
                
                # vcf-validator is not going to accept multiple sequences in the REF, as is the case in VarScan2's indel output:
                item[idx_ref] = re.sub(r'[^\w].*$', '', item[idx_ref])
                
                # Get rid of non-compliant characters in the ALT column:
                item[idx_alt] = re.sub(r'[^\w,.]', '', item[idx_alt])
                
                # Eliminate dupliate entries in ALT:
                item[idx_alt] = re.sub(r'(\w+),\1', r'\1', item[idx_alt])
                
                # Eliminate ALT entries when it matches with the REF column, to address vcf-validator complaints:
                if ',' in item[idx_alt]:
                    alt_item = item[idx_alt].split(',')
                    
                    if item[idx_ref] in alt_item:
                        
                        bad_idx = alt_item.index(item[idx_ref])
                        
                        alt_item.pop(bad_idx)
                        
                        item[idx_alt] = ','.join(alt_item)
                    
                    # To fix this vcf-validator complaints:
                    # Could not parse the allele(s) [GTC], first base does not match the reference
                    for n1,alt_i in enumerate(alt_item[1::]):
                        if not alt_i.startswith(item[idx_ref]):
                            
                            alt_item.pop(n1+1)
                            item[idx_alt] = ','.join(alt_item)
                
                
                # Reform the line:
                line_i = '\t'.join(item)
                
                # VarScan2 output a line with REF allele as "M". GATK CombineVariants complain about that.
                if not re.search(r'[^GCTAU]', item[idx_ref], re.I):
                    vcfout.write(line_i+'\n')
            

            
            # Next line:
            line_i = vcf.readline().rstrip()



# Close files and exit:
if args.call_method == 'VarDict':
    snpout.close()
    indelout.close()
else:
    vcfout.close()
