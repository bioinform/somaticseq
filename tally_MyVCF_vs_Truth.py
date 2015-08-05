#!/usr/bin/env python3
# Tally the same somatic mutation call from MuTect and VarScan2 (or any two VCF files)
# Take readings from refined_excelstyle_nonsynonymous.csv
# 3/5/2014
# Assumes Tumor and then Normal in MuTect
# Output will be Tumor and then Normal. 


import sys, argparse
import regex as re
import gzip
from os import sep


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-myvcf', '--my-vcf', type=str, help='My VCF', required=True, default=None)
parser.add_argument('-truth', '--truth-vcf', type=str, help='Truth', required=True, default=None)
parser.add_argument('-outfile', '--output-file', type=str, help='Out File Label', required=True, default=None)

parser.add_argument('-myid', '--my-id', type=str, help='FILE1 only tag', required=False, default='FalsePositive')
parser.add_argument('-truthid', '--truth-id', type=str, help='FILE2 only tag', required=False, default='FalseNegative')
parser.add_argument('-bothid', '--both-id', type=str, help='BOTH ID tag', required=False, default='Correct')


args = parser.parse_args()

# VCF File's index locations
idx_chrom,idx_pos,idx_id,idx_ref,idx_alt,idx_qual,idx_filter,idx_info,idx_format,idxN,idxT = 0,1,2,3,4,5,6,7,8,9,10


# The regular expression pattern for "chrXX 1234567" in both VarScan2 Output and VCF files:
pattern_chr_position = re.compile(r'^(?:chr)?(?:[0-9]+|[XYM]T?)\t[0-9]+\b')

###
# Define which chromosome coordinate is ahead:
chrom_sequence = [str(i) for i in range(1,23)]
chrom_sequence.append('X')
chrom_sequence.append('Y')
chrom_sequence.append('M')

def whoisbehind(coord_0, coord_1):
    '''coord_0 and coord_1 are two strings, specifying the chromosome, a (typically) tab, and then the location.'''
    
    if coord_0=='' and coord_1=='':
        return 10
        
    elif coord_1 == '':
        return 0
        
    elif coord_0 == '':
        return 1
    
    else:
    
        if isinstance(coord_0, str):
            chrom0, position0 = coord_0.split()
            chrom1, position1 = coord_1.split()
        elif isinstance(coord_0, list):
            chrom0, position0 = coord_0[0], coord_0[1]
            chrom1, position1 = coord_1[0], coord_1[1]
        
        chrom0 = chrom0.lstrip('chr').rstrip('T')
        chrom1 = chrom1.lstrip('chr').rstrip('T')
        
        if chrom_sequence.index(chrom0) < chrom_sequence.index(chrom1):
            return 0   # 1st coordinate is ahead
            
        elif chrom_sequence.index(chrom0) > chrom_sequence.index(chrom1):
            return 1   # 1st coordinate is ahead
        
        # Must be in the same chromosome
        else:
            
            position0 = int(position0)
            position1 = int(position1)
            
            if position0 < position1:
                return 0
            
            elif position0 > position1:
                return 1
            
            # Same chromosome, same position, then same coordinate:
            else:
                return 10



def only_care(nth_col, itsays, string_input):
    '''itsays is a regex pattern that matches string_input in the nth_col. Returns True if it matches.  Returns False it it does not.
    IGNORE CASE FOR NOW.'''
    
    itmatches = re.search(itsays, string_input.split()[nth_col], re.I)
    
    if itmatches:
        return True
    else:
        return False




def catch_up(line_1, line_2, file_1, file_2, output_vcf, id_1, id_2, id_12 ):
    
    id_1, id_2, id_12 = id_1, id_2, id_12
    
    # Read the first file until it hits data
    coord_1 = re.search( pattern_chr_position, line_1 )
        
    if coord_1:
        coord_1 = coord_1.group()
    else:
        coord_1 = ''

    # Read the second file until it hits data
    coord_2 = re.search( pattern_chr_position, line_2 )
    
    if coord_2:
        coord_2 = coord_2.group()
    else:
        coord_2 = ''
    
    is_behind = whoisbehind( coord_1, coord_2 )
    
    # As long as the coordinates are not the same, and both files are not finished:
    while is_behind != 10:
        
        # If 1st VCF is behind:
        if is_behind == 0:
            
            item_1 = line_1.rstrip('\n').split('\t')
            
            # Write, unless...
            if item_1[idx_filter] != 'PrintEmALL':
                
                #item_1[idx_id] = id_1
                id_item = item_1[idx_id].split(';')
                id_item.append( id_1 )
                item_1[idx_id] = ';'.join( id_item )
                item_1[idx_id] = re.sub(r'^\.;', '', item_1[idx_id])
                
                line_1 = '\t'.join( item_1 )
                
                output_vcf.write( line_1  + '\n')
            
            line_1 = file_1.readline()
            next_coord = re.search( pattern_chr_position, line_1 )
            
            if next_coord:
                coord_1 = next_coord.group()
            else:
                coord_1 = ''
        
        # If 2nd VCF is behind:    
        elif is_behind == 1:
            
            item_2 = line_2.rstrip('\n').split('\t')
            
            # Write, unless...
            #if item_2[idx_filter] != 'PrintEmALL':
                
            #IF 
            #item_2[idx_id] = id_2
            id_item = item_2[idx_id].split(';')
            id_item.append( id_2 )
            item_2[idx_id] = ';'.join( id_item )
            item_2[idx_id] = re.sub(r'^\.;', '', item_2[idx_id])

            
            line_2 = '\t'.join( item_2 )
            
            output_vcf.write( line_2 + '\n' )
            ## FI
            
            line_2 = file_2.readline()
            next_coord = re.search( pattern_chr_position, line_2 )
            
            if next_coord:
                coord_2 = next_coord.group()
            else:
                coord_2 = ''
        
        is_behind = whoisbehind( coord_1, coord_2 )
    
    
    # Returns the value of the function:
    if (coord_1=='') and (coord_2==''):
        result = 42
    else:
        
        item_1 = line_1.rstrip('\n').split('\t')
        item_2 = line_2.rstrip('\n').split('\t')
        
        #item_1[idx_id] = id_12
        id_item = item_1[idx_id].split(';')
        id_item.append( id_12 )
        item_1[idx_id] = ';'.join( id_item )
        item_1[idx_id] = re.sub(r'^.;', '', item_1[idx_id])
        
        
        line_1 = '\t'.join( item_1 )
        
        
        output_vcf.write( line_1  + '\n')
        
        result = (line_1, line_2,)
    
    return result
    



def open_my_vcf(file_name):
    
    # See if the input file is a .gz file:
    if file_name.lower().endswith('.gz'):
        return gzip.open(file_name, 'rt')
        
    else:
        return open(file_name)





# Read both files line by line:
with open_my_vcf(args.my_vcf) as myvcf, \
     open_my_vcf(args.truth_vcf)  as truthvcf, \
     open(args.output_file,'w') as vcfout:

    
    # For both VCF files, read until it hits data:
    double_palm_collector = []
    my_line  = myvcf.readline()
    while my_line.startswith('#'):
        
        if re.match(r'##fileformat=', my_line):
            header0 = my_line
        
        if re.match(r'##(INFO|FORMAT)', my_line):
            double_palm_collector.append( my_line )
            
        if re.match(r'#CH', my_line):
            headermain = my_line
        
        
        my_line = myvcf.readline()
    
    
    truth_line  = truthvcf.readline()
    while truth_line.startswith('#'):
        
        if re.match(r'##(INFO|FORMAT)', truth_line):
            
            double_palm_collector.append( truth_line )
                
        truth_line  = truthvcf.readline()
        
    
    double_palm_collector.sort()
    
    try:
        vcfout.write(header0)
    except NameError:
        pass
    
    try:
        [ vcfout.write(i) for i in double_palm_collector ]
    except NameError:
        pass
        
    try:
        vcfout.write(headermain)
    except NameError:
        pass
    
    # Start comparing:
    sidebyside = 0
        
    while sidebyside != 42:
        sidebyside = catch_up(my_line, truth_line, myvcf, truthvcf, vcfout, args.my_id, args.truth_id, args.both_id)
        my_line = myvcf.readline()
        truth_line  = truthvcf.readline()
