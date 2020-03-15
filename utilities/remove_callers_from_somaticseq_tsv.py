#!/usr/bin/env python3

import argparse, gzip
from math import nan, isnan

all_possible_callers = ('if_MuTect', 'if_VarScan2', 'if_JointSNVMix2', 'if_SomaticSniper', 'if_VarDict', 'MuSE_Tier', 'if_LoFreq', 'if_Scalpel', 'if_Strelka', 'if_TNscope', 'if_Platypus')


parser = argparse.ArgumentParser(description='In SomaticSeq TSV files, replace certain callers with nan and remove lines where they are only called by these. To mimic a TSV where only a subset of the callers were used.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',   '--infile',           type=str, help='input file', required=True)
parser.add_argument('-outfile',  '--outfile',          type=str, help='input file', required=True)
parser.add_argument('-subtract', '--subtract-callers', type=str, nargs='+', help='columns to make nan', required=True, choices=all_possible_callers)

args = parser.parse_args()




NA = NaN = NA = nan

for caller_i in args.subtract_callers:
    assert caller_i in all_possible_callers





def open_textfile(file_name):
    if file_name.lower().endswith('.gz'):
        return gzip.open(file_name, 'rt')
    else:
        return open(file_name)



def items_to_make_nan(callers_to_subtract):
    
    out_items = []
    
    for caller_i in callers_to_subtract:
        if caller_i == 'if_MuTect':
            out_items.append('M2_NLOD')
            out_items.append('M2_TLOD')
            out_items.append('M2_STR')
            out_items.append('M2_ECNT')
        
        elif caller_i == 'if_JointSNVMix2':
            out_items.append('SNVMix2_Score')
        
        elif caller_i == 'if_SomaticSniper':
            out_items.append('Sniper_Score')
        
        elif caller_i == 'if_VarDict':
            out_items.append('VarDict_Score')
            out_items.append('MSI')
            out_items.append('MSILEN')
            out_items.append('SHIFT3')

        elif caller_i == 'if_Strelka':
            out_items.append('Strelka_Score')
            out_items.append('Strelka_QSS')
            out_items.append('Strelka_TQSS')

    return out_items



with open_textfile(args.infile) as infile, open(args.outfile, 'w') as outfile:
    
    line_in = infile.readline().rstrip()
    item_in = line_in.split('\t')
    
    out_indices       = [ item_in.index(i) for i in args.subtract_callers ]
    remaining_indices = [ item_in.index(i) for i in all_possible_callers if i not in args.subtract_callers ]
    extra_nan_items   = items_to_make_nan( args.subtract_callers)
    extra_nan_indices = [ item_in.index(i) for i in extra_nan_items ]
    
    outfile.write( line_in + '\n' )
    
    line_in = infile.readline().rstrip()
    
    while line_in:
        
        item_in = line_in.split('\t')
                
        other_callers = 0
        for other_i in remaining_indices:
            classification_i = item_in[other_i]
            classification_i = eval( classification_i )
            if not isnan( classification_i ):
                other_callers += classification_i


        if other_callers > 0:
            for out_i in out_indices + extra_nan_indices:
                item_in[out_i] = 'nan'
                            
            line_out = '\t'.join( item_in )
            
            outfile.write( line_out + '\n' )

        line_in = infile.readline().rstrip()
