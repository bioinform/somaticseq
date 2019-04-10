#!/usr/bin/env python3

import sys, os, argparse, pysam, re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append( MY_DIR )

import refGene_reader


def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-genes',  '--genes', type=str, nargs='*', help='Input file',  required=True)

    args     = parser.parse_args()
    genes_in = args.genes

    if len(genes_in) == 1 and os.path.isfile(genes_in[0]):
        with open( genes_in[0] ) as fn:
            genes_in = []
            for line_i in fn:
                genes_in.append( line_i.rstrip() )
                
    return genes_in
        





def genes_to_beds(genes_in):
    regions = {}
    for gene_i in genes_in:
        
        gene_data_i = refGene_reader.refGene_data[gene_i]
        Gene_i = refGene_reader.Gene( gene_data_i )
        
        for txn_i in Gene_i.transcripts:
            if txn_i['cds_start_stat'] == txn_i['cds_end_stat'] == 'cmpl':
                
                if txn_i['strand'] == '+':
                    exon_num = 1
                elif txn_i['strand'] == '-':
                    exon_num = txn_i['exon_count']
                    
                for start_i, end_i in zip( txn_i['exon_starts'], txn_i['exon_ends'] ):
                    
                    if txn_i['cds_start'] <= end_i and txn_i['cds_end'] >= start_i:
                        
                        region_i = txn_i['chrom'], max(start_i, txn_i['cds_start']), min(end_i, txn_i['cds_end'])
                        
                        #print( '\t'.join(region_i), '{}:{}:{}'.format(txn_i['gene'], txn_i['name'], exon_num), sep='\t')
                        
                        if region_i not in regions:
                            regions[ region_i ] = []
                            
                        regions[ region_i ].append( '{}:{}:{}'.format(txn_i['gene'], txn_i['name'], exon_num) )
                            
                    if txn_i['strand'] == '+':
                        exon_num += 1
                    elif txn_i['strand'] == '-':
                        exon_num -= 1

    return regions
    


def getSequence(coordinates, fastaHandle, padding=0):
    
    seq = fastaHandle.fetch(coordinates[0], coordinates[1]-padding, coordinates[2]+padding)
    return seq


def getVariantSequence(coordinates, refbases, altbases, variant_start_coordinate, fastaHandle):
    
    assert coordinates[1] <= variant_start_coordinate <= coordinates[2]
    
    refseq = fastaHandle.fetch(coordinates[0], coordinates[1], coordinates[2])

    left_cutoff  = variant_start_coordinate - coordinates[1]
    right_cutoff = left_cutoff + len(refbases)
    
    refseq_left  = refseq[:left_cutoff]
    refseq_right = refseq[right_cutoff:]

    variant_seq = refseq_left + altbases + refseq_right
    
    return variant_seq



def gc_content(seq):
    
    gc = re.findall(r'[gcGC]', seq)
    return len(gc)/len(seq)


def longest_homopolymer(sequence):

    '''For a string, count the number of characters that appears in a row.
    E.g., for string "ABBCCCDDDDAAAAAAA", the function returns 1, 2, 3, 4, 7, because there is 1 A, 2 B's, 3 C's, 4 D's, and then 7 A's.
    '''
    counters = []
    previous_base = None

    for current_base in sequence:

        if current_base == previous_base:
            counters[-1] += 1
        else:
            counters.append(1)

        previous_base = current_base

    return max(counters)



def get_homopolymer_length(fasta, coordinate, refbases, altbases):

    covering = 100

    lseq  = fasta.fetch(coordinate[0], coordinate[1]-covering, coordinate[1])
    rseq  = fasta.fetch(coordinate[0], coordinate[1]+len(refbases), coordinate[1]+covering )
    
    ref_seq = lseq + refbases + rseq
    var_seq = lseq + altbases + rseq
    
    # Homopolymer spanning the variant site:
    ref_count_from_right = 1
    ref_count_from_left  = 1
    right_most_ref_base  = refbases[-1]
    left_most_ref_base   = refbases[0]
    
    ref_right = rseq
    ref_left  = lseq + refbases[:-1]
    for i in ref_right:
        if i == right_most_ref_base:
            ref_count_from_right += 1
        else:
            break
            
    for i in ref_left[::-1]:
        if i == right_most_ref_base:
            ref_count_from_right += 1
        else:
            break
    
    ref_right = refbases[1:] + rseq
    ref_left  = lseq
    for i in ref_right:
        if i == left_most_ref_base:
            ref_count_from_left += 1
        else:
            break
            
    for i in ref_left[::-1]:
        if i == left_most_ref_base:
            ref_count_from_left += 1
        else:
            break

    longest_ref_homopolymer = max(ref_count_from_left, ref_count_from_right)
    
    ##
    var_count_from_right = 1
    var_count_from_left  = 1
    right_most_var_base  = altbases[-1]
    left_most_var_base   = altbases[0]
    
    var_right = rseq
    var_left  = lseq + altbases[:-1]
    for i in var_right:
        if i == right_most_var_base:
            var_count_from_right += 1
        else:
            break
            
    for i in var_left[::-1]:
        if i == right_most_var_base:
            var_count_from_right += 1
        else:
            break
    
    var_right = altbases[1:] + rseq
    var_left  = lseq
    for i in var_right:
        if i == left_most_var_base:
            var_count_from_left += 1
        else:
            break
            
    for i in var_left[::-1]:
        if i == left_most_var_base:
            var_count_from_left += 1
        else:
            break

    longest_var_homopolymer = max(var_count_from_left, var_count_from_right)

    return longest_ref_homopolymer, longest_var_homopolymer



if __name__ == '__main__':
    
    genes_in = run()
    regions = genes_to_beds(genes_in)
    
    for region_i in regions:
        print( '\t'.join( [str(i) for i in region_i] ), ','.join(regions[region_i]), sep='\t')


