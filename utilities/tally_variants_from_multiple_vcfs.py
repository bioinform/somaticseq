#!/usr/bin/env python3

import argparse, gzip, re, sys, os, math
import genomicFileHandler.genomic_file_handlers as genome

def extract_snpEff(vcf_line):
    
    annGenes = []
    annAAs   = []
    annTxns  = []

    vcf_obj    = genome.Vcf_line(vcf_line)
    snpeff_ann = vcf_obj.get_info_value('ANN')
    
    if snpeff_ann:
        ann_items = snpeff_ann.split(',')
        
        for ann_i in ann_items:
        
            ann_item   = ann_i.split('|')
            gene_i     = ann_item[3]
            feature_i  = ann_item[6]
            ntchange_i = ann_item[9]
            aaChange_i = ann_item[10]
            
            if gene_i and aaChange_i:
                
                # Only do non-syn variants
                aa = re.search(r'p\.([a-zA-Z]+)[0-9]+([a-zA-Z]+)', aaChange_i)
                
                if aa and (aa.groups()[0] != aa.groups()[1]):
                    
                    annGenes.append( gene_i )
                    annAAs.append( aaChange_i )
                    annTxns.append( feature_i )

    return annGenes, annAAs, annTxns



def extract_dbsnp_cosmic(vcf_line):
    ids = []
    id_items = vcf_line.split('\t')[2]
    for item_i in re.split(r'[,;]', id_items):
        if re.match(r'rs[0-9]+', item_i):
            ids.append( item_i )
        elif re.match(r'COSM[0-9]+', item_i):
            ids.append( item_i )
            
    return ids



def vcfs2variants(vcf_files, sample_names, SampleColumnNumber, key='VAF'):
    
    assert len(vcf_files) == len(sample_names)
    
    variantDict = {}
    
    for vcf_file_i, sample_name_i in zip(vcf_files, sample_names):
        
        with genome.open_textfile(vcf_file_i) as vcf:
            line_i = vcf.readline().rstrip()
            while line_i.startswith('#'):
                line_i = vcf.readline().rstrip()
                
            while line_i:
                
                vcf_obj = genome.Vcf_line(line_i)
                item    = line_i.split('\t')

                contig_i   = item[0]
                pos_i      = int(item[1])
                refbase    = item[3]
                altbase    = item[4]
                ID_field   = item[2].split(';')
                filter_i   = item[6].split(';')
                
                genes, amino_acid_changes, txn_ids = extract_snpEff(line_i)
                dbsnp_cosmic_ids = extract_dbsnp_cosmic(line_i)
                
                variant_id = (contig_i, pos_i, refbase, altbase,)
                vaf_i      = vcf_obj.get_sample_value(key, SampleColumnNumber)

                if vaf_i:
                    try:
                        vaf_i = float(vaf_i)
                    except ValueError:
                        vaf_i = math.nan
                else:
                    vaf_i = math.nan

                if variant_id not in variantDict:
                    variantDict[ variant_id ] = {}
                    variantDict[ variant_id ]['GENES']      = genes
                    variantDict[ variant_id ]['AAChange']   = amino_acid_changes
                    variantDict[ variant_id ]['TRANSCRIPT'] = txn_ids
                    variantDict[ variant_id ]['DATABASE']   = dbsnp_cosmic_ids

                variantDict[ variant_id ][ sample_name_i ] = {'FILTER': filter_i, 'VAF': vaf_i,}
                line_i = vcf.readline().rstrip()

    return variantDict



def print_variantDict(variantDict, sample_names, filter_labels=['PASS',], min_number=0):
    
    line_out = 'CHROM\tPOS\tREF\tALT\tID\tAAChange\tNUM\t' + '\t'.join(sample_names)
    print( line_out )
    
    for variant_i in variantDict:
        
        if variantDict[variant_i]['DATABASE']:
            id_column = ','.join( variantDict[variant_i]['DATABASE'] )
        else:
            id_column = '.'
        
        if variantDict[variant_i]['GENES']:
            aa_column = '{}:{}'.format( variantDict[variant_i]['GENES'][0], ','.join( set(variantDict[variant_i]['AAChange'] )) )
        else:
            aa_column = '.'
        
        num_intersected_filters = 0
        data_text = []
        for sample_i in sample_names:
                                    
            if sample_i not in variantDict[variant_i]:
                text_i = 'NONE'
                
            else:
                if set(variantDict[variant_i][sample_i]['FILTER']) & set(filter_labels):
                    num_intersected_filters += 1

                text_i = '%s/%g' % ( ','.join(variantDict[variant_i][sample_i]['FILTER']), variantDict[variant_i][sample_i]['VAF'] )

            data_text.append(text_i)
        
        if num_intersected_filters >= min_number:
            
            line_out = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( variant_i[0], variant_i[1], variant_i[2], variant_i[3], id_column, aa_column, num_intersected_filters, '\t'.join(data_text) )
            
            print( line_out )

    return 0




def run():
    parser = argparse.ArgumentParser(description='Tally common variants in multiple VCF files, and also print out their FILTER label and VAF', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-vcfs',    '--vcf-files',     nargs='+',  type=str,  help='multiple vcf files',      required=True)
    parser.add_argument('-names',   '--sample-names',  nargs='+',  type=str,  help='names for the vcf files', required=True)
    parser.add_argument('-filters', '--filter-labels', nargs='+',  type=str,  help='Filter labels to count', default=['PASS',] )
    parser.add_argument('-key',     '--vaf-key',                   type=str,  help='What is the key for VAF', required=True)
    parser.add_argument('-column',  '--sample-column',             type=int,  help='Which sample column to get the VAF. If it is the first column, it is 0. If the 2nd column, it is 1.', required=True)
    parser.add_argument('-min',     '--minimum-samples',           type=int,  help='Print out only if at least this number of vcf files have the variant.', default=0)
    
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    args = run()
    variant_dict = vcfs2variants(args.vcf_files, args.sample_names, args.sample_column, args.vaf_key)
    print_variantDict(variant_dict, args.sample_names, args.filter_labels, args.minimum_samples)
