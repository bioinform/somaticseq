#!/usr/bin/env python3

import argparse, re, os, math, pysam, uuid
import multiprocessing
import somaticseq.genomicFileHandler.genomic_file_handlers as genome
import somaticseq.genomicFileHandler.read_info_extractor as read_info_extractor
import somaticseq.utilities.split_Bed_into_equal_regions as split_regions
import somaticseq.vcfModifier.vcfIntersector as vcfIntersector
import tempfile
from warnings import warn

def intersect_multiple_vcf_files(inVcfList_and_BedFile):
    
    vcf_files, bed_file = inVcfList_and_BedFile
    dirname = tempfile.gettempdir()
    
    out_files = []
    for vcf_i in vcf_files:
        outfile_i = os.path.join(dirname, uuid.uuid4().hex+'.vcf')
        vcfIntersector.bed_intersector(vcf_i, outfile_i, bed_file)
        out_files.append( outfile_i )

    return out_files


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
        elif re.match(r'(COSM|COSV)[0-9]+', item_i):
            ids.append( item_i )
            
    return ids



def vaf_from_bam(bam, my_coordinate, ref_base, first_alt, min_mq=1, ):

    '''
    bam is the opened file handle of bam file
    my_coordiate is a list or tuple of 0-based (contig, position)
    Returns: number of variant calls, reference calls, other calls, and total calls
    '''
    
    indel_length = len(first_alt) - len(ref_base)
    reads = bam.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] )

    dp          = 0
    var_calls   = 0
    ref_calls   = 0
    other_calls = 0
    for read_i in reads:
        if (not read_i.is_unmapped) and read_info_extractor.dedup_test(read_i) and read_i.mapping_quality >= min_mq:
            
            dp += 1
            
            code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = read_info_extractor.position_of_aligned_read(read_i, my_coordinate[1]-1 )

            # Reference calls:
            if code_i == 1 and base_call_i == ref_base[0]:
                ref_calls += 1
            
            # Alternate calls:
            # SNV, or Deletion, or Insertion where I do not check for matching indel length
            elif (indel_length == 0 and code_i == 1 and base_call_i == first_alt) or \
                 (indel_length < 0  and code_i == 2 and indel_length == indel_length_i) or \
                 (indel_length > 0  and code_i == 3):

                var_calls += 1
            
            # Inconsistent read or 2nd alternate calls:
            else:
                other_calls += 1
    
    return var_calls, ref_calls, other_calls, dp
    


def vcfs2variants(vcf_files, bam_files, sample_names):
    
    assert len(vcf_files) == len(sample_names) == len(bam_files)
    
    variantDict = {}
    i = 0
    for vcf_file_i, bam_file_i, sample_name_i in zip(vcf_files, bam_files, sample_names):
        
        with genome.open_textfile(vcf_file_i) as vcf, pysam.AlignmentFile(bam_file_i) as bam:
                        
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
                
                vdp, rdp, odp, totaldp = vaf_from_bam(bam, (contig_i, pos_i), refbase, altbase, 1)

                try:
                    vaf_i = vdp / totaldp
                except ZeroDivisionError:
                    vaf_i = math.nan

                if variant_id not in variantDict:
                    variantDict[ variant_id ] = {}
                    variantDict[ variant_id ]['GENES']      = genes
                    variantDict[ variant_id ]['AAChange']   = amino_acid_changes
                    variantDict[ variant_id ]['TRANSCRIPT'] = txn_ids
                    variantDict[ variant_id ]['DATABASE']   = dbsnp_cosmic_ids

                variantDict[ variant_id ][ sample_name_i ] = {'FILTER': filter_i, 'VAF': vaf_i, 'VDP': vdp, 'DP': totaldp}
                
                line_i = vcf.readline().rstrip()

        i += 1

    return variantDict



def fills_missing_vafs(variantDict, bam_files, sample_names):
    
    assert len(sample_names) == len(bam_files)

    bamDict = {}
    for bam_i, sample_i in zip(bam_files, sample_names):
        bamDict[ sample_i ] = pysam.AlignmentFile(bam_i)

    for variant_i in variantDict:
        
        for sample_i in sample_names:
            
            if sample_i not in variantDict[ variant_i ]:
                
                vdp, rdp, odp, totaldp = vaf_from_bam(bamDict[sample_i], (variant_i[0], variant_i[1]), variant_i[2], variant_i[3], 1)
                
                try:
                    vaf_i = vdp / totaldp
                except ZeroDivisionError:
                    vaf_i = math.nan

                variantDict[ variant_i ][ sample_i ] = {'FILTER': ['NONE',], 'VAF': vaf_i, 'VDP': vdp, 'DP': totaldp}

    for sample_i in bamDict:
        bamDict[sample_i].close()

    return variantDict



def make_variant_dict(inputListofLists):
    vcf_files, bam_files, sample_names = inputListofLists
    varDictWithMissings = vcfs2variants(vcf_files, bam_files, sample_names)
    completeVariantDict = fills_missing_vafs(varDictWithMissings, bam_files, sample_names)
    
    return completeVariantDict



def make_variant_dict_parallel(vcf_files, bam_files, sample_names, bed_region, nthreads):

    dirname         = tempfile.gettempdir() #os.curdir
    partial_regions = split_regions.split(bed_region, os.path.join(dirname, uuid.uuid4().hex+'.bed'), nthreads)

    pool = multiprocessing.Pool(nthreads)
    
    map_args              = map( lambda bed_i: (vcf_files, bed_i), partial_regions )
    partitioned_vcf_files = pool.map_async(intersect_multiple_vcf_files, map_args).get()
    
    map_args             = map( lambda partial_vcf_files: (partial_vcf_files, bam_files, sample_names), partitioned_vcf_files )
    variant_dictionaries = pool.map_async(make_variant_dict, map_args).get()

    pool.close()


    # partitioned_vcf_files is a list of a list of vcf files
    vcfs_to_be_deleted = []
    for list_of_vcfs in partitioned_vcf_files:
        vcfs_to_be_deleted = vcfs_to_be_deleted + list_of_vcfs

    for file_i in partial_regions + vcfs_to_be_deleted:
        os.remove( file_i )

    return variant_dictionaries


def print_variantDict(variantDict, sample_names, filter_labels=['PASS',], min_number=1, print_header=True):
    
    if print_header:
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
                                    
            if set(variantDict[variant_i][sample_i]['FILTER']) & set(filter_labels):
                
                num_intersected_filters += 1

            text_i = '%s:%i/%i=%g' % ( ','.join(variantDict[variant_i][sample_i]['FILTER']), variantDict[variant_i][sample_i]['VDP'], variantDict[variant_i][sample_i]['DP'], variantDict[variant_i][sample_i]['VAF'] )
            data_text.append(text_i)
        
        if num_intersected_filters >= min_number:
            
            line_out = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( variant_i[0], variant_i[1], variant_i[2], variant_i[3], id_column, aa_column, num_intersected_filters, '\t'.join(data_text) )
            
            print( line_out )

    return 0




def run():
    parser = argparse.ArgumentParser(description='Tally common variants in multiple VCF files, and also print out their FILTER label and VAF', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-vcfs',    '--vcf-files',     nargs='+', type=str,  help='multiple vcf files',      required=True)
    parser.add_argument('-bams',    '--bam-files',     nargs='+', type=str,  help='multiple bam files',      required=True)
    parser.add_argument('-names',   '--sample-names',  nargs='+', type=str,  help='names for the vcf files', required=True)
    parser.add_argument('-filters', '--filter-labels', nargs='+', type=str,  help='Filter labels to count', default=['PASS',] )
    parser.add_argument('-bed',     '--bed-inclusion',            type=str,  help='Bed file to include.')
    parser.add_argument('-nt',      '--num-threads',              type=int,  help="threads", default=1)
    parser.add_argument('-min',     '--minimum-samples',          type=int,  help='Print out only if at least this number of vcf files have the variant.', default=1)
    
    args = parser.parse_args()
    
    return args




if __name__ == '__main__':
    
    args = run()
    
    if args.num_threads > 1 and args.bed_inclusion:
        
        partitioned_variant_dicts = make_variant_dict_parallel(args.vcf_files, args.bam_files, args.sample_names, args.bed_inclusion, args.num_threads)
        
        for n, dict_i in enumerate(partitioned_variant_dicts):
            if n == 0:
                printHeader = True
            else:
                printHeader = False
                
            print_variantDict(dict_i, args.sample_names, args.filter_labels, args.minimum_samples, printHeader)
        
    else:
        
        if args.num_threads > 1:
            warn('This module is unable parallelize this task without bed file input.')
        
        variant_dict = make_variant_dict( (args.vcf_files, args.bam_files, args.sample_names) )
        print_variantDict(variant_dict, args.sample_names, args.filter_labels, args.minimum_samples)
