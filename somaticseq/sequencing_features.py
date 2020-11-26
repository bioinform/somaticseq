#!/usr/bin/env python3

import sys, os, re, pysam
import scipy.stats as stats
import somaticseq.genomicFileHandler.genomic_file_handlers as genome
from somaticseq.genomicFileHandler.read_info_extractor import * 

nan = float('nan')


def from_bam(bam, my_coordinate, ref_base, first_alt, min_mq=1, min_bq=10):

    '''
    bam is the opened file handle of bam file
    my_coordiate is a list or tuple of 0-based (contig, position)
    '''
    
    indel_length = len(first_alt) - len(ref_base)
    reads = bam.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] )
    
    ref_read_mq = []
    alt_read_mq = []
    ref_read_bq = []
    alt_read_bq = []
    ref_edit_distance = []
    alt_edit_distance = []
    
    ref_concordant_reads = alt_concordant_reads = ref_discordant_reads = alt_discordant_reads = 0
    ref_for = ref_rev = alt_for = alt_rev = dp = 0
    ref_SC_reads = alt_SC_reads = ref_notSC_reads = alt_notSC_reads = 0
    MQ0 = 0
    
    ref_pos_from_end = []
    alt_pos_from_end = []
    ref_flanking_indel = []
    alt_flanking_indel = []
    
    noise_read_count = poor_read_count  = 0
    
    qname_collector = {}
    
    for read_i in reads:
        if not read_i.is_unmapped and dedup_test(read_i):
            
            dp += 1
            
            code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(read_i, my_coordinate[1]-1 )
            
            if read_i.mapping_quality < min_mq and mean(read_i.query_qualities) < min_bq:
                poor_read_count += 1
            
            if read_i.mapping_quality == 0:
                MQ0 += 1
            
            # Reference calls:
            if code_i == 1 and base_call_i == ref_base[0]:

                try:
                    qname_collector[read_i.qname].append(0)
                except KeyError:
                    qname_collector[read_i.qname] = [0]
            
                ref_read_mq.append( read_i.mapping_quality )
                ref_read_bq.append( read_i.query_qualities[ith_base] )
                
                try:
                    ref_edit_distance.append( read_i.get_tag('NM') )
                except KeyError:
                    pass
                
                # Concordance
                if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    ref_concordant_reads += 1
                elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    ref_discordant_reads += 1
                
                # Orientation
                if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    ref_for += 1
                elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    ref_rev += 1
                
                # Soft-clipped reads?
                if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                    ref_SC_reads += 1
                else:
                    ref_notSC_reads += 1

                # Distance from the end of the read:
                if ith_base != None:
                    ref_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                    
                # Flanking indels:
                ref_flanking_indel.append( flanking_indel_i )

            
            # Alternate calls:
            # SNV, or Deletion, or Insertion where I do not check for matching indel length
            elif (indel_length == 0 and code_i == 1 and base_call_i == first_alt) or \
                 (indel_length < 0  and code_i == 2 and indel_length == indel_length_i) or \
                 (indel_length > 0  and code_i == 3):

                try:
                    qname_collector[read_i.qname].append(1)
                except KeyError:
                    qname_collector[read_i.qname] = [1]

                alt_read_mq.append( read_i.mapping_quality )
                alt_read_bq.append( read_i.query_qualities[ith_base] )
                
                try:
                    alt_edit_distance.append( read_i.get_tag('NM') )
                except KeyError:
                    pass
                
                # Concordance
                if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    alt_concordant_reads += 1
                elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    alt_discordant_reads += 1
                
                # Orientation
                if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    alt_for += 1
                elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    alt_rev += 1
                
                # Soft-clipped reads?
                if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                    alt_SC_reads += 1
                else:
                    alt_notSC_reads += 1

                # Distance from the end of the read:
                if ith_base != None:
                    alt_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                        
                # Flanking indels:
                alt_flanking_indel.append( flanking_indel_i )
            
            
            # Inconsistent read or 2nd alternate calls:
            else:
                
                try:
                    qname_collector[read_i.qname].append(2)
                except KeyError:
                    qname_collector[read_i.qname] = [2]
                
                noise_read_count += 1
    
    # Done extracting info from tumor BAM. Now tally them:
    ref_mq            = mean(ref_read_mq)
    alt_mq            = mean(alt_read_mq)
    
    try:
        p_mannwhitneyu_mq = stats.mannwhitneyu(alt_read_mq, ref_read_mq, use_continuity=True, alternative='less')[1]
    
    except ValueError:
        if len(alt_read_mq) > 0 and len(ref_read_mq) > 0:
            p_mannwhitneyu_mq = 0.5
        else:
            p_mannwhitneyu_mq = nan
    
    ref_bq            = mean(ref_read_bq)
    alt_bq            = mean(alt_read_bq)
    
    try:
        p_mannwhitneyu_bq = stats.mannwhitneyu(alt_read_bq, ref_read_bq, use_continuity=True, alternative='less')[1]
    
    except ValueError:
        if len(alt_read_bq) > 0 and len(ref_read_bq) > 0:
            p_mannwhitneyu_bq = 0.5
        else:
            p_mannwhitneyu_bq = nan
    
    ref_NM            = mean(ref_edit_distance)
    alt_NM            = mean(alt_edit_distance)
    NM_Diff           = alt_NM - ref_NM - abs(indel_length)
    
    concordance_fet = stats.fisher_exact(( (ref_concordant_reads, alt_concordant_reads), (ref_discordant_reads, alt_discordant_reads) ))[1]
    strandbias_fet  = stats.fisher_exact(( (ref_for, alt_for), (ref_rev, alt_rev) ))[1]
    clipping_fet    = stats.fisher_exact(( (ref_notSC_reads, alt_notSC_reads), (ref_SC_reads, alt_SC_reads) ))[1]
    
    try:
        p_mannwhitneyu_endpos = stats.mannwhitneyu(alt_pos_from_end, ref_pos_from_end, use_continuity=True, alternative='less')[1]
    
    except ValueError:
        if len(alt_pos_from_end) > 0 and len(ref_pos_from_end) > 0:
            p_mannwhitneyu_endpos = 0.5
        else:
            p_mannwhitneyu_endpos = nan
    
    ref_indel_1bp = ref_flanking_indel.count(1)
    ref_indel_2bp = ref_flanking_indel.count(2) + ref_indel_1bp
    ref_indel_3bp = ref_flanking_indel.count(3) + ref_indel_2bp
    alt_indel_1bp = alt_flanking_indel.count(1)
    alt_indel_2bp = alt_flanking_indel.count(2) + alt_indel_1bp
    alt_indel_3bp = alt_flanking_indel.count(3) + alt_indel_2bp
    
    consistent_mates = inconsistent_mates = 0
    for pairs_i in qname_collector:
        
        # Both are alternative calls:
        if qname_collector[pairs_i] == [1,1]:
            consistent_mates += 1
        
        # One is alternate call but the other one is not:
        elif len(qname_collector[pairs_i]) == 2 and 1 in qname_collector[pairs_i]:
            inconsistent_mates += 1

    return vars()



def from_genome_reference(ref_fa, my_coordinate, ref_base, first_alt):

    '''
    ref_fa is the opened reference fasta file handle
    my_coordiate is a list or tuple of 0-based (contig, position)
    '''

    # Homopolymer eval (Make sure to modify for INDEL):
    # The min and max is to prevent the +/- 20 bases from exceeding the ends of the reference sequence
    lseq  = ref_fa.fetch(my_coordinate[0], max(0, my_coordinate[1]-20), my_coordinate[1])
    rseq  = ref_fa.fetch(my_coordinate[0], my_coordinate[1]+1, min(ref_fa.get_reference_length(my_coordinate[0])+1, my_coordinate[1]+21) )
    
    # This is to get around buy in old version of pysam that reads the reference sequence in bytes instead of strings
    lseq = lseq.decode() if isinstance(lseq, bytes) else lseq
    rseq = rseq.decode() if isinstance(rseq, bytes) else rseq
    
    seq41_ref = lseq + ref_base  + rseq
    seq41_alt = lseq + first_alt + rseq
    
    ref_counts = genome.count_repeating_bases(seq41_ref)
    alt_counts = genome.count_repeating_bases(seq41_alt)
    
    homopolymer_length = max( max(ref_counts), max(alt_counts) )
    
    # Homopolymer spanning the variant site:
    ref_c = 0
    alt_c = 0
    for i in rseq:
        if i == ref_base:
            ref_c += 1
        else:
            break
            
    for i in lseq[::-1]:
        if i == ref_base:
            ref_c += 1
        else:
            break
    
    for i in rseq:
        if i == first_alt:
            alt_c += 1
        else:
            break
            
    for i in lseq[::-1]:
        if i == first_alt:
            alt_c += 1
        else:
            break

    site_homopolymer_length = max( alt_c+1, ref_c+1 )

    return homopolymer_length, site_homopolymer_length



def somaticOddRatio(n_ref, n_alt, t_ref, t_alt, max_value=100):

    # Odds Ratio just like VarDict's output
    sor_numerator   = n_alt * t_ref
    sor_denominator = n_ref * t_alt
    if sor_numerator == 0 and sor_denominator == 0:
        sor = nan
    elif sor_denominator == 0:
        sor = max_value
    else:
        sor = sor_numerator / sor_denominator
        if sor >= max_value:
            sor = max_value

    return sor




def max_vocabularies(seq_length):
    # According to:
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Assume 4 different nucleotides
    counts = 0
    k = 1
    while k <= seq_length:
        
        if 4**k < (seq_length - k + 1):
            counts = counts + 4**k
        else:
            counts = counts + (seq_length-k+1 + 1) * (seq_length-k+1 - 1 + 1)/2
            break
        
        k += 1
                
    return counts



def LC(sequence):
    # Calculate linguistic sequence complexity according to
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Assume 4 different nucleotides
    sequence = sequence.upper()
    
    if not 'N' in sequence:
        
        number_of_subseqs     = 0
        seq_length            = len(sequence)
        max_number_of_subseqs = max_vocabularies(seq_length)
    
        for i in range(1, seq_length+1):
            
            #max_vocab_1 = 4**i
            #max_vocab_2 = seq_length - i + 1
            set_of_seq_n = set()
    
            for n, nth_base in enumerate(sequence):
                
                if n+i <= len(sequence):
                    sub_seq = sequence[n:n+i]
                    set_of_seq_n.add( sub_seq )
    
            num_uniq_subseqs  = len(set_of_seq_n)
            number_of_subseqs = number_of_subseqs + num_uniq_subseqs
    
        lc = number_of_subseqs/max_number_of_subseqs
    
    else:
        lc = float('nan')

    return lc








def max_sub_vocabularies(seq_length, max_subseq_length):
    # According to:
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # capping the length of sub_string as an input parameter
    assert max_subseq_length <= seq_length
    
    counts = 0
    k = 1
    while k <= max_subseq_length:
        
        if 4**k < (seq_length - k + 1):
            counts = counts + 4**k
        else:
            counts = counts + (2*seq_length - k - max_subseq_length + 2) * (max_subseq_length - k + 1)/2
            break
        
        k += 1
                
    return counts




def subLC(sequence, max_substring_length=20):
    # Calculate linguistic sequence complexity according to
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Cut off substring at a fixed length
    sequence = sequence.upper()
    
    if not 'N' in sequence:
        
        number_of_subseqs     = 0
        seq_length            = len(sequence)
        max_number_of_subseqs = max_sub_vocabularies(seq_length, max_substring_length)
        
        set_of_seq_n = set()
        for i in range(1, min(max_substring_length+1, seq_length+1) ):
            set_of_seq_n.update((sequence[n: n+i] for n in range(len(sequence) - i + 1)))
        
        number_of_subseqs  = len(set_of_seq_n)
        lc = number_of_subseqs/max_number_of_subseqs
    
    else:
        lc = float('nan')

    return lc
