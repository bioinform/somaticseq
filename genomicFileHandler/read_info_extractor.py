#!/usr/bin/env python3

import re

cigar_aln_match    = 0
cigar_insertion    = 1
cigar_deletion     = 2
cigar_skip         = 3
cigar_soft_clip    = 4
cigar_hard_clip    = 5
cigar_padding      = 6
cigar_seq_match    = 7
cigar_seq_mismatch = 8

nan = float('nan')
inf = float('inf')

## Define functions:


### PYSAM ###
def position_of_aligned_read(read_i, target_position, win_size=3):
    '''
    Return the base call of the target position, and if it's a start of insertion/deletion.
    This target position follows pysam convension, i.e., 0-based.
    In VCF files, deletions/insertions occur AFTER the position.

    Return (Code, seq_i, base_at_target, indel_length, nearest insertion/deletion)

    The first number in result is a codeMatch to reference on CIGAR, which is either a reference read or a SNV (substitution counts as M in CIGAR) to reference, which is either a reference read or a SNV/SNP
        2: Deletion after the target position
        3: Insertion after the target position
        0: The target position does not match to reference, and may be discarded for "reference/alternate" read count purposes, but can be kept for "inconsistent read" metrics.
    '''

    flanking_deletion, flanking_insertion = nan, nan

    aligned_pairs = read_i.get_aligned_pairs()
    
    for i, align_i in enumerate(aligned_pairs):

        # If find a match:
        if align_i[1] == target_position:
            seq_i = align_i[0]
            idx_aligned_pair = i
            break

    # If the target position is aligned:
    try:
        if seq_i is not None:
            base_at_target = read_i.seq[seq_i]

            # Whether if it's a Deletion/Insertion depends on what happens after this position:
            # If the match (i.e., i, seq_i) is the final alignment, then you cannot know if it's an indel
            # if "i" is NOT the final alignment:
            if i != len(aligned_pairs) - 1:

                indel_length = 0
                # If the next alignment is the next sequenced base, then the target is either a reference read of a SNP/SNV:
                if aligned_pairs[i+1][0] == seq_i+1 and aligned_pairs[i+1][1] == target_position + 1:

                    code = 1 # Reference read for mismatch

                # If the next reference position has no read position to it, it is DELETED in this read:
                elif aligned_pairs[i+1][0] == None and aligned_pairs[i+1][1] == target_position + 1:

                    code = 2 # Deletion

                    for align_j in aligned_pairs[ i+1:: ]:
                        if align_j[0] == None:
                            indel_length -= 1
                        else:
                            break

                # Opposite of deletion, if the read position cannot be aligned to the reference, it can be an INSERTION.
                # Insertions sometimes show up wit soft-clipping at the end, if the inserted sequence is "too long" to align on a single read. In this case, the inserted length derived here is but a lower limit of the real inserted length.
                elif aligned_pairs[i+1][0] == seq_i+1 and aligned_pairs[i+1][1] == None:

                    code = 3 # Insertion or soft-clipping

                    for align_j in aligned_pairs[ i+1:: ]:
                        if align_j[1] == None:
                            indel_length += 1
                        else:
                            break

            # If "i" is the final alignment, cannt exam for indel:
            else:
                code = 1           # Assuming no indel
                indel_length = nan # Would be zero if certain no indel, but uncertain here

        # If the target position is deleted from the sequencing read (i.e., the deletion in this read occurs before the target position):
        else:
            code = 0
            base_at_target, indel_length, flanking_indel = None, None, None

        # See if there is insertion/deletion within 5 bp of "seq_i" on the query.
        # seq_i is the i_th aligned base
        if isinstance(indel_length, int):
            right_indel_flanks = inf
            left_indel_flanks  = inf
            left_side_start    = idx_aligned_pair - 1
            right_side_start   = idx_aligned_pair + abs(indel_length) + 1
                        
            #(i, None) = Insertion (or Soft-clips), i.e., means the i_th base in the query is not aligned to a reference
            #(None, coordinate) = Deletion, i.e., there is no base in it that aligns to this coordinate.
            # If those two scenarios occur right after an aligned base, that base position is counted as an indel.
            for step_right_i in range( min(win_size, len(aligned_pairs)-right_side_start-1 ) ):
                j = right_side_start + step_right_i
                
                if (aligned_pairs[j+1][1] == None or aligned_pairs[j+1][0] == None):
                    right_indel_flanks = step_right_i + 1
                    break
            
            for step_left_i in range( min(win_size, left_side_start) ):
                j = left_side_start - step_left_i
                
                if (aligned_pairs[j][1] == None or aligned_pairs[j][0] == None):
                    left_indel_flanks = step_left_i + 1
                    break
            
            flanking_indel = min(left_indel_flanks, right_indel_flanks)

        else:
            flanking_indel = None

        return code, seq_i, base_at_target, indel_length, flanking_indel

    # The target position does not exist in the read
    except UnboundLocalError:
        return None, None, None, None, None


## Dedup test for BAM file
def dedup_test(read_i, remove_dup_or_not=True):
    '''
    Return False (i.e., remove the read) if the read is a duplicate and if the user specify that duplicates should be removed.
    Else return True (i.e, keep the read)
    '''
    if read_i.is_duplicate and remove_dup_or_not:
        return False
    else:
        return True



### END OF PYSAM ###


# Useful to make BED region into an iterator of coordinates
def genomic_coordinates(contig_i, start, end):
    for pos_i in range(start, end+1):
        yield contig_i, pos_i




def mean(stuff):    
    return sum(stuff)/len(stuff) if stuff else nan



##### Extract Indel DP4 info from pileup files:
def pileup_indel_DP4(pileup_object, indel_pattern):
    if pileup_object.reads:
        ref_for = pileup_object.reads.count('.')
        ref_rev = pileup_object.reads.count(',')
        alt_for = pileup_object.reads.count( indel_pattern.upper() )
        alt_rev = pileup_object.reads.count( indel_pattern.lower() )

        dp4     = ref_for, ref_rev, alt_for, alt_rev

    else:
        dp4 = nan,nan,nan,nan

    return dp4


def pileup_DP4(pileup_object, ref_base, variant_call):

    base_calls = pileup_object.base_reads()

    if base_calls:

        # SNV
        if len(variant_call) == len(ref_base):

            ref_for,ref_rev,alt_for,alt_rev = base_calls[0], base_calls[1], base_calls[2].count(variant_call.upper()), base_calls[3].count(variant_call.lower())

        # Insertion:
        elif len(variant_call) > len(ref_base):

            inserted_sequence = variant_call[ len(ref_base):: ]

            ref_for,ref_rev,alt_for,alt_rev = base_calls[0], base_calls[1], base_calls[6].count(inserted_sequence.upper()), base_calls[7].count(inserted_sequence.lower())

        # Deletion:
        elif len(variant_call) < len(ref_base):

            deleted_sequence = ref_base[ len(variant_call):: ]

            ref_for,ref_rev,alt_for,alt_rev = base_calls[0], base_calls[1], base_calls[4].count(deleted_sequence.upper()), base_calls[5].count(deleted_sequence.lower())

    else:
        ref_for = ref_rev = alt_for = alt_rev = 0

    return ref_for, ref_rev, alt_for, alt_rev




def rescale(x, original='fraction', rescale_to=None, max_phred=1001):
    
    if ( rescale_to == None ) or ( original.lower() == rescale_to.lower() ):
        y = x if isinstance(x, int) else '%.2f' % x
    
    elif original.lower() == 'fraction' and rescale_to == 'phred':
        y = genome.p2phred(x, max_phred=max_phred)
        y = '%.2f' % y
    
    elif original.lower() == 'phred' and rescale_to == 'fraction':
        y = genome.phred2p(x)
        y = '%.2f' % y
    
    return y





##### Stuff from VarDict:
def find_MSI(vcf_object):

    msi = vcf_object.get_info_value('MSI')
    if msi:
        msi = float(msi)
    else:
        msi = nan
    return msi


def find_MSILEN(vcf_object):

    msilen = vcf_object.get_info_value('MSILEN')
    if msilen:
        msilen = float(msilen)
    else:
        msilen = nan
    return msilen


def find_SHIFT3(vcf_object):

    shift3 = vcf_object.get_info_value('SHIFT3')
    if shift3:
        shift3 = float(shift3)
    else:
        shift3 = nan
    return shift3



# MuTect2's Stuff:
def mutect2_nlod(vcf_object):
    nlod = vcf_object.get_info_value('NLOD')
    if nlod:
        return float(nlod)
    else:
        return nan


def mutect2_tlod(vcf_object):
    tlod = vcf_object.get_info_value('TLOD')
    if tlod:
        return float(tlod)
    else:
        return nan


def mutect2_STR(vcf_object):
    if vcf_object.get_info_value('STR'):
        return 1
    else:
        return 0


def mutect2_ECNT(vcf_object):
    ecnt = vcf_object.get_info_value('ECNT')
    if ecnt:
        try:
            ecnt = int( ecnt )
        except ValueError:
            ecnt = nan
    else:
        ecnt = nan

    return ecnt
