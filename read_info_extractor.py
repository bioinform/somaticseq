#!/usr/bin/env python3

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

## Define functions:
def position_of_aligned_read(read_i, target_position):
    '''
    Return the base call of the target position, or if it's a start of insertion/deletion.
    This target position follows pysam convension, i.e., 0-based. 
    In VCF files, deletions/insertions occur AFTER the position.
    
    Return (Code, seq_i, base_at_target, indel_length, nearest insertion/deletion)
    
    The first number in result is a code:
    1) Match to reference, which is either a reference read or a SNV/SNP
    2) Deletion after the target position
    3) Insertion after the target position
    0) The target position does not match to reference, and may be discarded for "reference/alternate" read count purposes, but can be kept for "inconsistent read" metrics. 
    '''
    
    flanking_deletion, flanking_insertion = nan, nan
    
    for i, align_i in enumerate(read_i.aligned_pairs):
        
        # If find a match:
        if align_i[1] == target_position:
            seq_i = align_i[0]
            break
    
    
    # If the target position is aligned:
    try:
        if isinstance(seq_i, int): 
            base_at_target = read_i.seq[seq_i]
            
            # Whether if it's a Deletion/Insertion depends on what happens after this position:
            # If the match (i.e., i, seq_i) is the final alignment, then you cannot know if it's an indel
            # if "i" is NOT the final alignment:
            if i != len(read_i.aligned_pairs) - 1:
                
                indel_length = 0
                # If the next alignment is the next sequenced base, then the target is either a reference read of a SNP/SNV:
                if read_i.aligned_pairs[i+1][0] == seq_i+1 and read_i.aligned_pairs[i+1][1] == target_position + 1:
                    
                    code = 1 # Reference read for mismatch
                
                # If the next reference position has no read position to it, it is DELETED in this read:
                elif read_i.aligned_pairs[i+1][0] == None and read_i.aligned_pairs[i+1][1] == target_position + 1:
                    
                    code = 2 # Deletion
                    
                    for align_j in read_i.aligned_pairs[ i+1:: ]:
                        if align_j[0] == None:
                            indel_length -= 1
                        else:
                            break
                        
                # Opposite of deletion, if the read position cannot be aligned to the reference, it can be an INSERTION.
                # Insertions sometimes show up wit soft-clipping at the end, if the inserted sequence is "too long" to align on a single read. In this case, the inserted length derived here is but a lower limit of the real inserted length. 
                elif read_i.aligned_pairs[i+1][0] == seq_i+1 and read_i.aligned_pairs[i+1][1] == None:
                    
                    code = 3 # Insertion or soft-clipping
                    
                    for align_j in read_i.aligned_pairs[ i+1:: ]:
                        if align_j[1] == None:
                            indel_length += 1
                        else:
                            break
            
            # If "i" is the final alignment, cannt exam for indel:
            else:
                code = 1           # Assuming no indel
                indel_length = nan # Would be zero if certain no indel, but uncertain here
        
        # If the target position is deleted from the sequencing read (i.e., the deletion in this read occurs before the target position):
        elif seq_i == None:
            code = 0
            base_at_target, indel_length, flanking_indel = None, None, None
        
        
        # See if there is insertion/deletion within 5 bp of "i":
        if isinstance(indel_length, int): 
            flanking_indel = inf
            left_side_start = seq_i
            right_side_start = seq_i + abs(indel_length) + 1
            switch = 1
            for j in (3,2,1):
                for indel_seeker_i in left_side_start, right_side_start:
                    
                    switch = switch * -1
                    displacement = j * switch
                    seq_j = indel_seeker_i + displacement
                                    
                    if 0 <= seq_j < len(read_i.aligned_pairs):
                    
                        # If the reference position has no base aligned to it, it's a deletion.
                        # On the other hand, if the base has no reference base aligned to it, it's an insertion.
                        if read_i.aligned_pairs[ seq_j ][1] == None or read_i.aligned_pairs[ seq_j ][0] == None:
                            flanking_indel = j
                            break
        else:
            flanking_indel = None
            
        return code, seq_i, base_at_target, indel_length, flanking_indel
    
    # The target position does not exist in the read
    except UnboundLocalError:
        return None, None, None, None, None






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
    
    
    


##### Extract information from external vcf files:
##### From Samtools vcf:
def sam_info_DP4(vcf_object):
    dp4_string = vcf_object.get_info_value('DP4')
    if dp4_string:
        dp4_string = dp4_string.split(',')
        dp4 = ( int(dp4_string[0]), int(dp4_string[1]), int(dp4_string[2]), int(dp4_string[3]) )
    else:
        dp4 = nan,nan,nan,nan
        
    return dp4
    


def sam_info_DP(vcf_object):
    result = vcf_object.get_info_value('DP')
    if result:
        return eval(result)
    else:
        return nan
    


def sam_info_MQ(vcf_object):
    result = vcf_object.get_info_value('MQ')
    if result:
        return eval(result)
    else:
        return nan



def sam_info_PV4(vcf_object):
    '''P-values for strand bias, baseQ bias, mapQ bias and tail distance bias'''
    pv4_string = vcf_object.get_info_value('PV4')
    if pv4_string:
        pv4_string = pv4_string.split(',')
        pv4 = ( float(pv4_string[0]), float(pv4_string[1]), float(pv4_string[2]), float(pv4_string[3]) )
    else:
        pv4 = nan,nan,nan,nan
    
    return pv4
    


##### From Haplotype caller vcf:
def haplo_MQ0(vcf_object):
    '''Total Mapping Quality Zero Reads'''
    
    mq0 = vcf_object.get_info_value('MQ0')
    if mq0:
        mq0 = eval(mq0)
    else:
        mq0 = nan
        
    return mq0



def haplo_MQ(vcf_object):
    '''RMS Mapping Quality'''
    result = vcf_object.get_info_value('MQ')
    if result:
        return eval(result)
    else:
        return nan
    
    
    
def haplo_MLEAF(vcf_object):
    '''Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed'''
    
    mleaf = vcf_object.get_info_value('MLEAF')
    
    if mleaf:
        mleaf = mleaf.split(',')
        mleaf = [eval(i) for i in mleaf]
        mleaf = max(mleaf)
        
    else:
        mleaf = nan
        
    return mleaf
    
    
    
def haplo_MLEAC(vcf_object):
    '''Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed'''
    
    mleac = vcf_object.get_info_value('MLEAC')
    
    if mleac:        
        mleac = mleac.split(',')
        mleac = [eval(i) for i in mleac]
        mleac = max(mleac)
        
    else:
        mleac = nan
        
    return mleac



def haplo_DP(vcf_object):
    result = vcf_object.get_sample_value('DP')
    if result:
        return eval(result)
    else:
        return nan


def haplo_BaseQRankSum(vcf_object):
    '''Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities'''
    result = vcf_object.get_info_value('BaseQRankSum')
    return eval(result) if result else nan
    
    
def haplo_ClippingRankSum(vcf_object):
    '''Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases'''
    result = vcf_object.get_info_value('ClippingRankSum')
    return eval(result) if result else nan
    
    
def haplo_LikelihoodRankSum(vcf_object):
    '''Z-score from Wilcoxon rank sum test of Alt Vs. Ref haplotype likelihoods'''
    result = vcf_object.get_info_value('LikelihoodRankSum')
    return eval(result) if result else nan
    
    
def haplo_ReadPosRankSum(vcf_object):
    '''Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias'''
    result = vcf_object.get_info_value('ReadPosRankSum')
    return eval(result) if result else nan
    
    
def haplo_MQRankSum(vcf_object):
    '''Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities'''
    result = vcf_object.get_info_value('MQRankSum')
    return eval(result) if result else nan


##### Stuff from my own vcf:
def calculate_baf(caf_string):
    try:
        caf = re.search(r'\[[0-9.,]+\]', caf_string)
        if caf:
            caf_match = re.sub(r'\.([^0-9])', r'0\g<1>', caf.group())
            caf = list( eval(caf_match) )
            caf.sort()
            baf = sum(caf[0:-1])  # Minor Allele Frequency
            return baf
            
    except TypeError:
        return nan

  

def find_AMQ(vcf_object, i):
    amq = vcf_object.get_sample_value('AMQ', idx=i)
    
    if amq:
        amq = amq.split(',')
        amq_ref = eval(amq[0])
        try:
            amq_alt = eval(amq[1])
        except IndexError:
            amq_alt = nan
    
    else:
        amq_ref, amq_alt = nan, nan
        
    return amq_ref, amq_alt



def find_BQ(vcf_object, i):
    bq = vcf_object.get_sample_value('BQ', idx=i)
    # If there are two numbers, it came from SomaticSniper. If there is one number, it came from MuTect. 
    
    if bq:
        
        if bq == '.':
            bq_ref, bq_alt = nan, nan
            
        elif ',' in bq:
            bq = bq.split(',')
            bq_ref = eval(bq[0])
            bq_alt = eval(bq[1])
            
        else:
            bq_ref, bq_alt = eval(bq), eval(bq)
            
    else:
        bq_ref, bq_alt = nan, nan
        
    return bq_ref, bq_alt
    
    

def find_SOR(vcf_object):
    # VarDict's odd ratio, could be Inf, but other than Inf max was 180, so I will convert Inf --> 200. Stored in the TUMOR sample. 
    
    sor = vcf_object.get_info_value('SOR')
    if sor:
        sor = float(sor) if sor != 'Inf' else 200
    else:
        sor = nan
        
    return sor
    


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
