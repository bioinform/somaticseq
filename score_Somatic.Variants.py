#!/usr/bin/env python3

# Add support for single-sample mode
# Extract more COSMIC features
# AD field now assumes Broad's convention 
# 4/12/2015
# Li Tai Fang


import sys, os, argparse, gzip, math, itertools
import regex as re

sys.path.append('/net/kodiak/volumes/lake/shared/opt/Bina_SomaticMerge')
import genomic_file_handlers as genome

# argparse Stuff
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-indications',   '--cancer-indications',                type=str,   help='Output VCF file', nargs='*', required=False, default=None)
parser.add_argument('-mincaller',     '--minimum-callers',                   type=int,   help='Minimum number of callers identify it as a SOMATIC for the site to be emitted', required=False, default=0)
parser.add_argument('-infile',        '--input-vcf',                         type=str,   help='Input VCF file', required=True, default=None)
parser.add_argument('-outfile',       '--output-vcf',                        type=str,   help='Output VCF file', required=True, default=None)

parser.add_argument('-genes',         '--cancer_gene_census',                type=str,   help='A file where the first column contain list of cancer census genes', required=False)

parser.add_argument('-cntlo' ,        '--low-threshold-CNT',                 type=int,   help='Lowest threshold for COSMIC +2 scoring', default=3)
parser.add_argument('-cntmed' ,       '--medium-threshold-CNT',              type=int,   help='Medium threshold for COSMIC +3 scoring', default=5)
parser.add_argument('-cnthi' ,        '--high-threshold-CNT',                type=int,   help='Highest threshold for COSMIC +5 scoring', default=10)
parser.add_argument('-cntex' ,        '--extreme-threshold-CNT',             type=int,   help='Highest threshold for COSMIC +6 scoring', default=100)
parser.add_argument('-cntri' ,        '--ridiculous-threshold-CNT',          type=int,   help='Highest threshold for COSMIC +7 scoring', default=1000)
parser.add_argument('-caf' ,          '--high-threshold-CAF',                type=float, help='Lowest threshold for common SNP -3 penalty', default=0.01)
parser.add_argument('-MSC',           '--minor-strand-count',                type=int,   help='Minor Strand Count, i.e,. filtering strand bias', default=1)
parser.add_argument('-DP' ,           '--depth-of-coverage',                 type=int,   help='Total Depth of Coverage', default=15)
parser.add_argument('-ALF',           '--alternative-allele-freq',           type=float, help='Alternative Allele Frequency', default=0.10)
parser.add_argument('-NAF',           '--normal-alternative-freq',           type=float, help='Normal Alternative Frequency (Upper Bound)', default=0.05)


parser.add_argument('-normal',        '--normal-sample-name',                type=str,   help='Normal Sample Name', required=False, default='NORMAL')
parser.add_argument('-tumor',         '--tumor-sample-name',                 type=str,   help='Tumor Sample Name', required=False, default='TUMOR')

parser.add_argument('-pass_score',    '--oncoscore-threshold-pass',          type=int,   help='Minimum OncoScore that generates a PASS filter', default=5)
parser.add_argument('-lowconf_score', '--oncoscore-threshold-lowconfidence', type=int,   help='Minimum OncoScore that generates a PASS filter', default=3)

parser.add_argument('-tools',         '--individual-mutation-tools',         type=str,   help='A list of all tools: have to match the annotated tool name in the input vcf files', nargs='*', required=False, default=('CGA', 'VarScan2', 'JointSNVMix2', 'SomaticSniper', 'VarDict'))

args = parser.parse_args()



# Rename input arguments:
tools = args.individual_mutation_tools
num_tools = len(tools)
min_tools = args.minimum_callers
assert args.oncoscore_threshold_pass > args.oncoscore_threshold_lowconfidence


# The regular expression pattern for "chrXX 1234567" in both VarScan2 Output and VCF files:
pattern_chr_position = re.compile(r'^(?:chr)?(?:[0-9]+|[XYM]|MT)\t[0-9]+\b')
pattern_COSM  = re.compile(r'COS[MN][0-9]+')
pattern_rsID = re.compile(r'rs[0-9]+')
pattern_CAF = re.compile(r'\[[0-9.,]+\]')

# snpEFF's fields
# EFF= Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_Length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] ),
# An example:
#EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|cGt/cAt|R229H|246|SLC35E2||CODING|NM_001199787.1|6|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|cGt/cAt|R229H|266|SLC35E2||CODING|NM_182838.2|5|1)'
pattern_EFF = re.compile(r'(EFF=)?(?<Effect>[^(]+)\((?<Effect_Impact>[^|]*)\|(?<Functional_Class>[^|]*)\|(?<Codon_Change>[^|]*)\|(?<Amino_Acid_Change>[^|]*)\|(?<Amino_Acid_Length>[^|]*)\|(?<Gene_Name>[^|]*)\|(?<Transcript_BioType>[^|]*)\|(?<Gene_Coding>[^|]*)\|(?<Transcript_ID>[^|]*)\|(?<Exon_Rank>[^|]*)\|(?<Genotype_Number>[^|]*)?.*\)' )


if not args.cancer_indications:
    print('ALERT: Cancer indications were empty. None will be scored for it!!!', file=sys.stderr)



# Obtain the cancer gene census set:
cancer_gene_set = {}
if args.cancer_gene_census:
    with open(args.cancer_gene_census) as cancer_gene_file:
        line_i = cancer_gene_file.readline().rstrip()  # Skip first line
        line_i = cancer_gene_file.readline().rstrip()
        
        while line_i:
            item_i = line_i.split('\t')
            
            # 0th column is gene symbol. 5th column is indications:
            cancer_gene_set[ item_i[0] ] = item_i[5]
            
            line_i = cancer_gene_file.readline().rstrip()


# Maximum scores for subEvidence and subKnowledge
max_subEvidence = 9
max_subKnowledge = 9
max_scale = 10



# All possible combination of tools:
MVJS_combinations = {}

idx_combo_total, idx_combo_dbsnp, idx_combo_common, idx_combo_cosmic = 0, 1, 2, 3

for combo_i in itertools.product( (True, False), repeat = len(tools) ):
    
    # The four zeros represent [Total, dbsnp, COMMON, COSMIC]
    MVJS_combinations[combo_i] = [0, 0, 0, 0]


# Keeping a tab on all those scores
bina_score_tally = {}
subscore_evidence_tally = {}
bonus_knowledge_tally = {}
penalty_knowledge_tally = {}
num_methods_tally = {}


with genome.open_textfile(args.input_vcf) as vcf, open(args.output_vcf, 'w') as vcf_out:
    
    line_i = vcf.readline().rstrip()
    
    while line_i.startswith('#'):
        # Read thru the headers and metadata:
        
        if line_i.startswith('#CHROM'):
            header_item = line_i.split('\t')
            
            if len(header_item) == 11:
                paired_mode = True
                idxN, idxT = 0,1
                
            elif len(header_item) == 10:
                paired_mode = False
                idxT = 0
                
            else:
                raise Exception('The number of samples is neither 1 nor 2')
        
        line_i = vcf.readline().rstrip()
    
    
    # Write the headers and meta-data into the output vcf file:
    vcf_out.write('##fileformat=VCFv4.1\n')
    
    vcf_out.write('##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation with ONCOSCORE at least {}">\n'.format(args.oncoscore_threshold_pass) )
    vcf_out.write('##FILTER=<ID=LowQual,Description="Lower confident somatic mutation with ONCOSCORE at least between {} and {}">\n'.format(args.oncoscore_threshold_lowconfidence ,args.oncoscore_threshold_pass-1) )
    vcf_out.write('##FILTER=<ID=REJECT,Description="Rejected as a confident somatic mutation with ONCOSCORE below {}">\n'.format(args.oncoscore_threshold_lowconfidence) )
    
    vcf_out.write('##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Passes as a SOMATIC Variant">\n')
    vcf_out.write('##INFO=<ID=COMMON,Number=1,Type=Integer,Description="RS is a common SNP.  A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.">\n')
    vcf_out.write('##INFO=<ID=CAF,Number=.,Type=String,Description="An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter">\n')
    vcf_out.write('##INFO=<ID=G5,Number=0,Type=Flag,Description=">5% minor allele frequency in 1+ populations">\n')
    vcf_out.write('##INFO=<ID=G5A,Number=0,Type=Flag,Description=">5% minor allele frequency in each and all populations">\n')
    
    vcf_out.write('##INFO=<ID=AA,Number=1,Type=String,Description="Peptide annotation">\n')
    vcf_out.write('##INFO=<ID=CDS,Number=1,Type=String,Description="CDS annotation">\n')
    vcf_out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">\n')
    vcf_out.write('##INFO=<ID=MUT,Number=0,Type=Flag,Description="Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources">\n')
    vcf_out.write('##INFO=<ID=CNT,Number=1,Type=Integer,Description="How many samples have this mutation">\n')

    vcf_out.write('##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: \'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )\'">\n')
    vcf_out.write('##INFO=<ID=subEvidence,Number=1,Type=Float,Description="subScore based on number of somatic mutation callers">\n')
    vcf_out.write('##INFO=<ID=subKnowledge,Number=2,Type=Float,Description="First number is Bonus subScore based on COSMIC and Cancer Gene Census membership plus SIFT prediction. Second number is penalty based on dbSNP membership and COMMON annotation.">\n')
    vcf_out.write('##INFO=<ID=SOURCES,Number=.,Type=String,Description="The Somatic Mutation Callers that called it">\n')
    vcf_out.write('##INFO=<ID=NUM_SMMETHODS,Number=1,Type=Integer,Description="Identified by number of Somatic Mutation Callers">\n')
    vcf_out.write('##INFO=<ID=OncoRK,Number=1,Type=Float,Description="BINA Somatic Mutation Score">\n')
    vcf_out.write('##INFO=<ID=LSEQ,Number=1,Type=String,Description="5\' flanking seq">\n')
    vcf_out.write('##INFO=<ID=RSEQ,Number=1,Type=String,Description="3\' flanking seq">\n')
    
    vcf_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_out.write('##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">\n')
    vcf_out.write('##FORMAT=<ID=SC,Number=4,Type=Integer,Description="Strandness: counts of forward-/reverse-aligned reference and indel-supporting reads (FwdRef,RevRef,FwdIndel,RevIndel)">\n')
    vcf_out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n')
    vcf_out.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
    vcf_out.write('##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">\n')
    vcf_out.write('##FORMAT=<ID=AMQ,Number=.,Type=Integer,Description="Average mapping quality for each allele present in the genotype">\n')
    vcf_out.write('##FORMAT=<ID=BQ,Number=.,Type=Float,Description="Average base quality for reads supporting alleles">\n')
    vcf_out.write('##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Average mapping quality across all reads">\n')
    vcf_out.write('##FORMAT=<ID=VAQ,Number=1,Type=Integer,Description="Variant allele quality">\n')
    vcf_out.write('##FORMAT=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">\n')
    vcf_out.write('##FORMAT=<ID=MSI,Number=1,Type=Float,Description="MicroSattelite. > 1 indicates MSI">\n')
    vcf_out.write('##FORMAT=<ID=MSILEN,Number=1,Type=Float,Description="MSI unit repeat length in bp">\n')
    vcf_out.write('##FORMAT=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">\n')
    vcf_out.write('##FORMAT=<ID=PMEAN,Number=1,Type=Float,Description="Mean position in reads">\n')
    vcf_out.write('##FORMAT=<ID=PSTD,Number=1,Type=String,Description="Position STD in reads. Use string because can be NA.">\n')
    vcf_out.write('##FORMAT=<ID=QSTD,Number=1,Type=String,Description="Quality score STD in reads. Use string because it can be NA.">\n')
    vcf_out.write('##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">\n')

    vcf_out.write('##FORMAT=<ID=SOR,Number=1,Type=String,Description="Odds ratio, could be Inf">\n')
    vcf_out.write('##FORMAT=<ID=SSSC,Number=1,Type=Integer,Description="SomaticSniper reported Somatic Score">\n')
    vcf_out.write('##FORMAT=<ID=VSSC,Number=1,Type=Integer,Description="VarScan reported Somatic Score">\n')
    vcf_out.write('##FORMAT=<ID=JSSC,Number=1,Type=Float,Description="JointSNVMix reported somatic probability">\n')
    vcf_out.write('##FORMAT=<ID=DSSF,Number=1,Type=Float,Description="VarDict reported variant p-value">\n')
    
    
    # Main Header:
    if paired_mode:
        vcf_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\t{}\n'.format(args.normal_sample_name, args.tumor_sample_name) )
    elif not paired_mode:
        vcf_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(args.tumor_sample_name) )
    

    # Start calculating scores:
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )
        
        # Initialize stuff:
        bina_score = 0
        subscore_evidence = 0
        bonus_knowledge = 0
        penalty_knowledge = 0
        snpEff_scores = [0,]
        
        # Apply scoring "matrix" here:
        info_item = vcf_i.get_info_items()
                
        ### Process the callers (score between 0 to +4):
        combo_i = [ vcf_i.get_info_value( tool_i ) for tool_i in tools ]
        combo_i = tuple(combo_i)
        
        MVJS_combinations[combo_i][idx_combo_total] += 1
        
        num_smmethods = sum(combo_i)
        
        
        # Evidence-based score: number of tools:
        # To make indel callers and snp caller the same scale:scored.snp2.vcf
        subscore_callers = num_smmethods*(5/num_tools)
        
        
        # Proceed only if there is AT LEAST "min_tools" CALLER calling it SOMATIC:
        if num_smmethods >= min_tools:
            
            subscore_filter = 0
            
            # DP4 (SC for Indellocator) has strand info and counts info, so highest precedence. Grab those if available. 
            got_dp4, got_sc = 'DP4' in vcf_i.get_sample_variable(), 'SC' in vcf_i.get_sample_variable()
            if (got_dp4 or got_sc):
                
                probe_dp4 = 'DP4' if got_dp4 else 'SC'
                
                tumor_ref_for,tumor_ref_rev,tumor_alt_for,tumor_alt_rev = vcf_i.get_sample_value(probe_dp4, idxT).split(',')
                tumor_ref_for = int( tumor_ref_for )
                tumor_ref_rev = int( tumor_ref_rev )
                tumor_alt_for = int( tumor_alt_for )
                tumor_alt_rev = int( tumor_alt_rev )
                
                tumor_alt_total = tumor_alt_for + tumor_alt_rev
                tumor_coverage_depth = sum( (tumor_ref_for, tumor_ref_rev, tumor_alt_for, tumor_alt_rev) )
                
                # For subScoring:
                try:
                    tumor_alt_allele_freq = tumor_alt_total/tumor_coverage_depth
                except ZeroDivisionError:
                    tumor_alt_allele_freq = 0
                    
                minor_strand_count = min(tumor_alt_for, tumor_alt_rev)


                if paired_mode:
                    normal_ref_for,normal_ref_rev,normal_alt_for,normal_alt_rev = vcf_i.get_sample_value(probe_dp4, idxN).split(',')
                    normal_ref_for = int( normal_ref_for )
                    normal_ref_rev = int( normal_ref_rev )
                    normal_alt_for = int( normal_alt_for )
                    normal_alt_rev = int( normal_alt_rev )
                    
                    normal_alt_total = normal_alt_for + normal_alt_rev
                    normal_coverage_depth = sum( (normal_ref_for, normal_ref_rev, normal_alt_for, normal_alt_rev) )

                    try:
                        normal_alt_allele_freq = normal_alt_total/normal_coverage_depth
                    except ZeroDivisionError:
                        normal_alt_allele_freq = 1
                    
            
            # MuTect and JointSNVMix2 has no strand info, so grab those if DP4/SC doesn't exist.
            # Following Broad's convention, the first number in AD is refernece, and the followings are alternative. We'll take the 2nd number as the alternative read counts. 
            elif 'AD' in vcf_i.get_sample_variable():
                
                
                try:
                    tumor_alt_total  = int( vcf_i.get_sample_value('AD', idxT).split(',')[1] )
                    tumor_ref_total  = int( vcf_i.get_sample_value('AD', idxT).split(',')[0] )
                    
                    # For subScoring:
                    tumor_coverage_depth = tumor_alt_total + tumor_ref_total            
                    
                    try:
                        tumor_alt_allele_freq = tumor_alt_total/tumor_coverage_depth
                    except ZeroDivisionError:
                        tumor_alt_allele_freq = 0
                
                except ValueError:
                    tumor_alt_allele_freq = 0


                if paired_mode:
                    
                    try:
                        normal_alt_total = int( vcf_i.get_sample_value('AD', idxN).split(',')[1] )
                        normal_ref_total = int( vcf_i.get_sample_value('AD', idxN).split(',')[0] )
        
                        try:
                            normal_alt_allele_freq = normal_alt_total/(normal_alt_total + normal_ref_total)
                        except ZeroDivisionError:
                            normal_alt_allele_freq = 0
                    
                    except ValueError:
                        tumor_alt_allele_freq = 0
            
            
            # Well then, grab DP and see if either RD or AD exists (not both, otherwise it wouldn't get there).
            elif 'DP' in vcf_i.get_sample_variable():
                
                tumor_coverage_depth  = int( vcf_i.get_sample_value('DP', idxT) )
                
                if paired_mode:
                    normal_coverage_depth = int( vcf_i.get_sample_value('DP', idxN) )
                
                if 'RD' in vcf_i.get_sample_variable():
                                        
                    tumor_ref_total = int( vcf_i.get_sample_value('RD', idxT) )
                    tumor_alt_total = tumor_coverage_depth - tumor_ref_total
                    
                    try:
                        tumor_alt_allele_freq = tumor_alt_total/tumor_coverage_depth
                    except ZeroDivisionError:
                        tumor_alt_allele_freq = 0
                    
                    
                    if paired_mode:
                        normal_ref_total = int( vcf_i.get_sample_value('RD', idxN) )
                        normal_alt_total = normal_coverage_depth - normal_ref_total
        
                        try:
                            normal_alt_allele_freq = normal_alt_total/normal_coverage_depth
                        except ZeroDivisionError:
                            normal_alt_allele_freq = 0
                            
                else:
                    tumor_alt_allele_freq = 0
                    normal_alt_allele_freq = 0
                
            
            ### myFilter Sub-Scoring
            ### Because some metrics in Filter are not always run, so make sure to delete them after each iteration. 
            try:
                if tumor_alt_allele_freq >= 2*args.alternative_allele_freq:
                    subscore_filter = subscore_filter + 1
                elif tumor_alt_allele_freq >= args.alternative_allele_freq:
                    subscore_filter = subscore_filter + .5
                    
                del tumor_alt_allele_freq
            except NameError:
                pass
            
            #
            try:
                if minor_strand_count >= 2*args.minor_strand_count:
                    subscore_filter = subscore_filter + 1
                elif minor_strand_count >= args.minor_strand_count:
                    subscore_filter = subscore_filter + .5
                elif minor_strand_count == 0:
                    subscore_filter = subscore_filter - 1
                    
                del minor_strand_count
            except NameError:
                pass
            
            #
            try:
                if tumor_coverage_depth >= 2*args.depth_of_coverage:
                    subscore_filter = subscore_filter + 1
                elif tumor_coverage_depth >= args.depth_of_coverage:
                    subscore_filter = subscore_filter + .5
                    
                del tumor_coverage_depth
            except NameError:
                pass
            
            #
            if paired_mode:
                try:
                    if normal_alt_allele_freq <= 0.5*args.normal_alternative_freq:
                        subscore_filter = subscore_filter + 1
                    elif normal_alt_allele_freq <= args.normal_alternative_freq:
                        subscore_filter = subscore_filter + .5
                        
                    del normal_alt_allele_freq
                except NameError:
                    pass    


            ### Evidence based score, scale it to 10, and then round to the nearest 0.5. 
            subscore_evidence = subscore_callers + subscore_filter
            subscore_evidence = subscore_evidence * (max_scale / max_subEvidence)   # Scale it to 10
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            
            
            ### Based on knowledge, i.e., COSMIC, and Cancer Gene Census:        
            ## Cosmic Test (score between 0 to +6):
            bonus_knowledge = 0
            cosm = re.search( pattern_COSM, vcf_i.identifier )
            if cosm:
                bonus_knowledge = bonus_knowledge + 1
                
                MVJS_combinations[combo_i][idx_combo_cosmic] += 1
                
                # Additional score based on number of samples
                cnt = vcf_i.get_info_value('CNT')
                if cnt:
                    
                    try:
                        cnt = int( cnt )
                    except ValueError:
                        cnt = 0
                    
                    if cnt >= args.extreme_threshold_CNT:
                        bonus_knowledge = bonus_knowledge + 5
                    elif cnt >= args.high_threshold_CNT:
                        bonus_knowledge = bonus_knowledge + 4
                    elif cnt >= args.medium_threshold_CNT:
                        bonus_knowledge = bonus_knowledge + 2
                    elif cnt >= args.low_threshold_CNT:
                        bonus_knowledge = bonus_knowledge + 1
            
    
            ### Cancer Gene Census Test (score either 0 or 1):
            # Look for Gene Name after snpEff in the EFF field. 
            snpEff = vcf_i.get_info_value('EFF')
            if snpEff:
                
                # Once you find where the EFF field is
                # There might be multiple annotations, one for each splice variant or whatever:
                snpEFFs = snpEff.split(',')
                snpEff_scores = []
                
                for IFF_i in snpEFFs:
                    
                    snp_score_i = 0
                    
                    snpeff_info = re.match(pattern_EFF, IFF_i)
                                    
                    # Following are annotated by snpEff:
                    gene = snpeff_info.groupdict()['Gene_Name']
                    effect_impact = snpeff_info.groupdict()['Effect_Impact']
                    
                    if gene in cancer_gene_set.keys():
                        snp_score_i = snp_score_i + 1
                        
                        # Only care about the SIFT prediction effect if the gene is in gene cancer census
                        if effect_impact.upper() == 'MODERATE' or effect_impact.upper() == 'HIGH':
                            snp_score_i = snp_score_i + 1
                        
                        if args.cancer_indications:
                            for indi_i in args.cancer_indications:
                                if re.search(indi_i, cancer_gene_set[gene], re.I):
                                    snp_score_i = snp_score_i + 1
                                    break
                                    
                    snpEff_scores.append( snp_score_i )
                        
                    # Break from the loop after done with EFF field:
                    break
                    
                    
            ### Knowledge based bonus, scale it to 10, and then round to the nearest 0.5. 
            bonus_knowledge = bonus_knowledge + max(snpEff_scores)
            bonus_knowledge = bonus_knowledge * (max_scale / max_subKnowledge)  # Scale it to 10
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            
            
            
            ## dbSNP Test (score between -3 to 0):
            penalty_knowledge = 0
            dbsnp = re.search(pattern_rsID, vcf_i.identifier)
            if dbsnp:
                
                MVJS_combinations[combo_i][idx_combo_dbsnp] += 1
                
                if vcf_i.get_info_value('COMMON') == '1' or vcf_i.get_info_value('G5') or vcf_i.get_info_value('G5A'):
                    penalty_knowledge = penalty_knowledge - 1
                    
                    MVJS_combinations[combo_i][idx_combo_common] += 1
                    
                    caf = vcf_i.get_info_value('CAF')
                    if caf:
                        caf_match = re.sub(r'\.([^0-9])', r'0\g<1>', caf )
                        
                        try:
                            caf = eval(caf_match)
                            caf.sort()
                            maf = sum(caf[0:-1])  # Minor Allele Frequency
                            
                            if maf >= args.high_threshold_CAF:
                                penalty_knowledge = penalty_knowledge - 2
                        
                        except SyntaxError:
                            pass
            
            
            # Calculate a "BINA_ONCOSCORE" only if at least one of the callers called it a somatic event:
            bina_score        = 0.6 * subscore_evidence + 0.4 * bonus_knowledge + penalty_knowledge
            
            bina_score        = 0.5 * math.ceil( 2 * bina_score )      # Round to nearest 0.5
            subscore_evidence = 0.5 * round( 2 * subscore_evidence )   # Round to nearest 0.5
            bonus_knowledge   = 0.5 * round( 2 * bonus_knowledge )     # Round to nearest 0.5
            
            if bina_score < 0: bina_score = 0
            
            
            ### FILTER Field String Determinant:
            if bina_score:
                if bina_score >= args.oncoscore_threshold_pass:
                    vcf_i.filters = 'PASS'
    
                elif bina_score >= args.oncoscore_threshold_lowconfidence:
                    vcf_i.filters = 'LowQual'
                    
                else:
                    vcf_i.filters = 'REJECT'
                    
            else:
                vcf_i.filters = 'REJECT'
            
            
            
            # Sort the output lines:
            
            ## RE-FORMAT INFO and SAMPLE (FORMAT) columns:
            new_info_field = []
            
            # Expected to be flags:
            for info_item in ('COMMON', 'CAF', 'G5', 'G5A', 'MUT'):
                
                info_value = vcf_i.get_info_value( info_item )
                
                if info_value == 'true' or info_value == '.':
                    new_info_field.append( info_item )

                elif info_value:
                    new_info_field.append( '{}={}'.format(info_item, info_value) )
            
            # Expected to have values:
            for info_item in ('GENE', 'AA', 'CDS', 'CNT', 'EFF', 'LSEQ', 'RSEQ'):
                
                info_value = vcf_i.get_info_value( info_item )
                
                if info_value == 'true':
                    new_info_field.append( info_item )
                    
                elif info_value == '.':
                    new_info_field.append( '{}={}'.format(info_item, '0') )
                
                elif info_value:
                    new_info_field.append( '{}={}'.format(info_item, info_value) )

            ##
            
            
            # Add entries:
            new_info_field.append( 'NUM_SMMETHODS={}'.format( num_smmethods) )
            
            if num_smmethods >= 1:
                
                sources = []
                for call_i, tool_i in zip(combo_i, tools):
                    
                    if call_i:
                        sources.append( tool_i )
                        
                new_info_field.append( 'SOURCES=' + ','.join(sources) )
    
                new_info_field.append( 'subEvidence=' + str(subscore_evidence) )
                new_info_field.append( 'subKnowledge={},{}'.format(bonus_knowledge, penalty_knowledge) )
            
            
            if bina_score:
                new_info_field.append( 'OncoRK=' + str(bina_score) )
    
    
    
            ### Replacing info field with this one:
            new_info_column = ';'.join(new_info_field)  # Modified INFO FIELD
            
            ### Format the sample info:
            new_format, new_normal_sample, new_tumor_sample = [], [], []
            
            # Genotype first, as always
            if 'GT' in vcf_i.get_sample_item().keys():
                new_format.append( 'GT' )
                new_tumor_sample.append( vcf_i.get_sample_value('GT', idxT) )
                if paired_mode: new_normal_sample.append( vcf_i.get_sample_value('GT', idxN) )
            
            # Followed by DP4, but only if it exists:
            if 'DP4' in vcf_i.get_sample_item().keys():
                new_format.append( 'DP4' )
                new_tumor_sample.append( vcf_i.get_sample_value('DP4', idxT) )
                if paired_mode: new_normal_sample.append( vcf_i.get_sample_value('DP4', idxN) )
                
            elif 'SC' in vcf_i.get_sample_item().keys():
                new_format.append( 'SC' )
                new_tumor_sample.append( vcf_i.get_sample_value('SC', idxT) )
                if paired_mode: new_normal_sample.append( vcf_i.get_sample_value('SC', idxN) )
            
            # If no DP4, get the depth info without strand info:
            else:
                if 'DP' in vcf_i.get_sample_item().keys():
                    new_format.append( 'DP' )
                    new_tumor_sample.append( vcf_i.get_sample_value('DP', idxT) )
                    if paired_mode: new_normal_sample.append( vcf_i.get_sample_value('DP', idxN) )
                    
                if 'AD' in vcf_i.get_sample_item().keys():
                    new_format.append( 'AD' )
                    new_tumor_sample.append( vcf_i.get_sample_value('AD', idxT) )
                    if paired_mode: new_normal_sample.append( vcf_i.get_sample_value('AD', idxN) )
                    
                if 'RD' in vcf_i.get_sample_item().keys():
                    new_format.append( 'RD' )
                    new_tumor_sample.append( vcf_i.get_sample_value('RD', idxT) )
                    if paired_mode: new_normal_sample.append( vcf_i.get_sample_value('RD', idxN) )
            
            # Remaining stuff
            remaining_variables = ('AMQ', 'BQ', 'MQ', 'VAQ', 'NM', 'PMEAN', 'QSTD', 'PSTD', 'QUAL')
            for remaining_i in remaining_variables:
                if remaining_i in vcf_i.get_sample_item().keys():
                    new_format.append( remaining_i )
                    new_tumor_sample.append( vcf_i.get_sample_value(remaining_i, idxT) )
                    if paired_mode: new_normal_sample.append( vcf_i.get_sample_value(remaining_i, idxN) )
            
            
            
            
            ## Move INFO's information to Sample Info:
            info2sample_variables = ('SHIFT3', 'MSI', 'MSILEN', 'SOR')
            for var_i in info2sample_variables:
                
                value_i = vcf_i.get_info_value( var_i )
                
                if value_i:
                    new_format.append( var_i )
                    new_tumor_sample.append( value_i )
                    if paired_mode: new_normal_sample.append( '.' )
            
            
            
            # Get the Somatic Scores to sample information:
            # SomaticSniper:
            if 'SSC' in vcf_i.get_sample_item().keys():
                new_format.append( 'SSSC' )
                new_tumor_sample.append( vcf_i.get_sample_value('SSC', idxT) )
                if paired_mode: new_normal_sample.append( vcf_i.get_sample_value('SSC', idxN) )
                
            # VarScan2:
            vssc = vcf_i.get_info_value('SSC')
            if vssc:
                new_format.append( 'VSSC' )
                new_tumor_sample.append(vssc)
                if paired_mode: new_normal_sample.append('.')
                
            # JointSNVMix2:
            aaab = vcf_i.get_info_value('AAAB')
            aabb = vcf_i.get_info_value('AABB')
            if aaab and aabb:
                aaab = float(aaab)
                aabb = float(aabb)
                jssc = str( aaab + aabb )
                new_format.append( 'JSSC' )
                new_tumor_sample.append(jssc)
                if paired_mode: new_normal_sample.append('.')
                
            # VarDict:
            dssf = vcf_i.get_info_value('SSF')
            if dssf:
                new_format.append( 'DSSF' )
                new_tumor_sample.append(dssf)
                if paired_mode: new_normal_sample.append('.')
            
            # Put together:
            new_format_string = ':'.join(new_format)
            new_tumor_string  = ':'.join(new_tumor_sample)
            if paired_mode: new_normal_string = ':'.join(new_normal_sample)
            
            
            ### Write:
            if paired_mode:
                modified_line = '\t'.join((vcf_i.chromosome, vcf_i.position, vcf_i.identifier, vcf_i.refbase, vcf_i.altbase, vcf_i.qual, vcf_i.filters, new_info_column, new_format_string, new_normal_string, new_tumor_string )) + '\n'
            else:
                modified_line = '\t'.join((vcf_i.chromosome, vcf_i.position, vcf_i.identifier, vcf_i.refbase, vcf_i.altbase, vcf_i.qual, vcf_i.filters, new_info_column, new_format_string, new_tumor_string )) + '\n'
            
            
            vcf_out.write( modified_line )
    
            
            ##### ##### For making histogram ##### #####
            try:
                bina_score_tally[bina_score] = bina_score_tally[bina_score] + 1
            except KeyError:
                bina_score_tally[bina_score] = 1
            
            try:
                subscore_evidence_tally[subscore_evidence] = subscore_evidence_tally[subscore_evidence] + 1
            except KeyError:
                subscore_evidence_tally[subscore_evidence] = 1
                            
            try:
                bonus_knowledge_tally[bonus_knowledge] = bonus_knowledge_tally[bonus_knowledge] + 1
            except KeyError:
                bonus_knowledge_tally[bonus_knowledge] = 1
    
            try:
                penalty_knowledge_tally[penalty_knowledge] = penalty_knowledge_tally[penalty_knowledge] + 1
            except KeyError:
                penalty_knowledge_tally[penalty_knowledge] = 1
    
            try:
                num_methods_tally[num_smmethods] = num_methods_tally[num_smmethods] + 1
            except KeyError:
                num_methods_tally[num_smmethods] = 1
            

        ## Continue:
        line_i = vcf.readline().rstrip()





# Print the histogram onto stderr:
print('BINA_OncoRK:')
for score_i in sorted( bina_score_tally.keys() ):
    print(score_i, bina_score_tally[score_i], sep='\t')

print('\nNUM_SMMETHODS:')
for score_i in sorted( num_methods_tally.keys() ):
    print(score_i, num_methods_tally[score_i], sep='\t')

print('\nsubEvidence:')
for score_i in sorted( subscore_evidence_tally.keys() ):
    print(score_i, subscore_evidence_tally[score_i], sep='\t')

print('\nbonusKnowledge:')
for score_i in sorted( bonus_knowledge_tally.keys() ):
    print(score_i, bonus_knowledge_tally[score_i], sep='\t')

print('\npenaltyKnowledge:')
for score_i in sorted( penalty_knowledge_tally.keys() ):
    print(score_i, penalty_knowledge_tally[score_i], sep='\t')



# Print combinations:
print('')
print('\t'.join(tools), 'Total', 'dbsnp', 'common', 'COSMIC', sep='\t' )
for mvjs_i in sorted( MVJS_combinations.keys() ):
    
    stats_i = MVJS_combinations[mvjs_i]
    
    for i in mvjs_i:
        print(i, end='\t')
        
    for i in stats_i:
        print(i, end='\t')
        
    print('')
