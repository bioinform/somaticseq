#!/usr/bin/env python3

# single-sample only

import argparse
import logging
import math
import os
import re
from copy import copy

import pysam
import scipy.stats as stats
import somaticseq.annotate_caller as annotate_caller
import somaticseq.genomicFileHandler.genomic_file_handlers as genome
import somaticseq.sequencing_features as sequencing_features
from somaticseq.genomicFileHandler.read_info_extractor import (
    genomic_coordinates,
    rescale,
)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger = logging.getLogger(os.path.basename(__file__))
logger.setLevel(logging.DEBUG)
logger.addHandler(ch)

# Header for the output data, created here so I won't have to indent this line:
out_header = "{CHROM}\t\
{POS}\t\
{ID}\t\
{REF}\t\
{ALT}\t\
{if_MuTect}\t\
{if_Strelka}\t\
{if_VarScan2}\t\
{if_VarDict}\t\
{if_LoFreq}\t\
{if_Scalpel}\t\
{VarScan2_Score}\t\
{if_dbsnp}\t\
{COMMON}\t\
{if_COSMIC}\t\
{COSMIC_CNT}\t\
{Consistent_Mates}\t\
{Inconsistent_Mates}\t\
{Seq_Complexity_Span}\t\
{Seq_Complexity_Adj}\t\
{M2_TLOD}\t\
{M2_ECNT}\t\
{MSI}\t\
{MSILEN}\t\
{SHIFT3}\t\
{MaxHomopolymer_Length}\t\
{SiteHomopolymer_Length}\t\
{T_DP}\t\
{tBAM_REF_MQ}\t\
{tBAM_ALT_MQ}\t\
{tBAM_p_MannWhitneyU_MQ}\t\
{tBAM_REF_BQ}\t\
{tBAM_ALT_BQ}\t\
{tBAM_p_MannWhitneyU_BQ}\t\
{tBAM_REF_NM}\t\
{tBAM_ALT_NM}\t\
{tBAM_NM_Diff}\t\
{tBAM_REF_Concordant}\t\
{tBAM_REF_Discordant}\t\
{tBAM_ALT_Concordant}\t\
{tBAM_ALT_Discordant}\t\
{tBAM_Concordance_FET}\t\
{T_REF_FOR}\t\
{T_REF_REV}\t\
{T_ALT_FOR}\t\
{T_ALT_REV}\t\
{tBAM_StrandBias_FET}\t\
{tBAM_p_MannWhitneyU_EndPos}\t\
{tBAM_REF_Clipped_Reads}\t\
{tBAM_ALT_Clipped_Reads}\t\
{tBAM_Clipping_FET}\t\
{tBAM_MQ0}\t\
{tBAM_Other_Reads}\t\
{tBAM_Poor_Reads}\t\
{tBAM_REF_InDel_3bp}\t\
{tBAM_REF_InDel_2bp}\t\
{tBAM_REF_InDel_1bp}\t\
{tBAM_ALT_InDel_3bp}\t\
{tBAM_ALT_InDel_2bp}\t\
{tBAM_ALT_InDel_1bp}\t\
{InDel_Length}"

extra_caller_header = ""
label_header = "{TrueVariant_or_False}"


def run():

    inputParameters = {}
    parser = argparse.ArgumentParser(
        description="This is a SomaticSeq subroutine to convert a VCF file into a TSV file with all the SomaticSeq features for tumor-only modes. Any VCF file can be used as the main input. The output will have the same variants. Also required is the BAM files, with additional optional inputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    input_sites = parser.add_mutually_exclusive_group()
    input_sites.add_argument(
        "-myvcf",
        "--vcf-format",
        type=str,
        help="Input file is VCF formatted.",
        required=False,
        default=None,
    )
    input_sites.add_argument(
        "-mybed",
        "--bed-format",
        type=str,
        help="Input file is BED formatted.",
        required=False,
        default=None,
    )
    input_sites.add_argument(
        "-mypos",
        "--positions-list",
        type=str,
        help="A list of positions: tab seperating contig and positions.",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-bam",
        "--in-bam",
        type=str,
        help="Tumor tBAM File",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-truth",
        "--ground-truth-vcf",
        type=str,
        help="VCF of true hits",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-dbsnp",
        "--dbsnp-vcf",
        type=str,
        help="dbSNP VCF: do not use if input VCF is annotated",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-cosmic",
        "--cosmic-vcf",
        type=str,
        help="COSMIC VCF: do not use if input VCF is annotated",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-mutect",
        "--mutect-vcf",
        type=str,
        help="MuTect VCF",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-varscan",
        "--varscan-vcf",
        type=str,
        help="VarScan2 VCF",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-vardict",
        "--vardict-vcf",
        type=str,
        help="VarDict VCF",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-lofreq",
        "--lofreq-vcf",
        type=str,
        help="LoFreq VCF",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-scalpel",
        "--scalpel-vcf",
        type=str,
        help="Scalpel VCF",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-strelka",
        "--strelka-vcf",
        type=str,
        help="Strelka VCF",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-arbvcfs",
        "--arbitrary-vcfs",
        type=str,
        help="Arbitrary extra VCFs",
        nargs="*",
        default=[],
    )
    parser.add_argument(
        "-ref",
        "--genome-reference",
        type=str,
        help=".fasta.fai file to get the contigs",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-dedup",
        "--deduplicate",
        action="store_true",
        help="Do not consider duplicate reads from tBAM files. Default is to count everything",
        required=False,
        default=False,
    )
    parser.add_argument(
        "-minMQ",
        "--minimum-mapping-quality",
        type=float,
        help="Minimum mapping quality below which is considered poor",
        required=False,
        default=1,
    )
    parser.add_argument(
        "-minBQ",
        "--minimum-base-quality",
        type=float,
        help="Minimum base quality below which is considered poor",
        required=False,
        default=5,
    )
    parser.add_argument(
        "-mincaller",
        "--minimum-num-callers",
        type=float,
        help="Minimum number of tools to be considered",
        required=False,
        default=0,
    )
    parser.add_argument(
        "-scale",
        "--p-scale",
        type=str,
        help="phred, fraction, or none",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-outfile",
        "--output-tsv-file",
        type=str,
        help="Output TSV Name",
        required=False,
        default=os.sys.stdout,
    )
    args = parser.parse_args()
    inputParameters = vars(args)
    return inputParameters


def vcf2tsv(
    is_vcf=None,
    is_bed=None,
    is_pos=None,
    bam_fn=None,
    truth=None,
    cosmic=None,
    dbsnp=None,
    mutect=None,
    varscan=None,
    vardict=None,
    lofreq=None,
    scalpel=None,
    strelka=None,
    arbitrary_vcfs=None,
    dedup=True,
    min_mq=1,
    min_bq=5,
    min_caller=0,
    ref_fa=None,
    p_scale=None,
    outfile=None,
):

    if arbitrary_vcfs is None:
        arbitrary_vcfs = []

    # Convert contig_sequence to chrom_seq dict:
    fai_file = ref_fa + ".fai"
    chrom_seq = genome.faiordict2contigorder(fai_file, "fai")

    # Determine input format:
    if is_vcf:
        mysites = is_vcf
    elif is_bed:
        mysites = is_bed
    elif is_pos:
        mysites = is_pos
    else:
        mysites = fai_file
        logger.info("No position supplied. Will evaluate the whole genome.")

    # Re-scale output or not:
    if p_scale == None:
        logger.info("NO RE-SCALING")
    elif p_scale.lower() == "phred":
        p_scale = "phred"
    elif p_scale.lower() == "fraction":
        p_scale = "fraction"
    else:
        p_scale = None
        logger.info("NO RE-SCALING")

    # Define NaN and Inf:
    nan = float("nan")
    inf = float("inf")
    pattern_chr_position = genome.pattern_chr_position

    ## Running
    with genome.open_textfile(mysites) as my_sites, open(outfile, "w") as outhandle:

        my_line = my_sites.readline().rstrip()
        bam = pysam.AlignmentFile(bam_fn, reference_filename=ref_fa)
        ref_fa = pysam.FastaFile(ref_fa)

        if truth:
            truth = genome.open_textfile(truth)
            truth_line = genome.skip_vcf_header(truth)
        if cosmic:
            cosmic = genome.open_textfile(cosmic)
            cosmic_line = genome.skip_vcf_header(cosmic)
        if dbsnp:
            dbsnp = genome.open_textfile(dbsnp)
            dbsnp_line = genome.skip_vcf_header(dbsnp)
        # 6 Incorporate callers: get thru the #'s
        if mutect:
            mutect = genome.open_textfile(mutect)
            mutect_line = genome.skip_vcf_header(mutect)
        if varscan:
            varscan = genome.open_textfile(varscan)
            varscan_line = genome.skip_vcf_header(varscan)
        if vardict:
            vardict = genome.open_textfile(vardict)
            vardict_line = genome.skip_vcf_header(vardict)
        if lofreq:
            lofreq = genome.open_textfile(lofreq)
            lofreq_line = genome.skip_vcf_header(lofreq)
        if scalpel:
            scalpel = genome.open_textfile(scalpel)
            scalpel_line = genome.skip_vcf_header(scalpel)
        if strelka:
            strelka = genome.open_textfile(strelka)
            strelka_line = genome.skip_vcf_header(strelka)

        arbitrary_file_handle = {}
        arbitrary_line = {}
        for ith_arbi, arbitrary_vcf_i in enumerate(arbitrary_vcfs):
            arbitrary_file_handle[ith_arbi] = genome.open_textfile(arbitrary_vcf_i)
            arbitrary_line[ith_arbi] = genome.skip_vcf_header(
                arbitrary_file_handle[ith_arbi]
            )

        # Get through all the headers:
        while my_line.startswith("#") or my_line.startswith("track="):
            my_line = my_sites.readline().rstrip()

        # First coordinate, for later purpose of making sure the input is sorted properly
        coordinate_i = re.match(genome.pattern_chr_position, my_line)
        coordinate_i = coordinate_i.group() if coordinate_i else ""

        # First line:
        header_part_1 = out_header.replace("{", "").replace("}", "")

        additional_arbi_caller_numbers = sorted(arbitrary_file_handle.keys())
        for arbi_caller_num in additional_arbi_caller_numbers:
            header_part_1 = (
                header_part_1 + "\t" + "if_Caller_{}".format(arbi_caller_num)
            )
        header_last_part = label_header.replace("{", "").replace("}", "")
        outhandle.write("\t".join((header_part_1, header_last_part)) + "\n")
        while my_line:
            # If VCF, get all the variants with the same coordinate into a list:
            if is_vcf:
                my_vcf = genome.VcfLine(my_line)
                my_coordinates = [(my_vcf.chromosome, my_vcf.position)]
                variants_at_my_coordinate = []
                alt_bases = my_vcf.altbase.split(",")
                for alt_i in alt_bases:
                    vcf_i = copy(my_vcf)
                    vcf_i.altbase = alt_i
                    variants_at_my_coordinate.append(vcf_i)

                # As long as the "coordinate" stays the same, it will keep reading until it's different.
                while my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):
                    my_line = my_sites.readline().rstrip()
                    my_vcf = genome.VcfLine(my_line)

                    ########## This block is code is to ensure the input VCF file is properly sorted ##
                    coordinate_j = re.match(genome.pattern_chr_position, my_line)
                    coordinate_j = coordinate_j.group() if coordinate_j else ""
                    if genome.whoisbehind(coordinate_i, coordinate_j, chrom_seq) == 1:
                        raise Exception(
                            "{} does not seem to be properly sorted.".format(mysites)
                        )

                    coordinate_i = coordinate_j
                    ###################################################################################
                    if my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):
                        alt_bases = my_vcf.altbase.split(",")
                        for alt_i in alt_bases:

                            vcf_i = copy(my_vcf)
                            vcf_i.altbase = alt_i
                            variants_at_my_coordinate.append(vcf_i)
            elif is_bed:
                bed_item = my_line.split("\t")
                my_coordinates = genomic_coordinates(
                    bed_item[0], int(bed_item[1]) + 1, int(bed_item[2])
                )
            elif is_pos:
                pos_item = my_line.split("\t")
                my_coordinates = genomic_coordinates(
                    pos_item[0], int(pos_item[1]), int(pos_item[1])
                )
            elif fai_file:
                fai_item = my_line.split("\t")
                my_coordinates = genomic_coordinates(fai_item[0], 1, int(fai_item[1]))

            ##### ##### ##### ##### ##### #####
            for my_coordinate in my_coordinates:
                ######## If VCF, can get ref base, variant base, as well as other identifying information ########
                if is_vcf:
                    ref_bases = []
                    alt_bases = []
                    indel_lengths = []
                    all_my_identifiers = []
                    for variant_i in variants_at_my_coordinate:
                        ref_base = variant_i.refbase
                        first_alt = variant_i.altbase.split(",")[0]
                        indel_length = len(first_alt) - len(ref_base)
                        ref_bases.append(ref_base)
                        alt_bases.append(first_alt)
                        indel_lengths.append(indel_length)

                        # Extract these information if they exist in the VCF file, but they could be re-written if dbSNP/COSMIC are supplied.
                        if_dbsnp = (
                            1 if re.search(r"rs[0-9]+", variant_i.identifier) else 0
                        )
                        if_cosmic = (
                            1
                            if re.search(r"COS[MN][0-9]+", variant_i.identifier)
                            else 0
                        )
                        if_common = (
                            1 if variant_i.get_info_value("COMMON") == "1" else 0
                        )
                        num_cases = (
                            variant_i.get_info_value("CNT")
                            if variant_i.get_info_value("CNT")
                            else nan
                        )
                        if variant_i.identifier == ".":
                            my_identifier_i = set()
                        else:
                            my_identifier_i = variant_i.identifier.split(";")
                            my_identifier_i = set(my_identifier_i)

                        all_my_identifiers.append(my_identifier_i)

                ## If not, 1) get ref_base, first_alt from other VCF files.
                #          2) Create placeholders for dbSNP and COSMIC that can be overwritten with dbSNP/COSMIC VCF files (if provided)
                else:
                    variants_at_my_coordinate = [
                        None
                    ]  # Just to have something to iterate
                    ref_base = first_alt = indel_length = None

                    # Could be re-written if dbSNP/COSMIC are supplied. If not, they will remain NaN.
                    if_dbsnp = if_cosmic = if_common = num_cases = nan

                #################################### Find the same coordinate in those VCF files ####################################
                if mutect:
                    (
                        got_mutect,
                        mutect_variants,
                        mutect_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, mutect_line, mutect, chrom_seq
                    )
                if varscan:
                    (
                        got_varscan,
                        varscan_variants,
                        varscan_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, varscan_line, varscan, chrom_seq
                    )
                if vardict:
                    (
                        got_vardict,
                        vardict_variants,
                        vardict_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, vardict_line, vardict, chrom_seq
                    )
                if lofreq:
                    (
                        got_lofreq,
                        lofreq_variants,
                        lofreq_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, lofreq_line, lofreq, chrom_seq
                    )
                if scalpel:
                    (
                        got_scalpel,
                        scalpel_variants,
                        scalpel_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, scalpel_line, scalpel, chrom_seq
                    )
                if strelka:
                    (
                        got_strelka,
                        strelka_variants,
                        strelka_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, strelka_line, strelka, chrom_seq
                    )
                if truth:
                    (
                        got_truth,
                        truth_variants,
                        truth_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, truth_line, truth, chrom_seq
                    )
                if dbsnp:
                    (
                        got_dbsnp,
                        dbsnp_variants,
                        dbsnp_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, dbsnp_line, dbsnp, chrom_seq
                    )
                if cosmic:
                    (
                        got_cosmic,
                        cosmic_variants,
                        cosmic_line,
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate, cosmic_line, cosmic, chrom_seq
                    )

                got_arbitraries = {}
                arbitrary_variants = {}
                for ith_arbi in arbitrary_file_handle:
                    (
                        got_arbitraries[ith_arbi],
                        arbitrary_variants[ith_arbi],
                        arbitrary_line[ith_arbi],
                    ) = genome.find_vcf_at_coordinate(
                        my_coordinate,
                        arbitrary_line[ith_arbi],
                        arbitrary_file_handle[ith_arbi],
                        chrom_seq,
                    )

                # Now, use pysam to look into the tBAM file(s), variant by variant from the input:
                for ith_call, my_call in enumerate(variants_at_my_coordinate):
                    if is_vcf:
                        # The particular line in the input VCF file:
                        variant_id = (
                            (my_call.chromosome, my_call.position),
                            my_call.refbase,
                            my_call.altbase,
                        )
                        ref_base = ref_bases[ith_call]
                        first_alt = alt_bases[ith_call]
                        indel_length = indel_lengths[ith_call]
                        my_identifiers = all_my_identifiers[ith_call]
                    else:
                        variant_id = (
                            (my_coordinate[0], my_coordinate[1]),
                            ref_base,
                            first_alt,
                        )
                    # Reset num_caller to 0 for each variant in the same coordinate
                    num_callers = 0

                    #################### Collect Caller Vcf ####################:
                    if mutect:
                        mutect_classification, tlod, ecnt = annotate_caller.ssMuTect(
                            variant_id, mutect_variants
                        )
                        num_callers += mutect_classification
                    else:
                        mutect_classification = tlod = ecnt = nan

                    if varscan:
                        (
                            varscan_classification,
                            score_varscan2,
                        ) = annotate_caller.ssVarScan(variant_id, varscan_variants)
                        num_callers += varscan_classification
                    else:
                        varscan_classification = score_varscan2 = nan

                    if vardict:
                        (
                            vardict_classification,
                            msi,
                            msilen,
                            shift3,
                            t_pmean,
                            t_pstd,
                            t_qstd,
                        ) = annotate_caller.ssVarDict(variant_id, vardict_variants)
                        num_callers += vardict_classification
                    else:
                        vardict_classification = (
                            msi
                        ) = msilen = shift3 = t_pmean = t_pstd = t_qstd = nan

                    if lofreq:
                        lofreq_classification = annotate_caller.ssLoFreq(
                            variant_id, lofreq_variants
                        )
                        num_callers += lofreq_classification
                    else:
                        lofreq_classification = nan

                    if scalpel:
                        scalpel_classification = annotate_caller.ssScalpel(
                            variant_id, scalpel_variants
                        )
                        num_callers += scalpel_classification
                    else:
                        scalpel_classification = nan

                    if strelka:
                        strelka_classification = annotate_caller.ssStrelka(
                            variant_id, strelka_variants
                        )
                        num_callers += strelka_classification
                    else:
                        strelka_classification = nan

                    arbitrary_classifications = {}
                    for ith_arbi_var in arbitrary_file_handle:
                        arbi_classification_i = annotate_caller.anyInputVcf(
                            variant_id, arbitrary_variants[ith_arbi_var]
                        )
                        arbitrary_classifications[ith_arbi_var] = arbi_classification_i
                        num_callers += arbi_classification_i

                    # Potentially write the output only if it meets this threshold:
                    if num_callers >= min_caller:

                        ########## Ground truth file ##########
                        if truth:
                            if variant_id in truth_variants.keys():
                                judgement = 1
                                my_identifiers.add("TruePositive")
                            else:
                                judgement = 0
                                my_identifiers.add("FalsePositive")
                        else:
                            judgement = nan

                        ########## dbSNP ########## Will overwrite dbSNP info from input VCF file
                        if dbsnp:
                            if_dbsnp, if_common, rsID = annotate_caller.dbSNP(
                                variant_id, dbsnp_variants
                            )
                            for ID_i in rsID:
                                my_identifiers.add(ID_i)

                        ########## COSMIC ########## Will overwrite COSMIC info from input VCF file
                        if cosmic:
                            if_cosmic, num_cases, cosmicID = annotate_caller.COSMIC(
                                variant_id, cosmic_variants
                            )
                            for ID_i in cosmicID:
                                my_identifiers.add(ID_i)

                        ########## ######### INFO EXTRACTION FROM BAM FILES ########## #########
                        # Tumor tBAM file:
                        tBamFeatures = sequencing_features.from_bam(
                            bam, my_coordinate, ref_base, first_alt, min_mq, min_bq
                        )

                        # Homopolymer eval:
                        (
                            homopolymer_length,
                            site_homopolymer_length,
                        ) = sequencing_features.from_genome_reference(
                            ref_fa, my_coordinate, ref_base, first_alt
                        )

                        # Linguistic sequence complexity in a +/-80bp window, but substring calculation stops at 20-bp substring.
                        seq_span_80bp = ref_fa.fetch(
                            my_coordinate[0],
                            max(0, my_coordinate[1] - 41),
                            my_coordinate[1] + 40,
                        )
                        seq_left_80bp = ref_fa.fetch(
                            my_coordinate[0],
                            max(0, my_coordinate[1] - 81),
                            my_coordinate[1],
                        )
                        seq_right_80bp = ref_fa.fetch(
                            my_coordinate[0], my_coordinate[1], my_coordinate[1] + 81
                        )

                        if len(seq_span_80bp) > 20:
                            LC_spanning = sequencing_features.subLC(seq_span_80bp, 20)
                        else:
                            LC_spanning = math.nan

                        if len(seq_left_80bp) > 20:
                            left_LC = sequencing_features.subLC(seq_left_80bp, 20)
                        else:
                            left_LC = math.nan

                        if len(seq_right_80bp) > 20:
                            right_LC = sequencing_features.subLC(seq_right_80bp, 20)
                        else:
                            right_LC = math.nan

                        LC_adjacent = min(left_LC, right_LC)
                        LC_spanning_phred = genome.p2phred(1 - LC_spanning, 40)
                        LC_adjacent_phred = genome.p2phred(1 - LC_adjacent, 40)

                        # Fill the ID field of the TSV/VCF
                        my_identifiers = (
                            ";".join(my_identifiers) if my_identifiers else "."
                        )
                        ###
                        out_line_part_1 = out_header.format(
                            CHROM=my_coordinate[0],
                            POS=my_coordinate[1],
                            ID=my_identifiers,
                            REF=ref_base,
                            ALT=first_alt,
                            if_MuTect=mutect_classification,
                            if_Strelka=strelka_classification,
                            if_VarScan2=varscan_classification,
                            if_VarDict=vardict_classification,
                            if_LoFreq=lofreq_classification,
                            if_Scalpel=scalpel_classification,
                            VarScan2_Score=rescale(
                                score_varscan2, "phred", p_scale, 1001
                            ),
                            if_dbsnp=if_dbsnp,
                            COMMON=if_common,
                            if_COSMIC=if_cosmic,
                            COSMIC_CNT=num_cases,
                            Consistent_Mates=tBamFeatures["consistent_mates"],
                            Inconsistent_Mates=tBamFeatures["inconsistent_mates"],
                            Seq_Complexity_Span=LC_spanning_phred,
                            Seq_Complexity_Adj=LC_adjacent_phred,
                            M2_TLOD=tlod,
                            M2_ECNT=ecnt,
                            MSI=msi,
                            MSILEN=msilen,
                            SHIFT3=shift3,
                            MaxHomopolymer_Length=homopolymer_length,
                            SiteHomopolymer_Length=site_homopolymer_length,
                            T_DP=tBamFeatures["dp"],
                            tBAM_REF_MQ="%g" % tBamFeatures["ref_mq"],
                            tBAM_ALT_MQ="%g" % tBamFeatures["alt_mq"],
                            tBAM_p_MannWhitneyU_MQ="%g"
                            % tBamFeatures["p_mannwhitneyu_mq"],
                            tBAM_REF_BQ="%g" % tBamFeatures["ref_bq"],
                            tBAM_ALT_BQ="%g" % tBamFeatures["alt_bq"],
                            tBAM_p_MannWhitneyU_BQ="%g"
                            % tBamFeatures["p_mannwhitneyu_bq"],
                            tBAM_REF_NM="%g" % tBamFeatures["ref_NM"],
                            tBAM_ALT_NM="%g" % tBamFeatures["alt_NM"],
                            tBAM_NM_Diff="%g" % tBamFeatures["NM_Diff"],
                            tBAM_REF_Concordant=tBamFeatures["ref_concordant_reads"],
                            tBAM_REF_Discordant=tBamFeatures["ref_discordant_reads"],
                            tBAM_ALT_Concordant=tBamFeatures["alt_concordant_reads"],
                            tBAM_ALT_Discordant=tBamFeatures["alt_discordant_reads"],
                            tBAM_Concordance_FET=rescale(
                                tBamFeatures["concordance_fet"],
                                "fraction",
                                p_scale,
                                1001,
                            ),
                            T_REF_FOR=tBamFeatures["ref_for"],
                            T_REF_REV=tBamFeatures["ref_rev"],
                            T_ALT_FOR=tBamFeatures["alt_for"],
                            T_ALT_REV=tBamFeatures["alt_rev"],
                            tBAM_StrandBias_FET=rescale(
                                tBamFeatures["strandbias_fet"],
                                "fraction",
                                p_scale,
                                1001,
                            ),
                            tBAM_p_MannWhitneyU_EndPos="%g"
                            % tBamFeatures["p_mannwhitneyu_endpos"],
                            tBAM_REF_Clipped_Reads=tBamFeatures["ref_SC_reads"],
                            tBAM_ALT_Clipped_Reads=tBamFeatures["alt_SC_reads"],
                            tBAM_Clipping_FET=rescale(
                                tBamFeatures["clipping_fet"], "fraction", p_scale, 1001
                            ),
                            tBAM_MQ0=tBamFeatures["MQ0"],
                            tBAM_Other_Reads=tBamFeatures["noise_read_count"],
                            tBAM_Poor_Reads=tBamFeatures["poor_read_count"],
                            tBAM_REF_InDel_3bp=tBamFeatures["ref_indel_3bp"],
                            tBAM_REF_InDel_2bp=tBamFeatures["ref_indel_2bp"],
                            tBAM_REF_InDel_1bp=tBamFeatures["ref_indel_1bp"],
                            tBAM_ALT_InDel_3bp=tBamFeatures["alt_indel_3bp"],
                            tBAM_ALT_InDel_2bp=tBamFeatures["alt_indel_2bp"],
                            tBAM_ALT_InDel_1bp=tBamFeatures["alt_indel_1bp"],
                            InDel_Length=indel_length,
                        )

                        additional_caller_columns = []
                        for arbi_key_i in additional_arbi_caller_numbers:
                            additional_caller_columns.append(
                                str(arbitrary_classifications[arbi_key_i])
                            )
                        additional_caller_columns = "\t".join(additional_caller_columns)

                        label_column = label_header.format(
                            TrueVariant_or_False=judgement
                        )

                        if len(additional_arbi_caller_numbers) > 0:
                            out_line = "\t".join(
                                (
                                    out_line_part_1,
                                    additional_caller_columns,
                                    label_column,
                                )
                            )
                        else:
                            out_line = "\t".join((out_line_part_1, label_column))

                        # Print it out to stdout:
                        outhandle.write(out_line + "\n")

            # Read into the next line:
            if not is_vcf:
                my_line = my_sites.readline().rstrip()

        ##########  Close all open files if they were opened  ##########
        opened_files = [
            ref_fa,
            bam,
            truth,
            cosmic,
            dbsnp,
            mutect,
            varscan,
            vardict,
            lofreq,
            scalpel,
            strelka,
        ]
        [
            opened_files.append(extra_opened_file)
            for extra_opened_file in arbitrary_file_handle.values()
        ]
        [opened_file.close() for opened_file in opened_files if opened_file]


if __name__ == "__main__":
    runParameters = run()

    vcf2tsv(
        is_vcf=runParameters["vcf_format"],
        is_bed=runParameters["bed_format"],
        is_pos=runParameters["positions_list"],
        bam_fn=runParameters["in_bam"],
        truth=runParameters["ground_truth_vcf"],
        cosmic=runParameters["cosmic_vcf"],
        dbsnp=runParameters["dbsnp_vcf"],
        mutect=runParameters["mutect_vcf"],
        varscan=runParameters["varscan_vcf"],
        vardict=runParameters["vardict_vcf"],
        lofreq=runParameters["lofreq_vcf"],
        scalpel=runParameters["scalpel_vcf"],
        strelka=runParameters["strelka_vcf"],
        arbitrary_vcfs=runParameters["arbitrary_vcfs"],
        dedup=runParameters["deduplicate"],
        min_mq=runParameters["minimum_mapping_quality"],
        min_bq=runParameters["minimum_base_quality"],
        min_caller=runParameters["minimum_num_callers"],
        ref_fa=runParameters["genome_reference"],
        p_scale=runParameters["p_scale"],
        outfile=runParameters["output_tsv_file"],
    )
