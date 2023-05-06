#!/usr/bin/env python3

# Quick and dirty, hacky way to do it.

import argparse
import itertools
import os
import sys

import pysam
import somaticseq.genomicFileHandler.genomic_file_handlers as genome
from somaticseq.genomicFileHandler.read_info_extractor import position_of_aligned_read

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "-infile", "--input-vcf-file", type=str, help="Input VCF file", required=True
)
parser.add_argument("-bam", "--bam-file", type=str, help="BAM file", required=True)
parser.add_argument(
    "-ref",
    "--genome-reference",
    type=str,
    help=".fasta file to get the ref base",
    required=True,
    default=None,
)
parser.add_argument(
    "-outfile", "--output-vcf-file", type=str, help="Output VCF file", required=True
)
parser.add_argument(
    "-threshold",
    "--phasing-threshold",
    type=int,
    help="How far apart do we try to phase",
    required=False,
    default=1,
)

args = parser.parse_args()

infile = args.input_vcf_file
bam = args.bam_file
ref_fa = args.genome_reference
outfile = args.output_vcf_file
threshold = args.phasing_threshold

with genome.open_textfile(infile) as infile, pysam.AlignmentFile(bam) as bam, open(
    outfile, "w"
) as outfile, pysam.FastaFile(ref_fa) as ref_fa:

    my_line = infile.readline().rstrip()

    while my_line.startswith("##"):
        outfile.write(my_line + "\n")
        my_line = infile.readline().rstrip()

    # This is to read through and copy the #CHROM line
    assert my_line.startswith("#CHROM")
    outfile.write(
        '##INFO=<ID=COORDINATES,Number=.,Type=Integer,Description="Coordinates of the bases">\n'
    )
    outfile.write(
        '##INFO=<ID=PDP,Number=.,Type=Integer,Description="Phased DP, one for reference, and each of the variant calls.">\n'
    )
    outfile.write(my_line + "\n")
    my_line = infile.readline().rstrip()

    # Get into the bulk of the VCF file
    while my_line:
        my_vcf = genome.VcfLine(my_line)

        if len(my_vcf.refbase) == 1:
            my_coordinates = [(my_vcf.chromosome, my_vcf.position)]
            base_options = [[my_vcf.refbase]]
            base_options[0].extend(my_vcf.altbase.split(","))
            vcf_lines = [my_line]
            filter_status = [my_vcf.filters]

        elif len(my_vcf.refbase) == len(my_vcf.altbase):
            pass

        while (my_coordinates[-1][0] == my_vcf.chromosome) and (
            my_vcf.position - my_coordinates[-1][1] <= threshold
        ):

            my_line = infile.readline().rstrip()
            my_vcf = genome.VcfLine(my_line)

            # If the next vcf line is the same coordinate
            if (my_coordinates[-1][0] == my_vcf.chromosome) and (
                my_vcf.position - my_coordinates[-1][1] == 0
            ):
                base_options[-1].extend(my_vcf.altbase.split(","))
                vcf_lines.append(my_line)
            # If the next vcf line is still within the threshold:
            elif (my_coordinates[-1][0] == my_vcf.chromosome) and (
                my_vcf.position - my_coordinates[-1][1] <= threshold
            ):
                # For the "missing" reference bases:
                missing_bases = ref_fa.fetch(
                    my_vcf.chromosome, my_coordinates[-1][1], my_vcf.position - 1
                )
                for i, missing_coordinate_i in enumerate(
                    range(my_coordinates[-1][1] + 1, my_vcf.position)
                ):
                    my_coordinates.append((my_vcf.chromosome, missing_coordinate_i))
                    base_options.append([missing_bases[i]])
                # The next vcf coordinate
                my_coordinates.append((my_vcf.chromosome, my_vcf.position))
                base_options.append([my_vcf.refbase])
                base_options[-1].extend(my_vcf.altbase.split(","))
                vcf_lines.append(my_line)
                filter_status.append(my_vcf.filters)

        # Trigger the procedure:
        # If the number is too large, 2^n possibilities and will run out of memory (hence less than 6)
        if 6 > len(my_coordinates) > 1:
            ref_string = ref_fa.fetch(
                my_coordinates[0][0], my_coordinates[0][1] - 1, my_coordinates[-1][1]
            )
            mnp_list = list(itertools.product(*base_options))

            # This also serves to "unique-fy" duplicate strings from mnp_list
            mnp_tally = {}
            for string_i in mnp_list:
                mnp_tally["".join(string_i)] = 0

            # Grab the reads from the first coordinate to the last coordinate
            reads = bam.fetch(
                my_coordinates[0][0], my_coordinates[0][1] - 1, my_coordinates[-1][1]
            )
            for read_i in reads:
                if not (read_i.is_unmapped or read_i.is_duplicate):
                    mnp_call = ""
                    for coordinate_i in my_coordinates:
                        (
                            code_i,
                            ith_base,
                            base_call_i,
                            indel_length_i,
                            flanking_indel_i,
                        ) = position_of_aligned_read(read_i, coordinate_i[1] - 1)
                        # The position is matched:
                        if (code_i == 0 or code_i == 1 or code_i == 2) and base_call_i:
                            mnp_call = mnp_call + base_call_i
                    if mnp_call in mnp_tally:
                        mnp_tally[mnp_call] += 1

            mnp_calls = []
            mnp_depths = []
            for call_i in mnp_tally:
                if call_i != ref_string and mnp_tally[call_i] != 0:
                    mnp_calls.append(call_i)
                    mnp_depths.append(str(mnp_tally[call_i]))

            ref_depth = str(mnp_tally[ref_string])
            alt_depths = ",".join(mnp_depths)
            pdp = ",".join((ref_depth, alt_depths))

            # phased_coordinates = [ str(item[1]) for item in my_coordinates]
            # info_line = 'PDP={};COORDINATES={}'.format(pdp, ','.join(phased_coordinates) )
            # How to label PASS/LowQual/REJECT for the phased calls:
            status_score = 0
            for i_th_status in filter_status:
                if "PASS" in i_th_status:
                    status_score += 1
                elif "LowQual" in i_th_status:
                    status_score += 0.5

            if status_score > (2 / 3) * len(filter_status):
                filters_i = "PASS"
            elif status_score >= 0.5 * len(filter_status):
                filters_i = "LowQual"
            else:
                filters_i = "REJECT"

            info_line = "PDP={}".format(pdp)
            mnp_line = "{CHR}\t{POS}\t.\t{REF}\t{ALT}\t.\t{FILTER}\t{INFO}\t{FORMAT}\t{SAMPLE}".format(
                CHR=my_coordinates[0][0],
                POS=my_coordinates[0][1],
                REF=ref_string,
                ALT=",".join(mnp_calls),
                FILTER=filters_i,
                INFO=info_line,
                FORMAT="GT",
                SAMPLE="0/1",
            )
            outfile.write(mnp_line + "\n")

        for line_j in vcf_lines:
            outfile.write(line_j + "\n")
