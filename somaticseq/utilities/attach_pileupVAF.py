#!/usr/bin/env python3
# Supports Insertion/Deletion as well as SNVs
# Last updated: 8/29/2015

import argparse
import gzip
import math
import os
import re
import sys

import somaticseq.genomicFileHandler.genomic_file_handlers as genome
import somaticseq.genomicFileHandler.pileup_reader as pileup

nan = float("nan")
inf = float("inf")

parser = argparse.ArgumentParser(
    description="Given either a tumor-only or tumor-normal VCF file (requires SAMPLE NAME specified), and pileup file, it will attach VAF calculated from pileup file to the VCF file. The pileup file can also be streamed in.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-myvcf", "--my-vcf-file", type=str, help="My VCF", required=True, default=None
)

parser.add_argument(
    "-normal",
    "--normal-sample-name",
    type=str,
    help="Normal Sample Name",
    required=False,
    default="NORMAL",
)
parser.add_argument(
    "-tumor",
    "--tumor-sample-name",
    type=str,
    help="Tumor Sample Name",
    required=False,
    default="TUMOR",
)
parser.add_argument(
    "-Npileup",
    "--normal-pileup-file",
    type=str,
    help="Normal VCF File",
    required=False,
    default=None,
)
parser.add_argument(
    "-Tpileup", "--tumor-pileup-file", type=str, help="Tumor VCF File", required=True
)
parser.add_argument(
    "-fai",
    "--reference-fasta-fai",
    type=str,
    help="Use the fasta.fai file to get the valid contigs",
    required=False,
    default=None,
)
parser.add_argument(
    "-dict",
    "--reference-fasta-dict",
    type=str,
    help="Use the reference dict file to get the valid contigs",
    required=False,
    default=None,
)

# From pileup:
parser.add_argument(
    "-plVAF",
    "--pileup-variant-allele-frequency",
    action="store_true",
    help="Variant Allele Frequency calculated from pileup file",
    required=False,
)
parser.add_argument(
    "-plDP4",
    "--pileup-DP4",
    action="store_true",
    help="DP4 from pileup file",
    required=False,
)

# output file
parser.add_argument(
    "-outfile", "--output-file", type=str, help="Output File Name", required=True
)

args = parser.parse_args()


##
my_vcf = args.my_vcf_file
Tpileup = args.tumor_pileup_file
Npileup = args.normal_pileup_file
tumor_name = args.tumor_sample_name
normal_name = args.normal_sample_name
fai_file = args.reference_fasta_fai
dict_file = args.reference_fasta_dict
outfile = args.output_file

nan = float("nan")

#### Append headers according to user selection ####
header_append = []
format_append = []

if args.pileup_DP4:
    header_append.append(
        '##FORMAT=<ID=plDP4,Number=4,Type=Integer,Description="DP4 from pileup: ref forward, ref reverse, alt forward, alt reverse">'
    )
    format_append.append("plDP4")

if args.pileup_variant_allele_frequency:
    header_append.append(
        '##FORMAT=<ID=plVAF,Number=1,Type=Float,Description="Variant allele frequency calculated from pileup">'
    )
    format_append.append("plVAF")


# Start Working by opening files:
try:
    my_vcf = genome.open_textfile(my_vcf)
    Tpileup = genome.open_textfile(Tpileup)
    outhandle = open(outfile, "w")
    Npileup = genome.open_textfile(Npileup)
except AttributeError:
    pass

if Npileup:
    npileup_line = Npileup.readline().rstrip("\n")

if Tpileup:
    tpileup_line = Tpileup.readline().rstrip("\n")

# Add the extra headers:
out_vcf_headers = genome.vcf_header_modifier(my_vcf, addons=header_append)

# Find out where the tumor and normal samples are in the vcf files, i.e., which column.
# Then, Assuming there are two sample names in "my vcf," the one that appears first should have an index of 0, and the next one is 1:
main_header = out_vcf_headers[3].split("\t")
vcf_idxT = main_header.index(tumor_name)
idxT = vcf_idxT - 9

try:
    vcf_idxN = main_header.index(normal_name)
    idxN = vcf_idxN - 9
except ValueError:
    vcf_idxN = None
    idxN = None


# Write the headers to the output vcf file:
outhandle.write(out_vcf_headers[0] + "\n")  ##fileformat=VCFv4.1
[outhandle.write(out_i + "\n") for out_i in out_vcf_headers[1]]
[outhandle.write(out_i + "\n") for out_i in out_vcf_headers[2]]
outhandle.write(out_vcf_headers[3] + "\n")  # CHROM...


# Convert contig_sequence to chrom_seq dict:
if dict_file:
    chrom_seq = genome.faiordict2contigorder(dict_file, "dict")
elif fai_file:
    chrom_seq = genome.faiordict2contigorder(fai_file, "fai")
else:
    raise Exception(
        "I need a fai or dict file, or else I do not know the contig order."
    )


pattern_chrom = r"|".join(chrom_seq)
r_chrom = r"(" + pattern_chrom + r")"
pattern_chr_position = r_chrom + r"\t[0-9]+"


# Figure out the order of NORMAL and TUMOR
if idxN != None:
    if Npileup and idxN == 0:
        external_pileups = [[Npileup, Tpileup], [npileup_line, tpileup_line]]
    elif Npileup and idx == 1:
        external_pileups = [[Tpileup, Npileup], [tpileup_line, npileup_line]]
    elif not Npileup:
        external_pileups = [[Tpileup], [tpileup_line]]
else:
    external_pileups = [[Tpileup], [tpileup_line]]


line_i = my_vcf.readline().rstrip("\n")
while line_i:

    my_coordinate = re.search(pattern_chr_position, line_i)
    if my_coordinate:
        my_coordinate = my_coordinate.group()
    else:
        print(line_i, file=sys.stderr)
        raise Exception("Your VCF file has a contig that does not exist.")

    # my_vcf:
    vcf_i = genome.VcfLine(line_i)

    # Modify the FORMAT column:
    field_items = vcf_i.get_sample_variable()
    field_items.extend(format_append)
    field_format_line = ":".join(field_items)

    ###########################################################################################
    ###################### Find the same coordinate in the pileup file ########################
    # Line up the order of reading the two files the same order as the sample columns in my_vcf:
    samples_collect = []
    for SM_idx, current_vcf in enumerate(external_pileups[0]):

        latest_pileup_run = genome.catchup(
            my_coordinate, external_pileups[1][SM_idx], current_vcf, chrom_seq
        )
        latest_sample = pileup.Pileup_line(latest_pileup_run[1])

        sample_append = []

        # If the position exists in this samtools generated vcf file:
        if latest_pileup_run[0]:

            assert vcf_i.position == latest_sample.position

            # Figure out alternate pattern:
            first_alt_call = vcf_i.altbase.split(",")[0]

            base_calls = latest_sample.base_reads()

            if base_calls:
                # SNV
                if len(first_alt_call) == len(vcf_i.refbase):

                    ref_for, ref_rev, alt_for, alt_rev = (
                        base_calls[0],
                        base_calls[1],
                        base_calls[2].count(first_alt_call.upper()),
                        base_calls[3].count(first_alt_call.lower()),
                    )

                # Insertion:
                elif len(first_alt_call) > len(vcf_i.refbase):

                    inserted = first_alt_call[len(vcf_i.refbase) : :]

                    ref_for, ref_rev, alt_for, alt_rev = (
                        base_calls[0],
                        base_calls[1],
                        base_calls[6].count(inserted.upper()),
                        base_calls[7].count(inserted.lower()),
                    )

                # Deletion:
                elif len(first_alt_call) < len(vcf_i.refbase):

                    deleted = vcf_i.refbase[len(first_alt_call) : :]

                    ref_for, ref_rev, alt_for, alt_rev = (
                        base_calls[0],
                        base_calls[1],
                        base_calls[4].count(deleted.upper()),
                        base_calls[5].count(deleted.lower()),
                    )

            else:
                ref_for = ref_rev = alt_for = alt_rev = 0

            ### Pre-defined material ###
            ### If user wants DP4 ###
            if args.pileup_DP4:

                pl_DP4 = "{},{},{},{}".format(ref_for, ref_rev, alt_for, alt_rev)
                sample_append.append(pl_DP4)

            ### If user wants VAF ###
            if args.pileup_variant_allele_frequency:

                try:
                    pl_vaf = (alt_for + alt_rev) / (
                        alt_for + alt_rev + ref_for + ref_rev
                    )
                except ZeroDivisionError:
                    pl_vaf = 0

                pl_vaf = "%.3g" % pl_vaf
                sample_append.append(pl_vaf)

            # Reset the current line:
            sample_items = list(vcf_i.get_sample_item(idx=SM_idx, out_type="l")[1])
            sample_items.extend(sample_append)
            sample_out = ":".join(sample_items)

            # Reset the current line:
            external_pileups[1][SM_idx] = latest_sample.pileup_line

            # New format and sample columns:
            samples_collect.append(sample_out)

        # If the position does not exist in pileup file:
        else:

            sample_items = list(vcf_i.get_sample_item(idx=SM_idx, out_type="l")[1])
            sample_append = ["." if i != "plDP4" else ".,.,.,." for i in format_append]
            sample_items.extend(sample_append)

            sample_out = ":".join(sample_items)

            samples_collect.append(sample_out)
            external_pileups[1][SM_idx] = latest_sample.pileup_line

    ### Write out will have a few different possible situations ###
    # If NORMAL and TUMOR both exist in the designated VCF file:
    if vcf_idxT and vcf_idxN:

        # But the Nvcf is not supplied, modified the NORMAL column to reflect change in FORMAT column:
        if not Npileup:
            normal_items = list(vcf_i.get_sample_item(idx=idxN, out_type="l")[1])
            extra_normal_items = [
                "." if i != "plDP4" else ".,.,.,." for i in format_append
            ]
            normal_out = ":".join(extra_normal_items)
            samples_collect.append(normal_out)

        # Write out:
        out_i = "\t".join(
            (
                vcf_i.chromosome,
                str(vcf_i.position),
                vcf_i.identifier,
                vcf_i.refbase,
                vcf_i.altbase,
                vcf_i.qual,
                vcf_i.filters,
                vcf_i.info,
                field_format_line,
                samples_collect[0],
                samples_collect[1],
            )
        )
        outhandle.write(out_i + "\n")

    # Only TUMOR exists in the designated VCF file:
    if not vcf_idxN:

        # Write out:
        out_i = "\t".join(
            (
                vcf_i.chromosome,
                str(vcf_i.position),
                vcf_i.identifier,
                vcf_i.refbase,
                vcf_i.altbase,
                vcf_i.qual,
                vcf_i.filters,
                vcf_i.info,
                field_format_line,
                samples_collect[0],
            )
        )
        outhandle.write(out_i + "\n")

    # Read the next line in the designated VCF file:
    line_i = my_vcf.readline().rstrip("\n")


# Close files:
my_vcf.close()
Tpileup.close()
outhandle.close()
if Npileup != None:
    Npileup.close()
