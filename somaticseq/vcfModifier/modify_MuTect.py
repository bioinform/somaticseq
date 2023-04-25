#!/usr/bin/env python3

# For single-sample mode, remove the "none" column from the VCF file to be consistent with VCF files from other single-sample tools.
# Keep Broad's convention for AD, i.e., allelic depths for the ref and alt alleles in the order listed
# 4/12/2015

import argparse
import sys

import somaticseq.genomicFileHandler.genomic_file_handlers as genome


def run():

    # argparse Stuff
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Variant Call Type, i.e., snp or indel
    parser.add_argument(
        "-infile",
        "--input-vcf",
        type=str,
        help="Input VCF file",
        required=True,
    )
    parser.add_argument(
        "-outfile", "--output-vcf", type=str, help="Output VCF file", default=sys.stdout
    )
    parser.add_argument(
        "-tbam",
        "--tumor-bam",
        type=str,
        help="A tumor bam file for sample name identification.",
    )
    parser.add_argument(
        "-nbam",
        "--normal-bam",
        type=str,
        help="A normal bam file for sample name identification.",
    )

    # Parse the arguments:
    args = parser.parse_args()

    in_vcf = args.input_vcf
    out_vcf = args.output_vcf
    tbam = args.tumor_bam
    nbam = args.normal_bam

    return in_vcf, out_vcf, tbam, nbam


def convert(infile, outfile, tbam, nbam=None):

    paired_mode = True if nbam else False

    # Get tumor and normal sample names from the bam files:
    nbam_header = genome.pysam_header(nbam) if nbam else None
    tbam_header = genome.pysam_header(tbam)

    # When MuTect is run in a "single sample mode," the "normal" will be named "none."
    n_samplename = nbam_header.SM() if nbam else ["none"]
    t_samplename = tbam_header.SM()

    if not (len(n_samplename) == 1 and len(t_samplename) == 1):
        sys.stderr.write("There are multiple Sample Names present in the BAM file!")

    n_samplename = n_samplename[0]
    t_samplename = t_samplename[0]

    assert t_samplename or n_samplename

    if t_samplename and n_samplename:
        paired_mode = True
    else:
        paired_mode = False

    (
        idx_chrom,
        idx_pos,
        idx_id,
        idx_ref,
        idx_alt,
        idx_qual,
        idx_filter,
        idx_info,
        idx_format,
    ) = (0, 1, 2, 3, 4, 5, 6, 7, 8)
    idx_SM1, idx_SM2 = 9, 10

    with genome.open_textfile(infile) as vcf, open(outfile, "w") as vcfout:

        line_i = vcf.readline().rstrip()

        while line_i.startswith("#"):

            if line_i.startswith("##"):
                vcfout.write(line_i + "\n")

            elif line_i.startswith("#CHROM"):
                header_items = line_i.rstrip().split("\t")

                idxN = header_items.index(n_samplename)
                idxT = header_items.index(t_samplename)

                if paired_mode:
                    header_items[idx_SM1] = "NORMAL"
                    header_items[idx_SM2] = "TUMOR"

                else:

                    # Keep up to the first sample column, then make sure it's labeled the TUMOR sample name
                    header_items = header_items[: idx_SM1 + 1]
                    header_items[idx_SM1] = "TUMOR"

                replaced_header = "\t".join(header_items)
                vcfout.write(replaced_header + "\n")

            line_i = vcf.readline().rstrip()

        while line_i:

            items_i = line_i.split("\t")

            if paired_mode:
                items_i[idx_SM1], items_i[idx_SM2] = items_i[idxN], items_i[idxT]

            else:
                items_i = items_i[:idx_SM1] + [items_i[idxT]]

            # Print the new stuff:
            new_line = "\t".join(items_i)

            # Have to get rid of "N" in REF, because after snpSift annotation, it changes the ALT and vcf-validator will complain.
            if not ("N" in items_i[idx_ref]):
                vcfout.write(new_line + "\n")

            line_i = vcf.readline().rstrip()


if __name__ == "__main__":
    infile, outfile, tbam, nbam = run()
    convert(infile, outfile, tbam, nbam)
