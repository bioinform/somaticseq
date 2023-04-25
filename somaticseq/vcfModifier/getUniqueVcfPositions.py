#!/usr/bin/env python3

# A simple and quick way to replace GATK3 CombineVariants

import argparse
import gzip
import re


def open_textfile(file_name):

    # See if the input file is a .gz file:
    if file_name.lower().endswith(".gz"):
        return gzip.open(file_name, "rt")

    else:
        return open(file_name)


def run():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-vcfs",
        "--input-vcfs",
        nargs="*",
        type=str,
        help="Input VCF file",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-out", "--output-vcf", type=str, help="Output VCF file", required=True
    )

    args = parser.parse_args()

    infiles = args.input_vcfs
    outfile = args.output_vcf

    return infiles, outfile


def combine(infiles, outfile):

    variant_positions = set()

    for file_i in infiles:

        with open_textfile(file_i) as vcf:

            line_i = vcf.readline().rstrip()

            while line_i.startswith("#"):
                line_i = vcf.readline().rstrip()

            while line_i:

                item = line_i.split("\t")

                chromosome = item[0]
                position = int(item[1])
                refbase = item[3]
                altbases = re.split(r"[,/]", item[4])

                for altbase_i in altbases:
                    variant_positions.add((chromosome, position, refbase, altbase_i))

                line_i = vcf.readline().rstrip()

    with open(outfile, "w") as vcf_out:
        vcf_out.write("##fileformat=VCFv4.1\n")
        vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for variant_position_i in sorted(variant_positions):

            vcf_out.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    variant_position_i[0],
                    variant_position_i[1],
                    ".",
                    variant_position_i[2],
                    variant_position_i[3],
                    ".",
                    "PASS",
                    ".",
                )
            )


if __name__ == "__main__":
    infiles, outfile = run()
    combine(infiles, outfile)
