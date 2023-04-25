#!/usr/bin/env python3

# Add GT to Strelka's samples to make compatible with GATK CombineVariants, so don't care about the content. Just 0/1 for everyone.

import argparse

import somaticseq.genomicFileHandler.genomic_file_handlers as genome


def run():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Variant Call Type, i.e., snp or indel
    parser.add_argument(
        "-infile", "--input-vcf", type=str, help="Input VCF file", required=True
    )
    parser.add_argument(
        "-outfile", "--output-vcf", type=str, help="Output VCF file", required=True
    )

    # Parse the arguments:
    args = parser.parse_args()
    infile = args.input_vcf
    outfile = args.output_vcf

    return infile, outfile


def convert(infile, outfile):
    with genome.open_textfile(infile) as vcf_in, open(outfile, "w") as vcf_out:

        line_i = vcf_in.readline().rstrip()

        while line_i.startswith("##"):

            vcf_out.write(line_i + "\n")
            line_i = vcf_in.readline().rstrip()

        # This is the #CHROM line:
        headers = line_i.split("\t")
        num_columns = len(headers)
        vcf_out.write(line_i + "\n")

        line_i = vcf_in.readline().rstrip()
        while line_i:

            items = line_i.split("\t")

            items[8] = "GT:" + items[8]

            for i in range(9, num_columns):
                items[i] = "0/1:" + items[i]

            line_out = "\t".join(items)
            vcf_out.write(line_out + "\n")

            line_i = vcf_in.readline().rstrip()


if __name__ == "__main__":
    infile, outfile = run()
    convert(infile, outfile)
