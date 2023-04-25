#!/usr/bin/env python3

import argparse

import somaticseq.genomicFileHandler.genomic_file_handlers as genome


def run():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Variant Call Type, i.e., snp or indel
    parser.add_argument(
        "-infile", "--input-file", type=str, help="Input VCF file", required=True
    )
    parser.add_argument(
        "-outfile", "--output-file", type=str, help="Output VCF file", required=True
    )

    # Parse the arguments:
    args = parser.parse_args()
    infile = args.input_file
    outfile = args.output_file

    return infile, outfile


def copy(infile, outfile):

    with genome.open_textfile(infile) as filein, open(outfile, "w") as fileout:
        line_i = filein.readline()
        while line_i:
            fileout.write(line_i)
            line_i = filein.readline()


if __name__ == "__main__":
    infile, outfile = run()
    copy(infile, outfile)
