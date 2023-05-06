#!/usr/bin/env python3

import argparse
import os
from multiprocessing import Pool

import pysam
import somaticseq.genomicFileHandler.genomic_file_handlers as genome


def bgzip_compress(infile, remove_infile=True):
    pysam.tabix_compress(infile, infile + ".gz", force=True)
    os.remove(infile)
    return infile + ".gz"


def vcf(infileList, outfile, bgzip=False):
    with open(outfile, "w") as vcfout:
        headerWritten = False
        for file_i in infileList:
            with genome.open_textfile(file_i) as vcfin:
                line_i = vcfin.readline()
                while line_i.startswith("#"):
                    if not headerWritten:
                        vcfout.write(line_i)

                    line_i = vcfin.readline()

                # Turn off header writing from now on:
                headerWritten = True
                while line_i:
                    vcfout.write(line_i)
                    line_i = vcfin.readline()

    if bgzip:
        bgzip_compress(outfile, True)
        actual_outfile = outfile + ".gz"
    else:
        actual_outfile = outfile

    return actual_outfile


def tsv(infileList, outfile, bgzip=False):
    with open(outfile, "w") as tsvout:
        headerWritten = False
        for file_i in infileList:
            with genome.open_textfile(file_i) as tsvin:
                # First line is a header
                line_i = tsvin.readline()
                if not headerWritten:
                    tsvout.write(line_i)

                # Turn off header writing from now on:
                headerWritten = True
                line_i = tsvin.readline()
                while line_i:
                    tsvout.write(line_i)
                    line_i = tsvin.readline()
    if bgzip:
        bgzip_compress(outfile, True)
        actual_outfile = outfile + ".gz"
    else:
        actual_outfile = outfile

    return actual_outfile


def bed(infileList, outfile, bgzip=False):
    with open(outfile, "w") as bedout:
        for file_i in infileList:
            with genome.open_textfile(file_i) as bedin:
                for line_i in bedin:
                    bedout.write(line_i)
    if bgzip:
        bgzip_compress(outfile, True)
        actual_outfile = outfile + ".gz"
    else:
        actual_outfile = outfile

    return actual_outfile


def spreader(infileList, outfiles, chunk=4, bgzip=False, threads=1):
    """
    Given an infile, it will spread its content into the outfiles "chunk" at a time, e.g,.
    If infile is a fastq file, and output is 3 fastq files, then the first 4 lines will go to the 1st output, the next 4 lines to go the 2nd output, the next 4 lines go to the 3rd output, and then the next 4 lines will go back to the 1st output, so on and so forth.
    """

    outs = [open(out_i, "w") for out_i in outfiles]
    for infile in infileList:
        with genome.open_textfile(infile) as text_in:
            line_i = text_in.readline()
            while line_i:
                for out_i in outs:
                    for i in range(chunk):
                        out_i.write(line_i)
                        line_i = text_in.readline()

    [out_i.close() for out_i in outs]

    if bgzip:
        pool = Pool(processes=threads)
        bash_async = pool.map_async(bgzip_compress, outfiles)
        actual_outfiles = bash_async.get()
        pool.close()

    else:
        actual_outfiles = outfiles

    return actual_outfiles


def run():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Variant Call Type, i.e., snp or indel
    parser.add_argument(
        "-infiles", "--input-files", type=str, nargs="*", help="Input files"
    )
    parser.add_argument("-outfile", "--output-file", type=str, help="Output file")
    parser.add_argument(
        "-outfiles",
        "--output-files",
        type=str,
        nargs="*",
        help="Output files for spreader",
    )
    parser.add_argument(
        "-chunk",
        "--chunk-size",
        type=int,
        help="In --spread mode, the number of lines to be written into the output file each time. By default chunk=4 by default for fastq files, i.e., every 4 lines make up one read record.",
        default=4,
    )
    parser.add_argument(
        "-nt",
        "--threads",
        type=int,
        help="only invoked in -spread -bgzip when bgzip compress of output files can be parallelized",
    )
    parser.add_argument(
        "-spread",
        "--spread",
        action="store_true",
        help="Spread content into multiple files.",
    )
    parser.add_argument(
        "-bgzip",
        "--bgzip-output",
        action="store_true",
        help="compress the output files",
    )

    # Parse the arguments:
    args = parser.parse_args()
    if args.spread:
        filetype = "spread"
    elif args.input_files[0].lower().endswith(".vcf") or args.input_files[
        0
    ].lower().endswith(".vcf.gz"):
        filetype = "vcf"
    elif args.input_files[0].lower().endswith(".tsv") or args.input_files[
        0
    ].lower().endswith(".tsv.gz"):
        filetype = "tsv"
    elif args.input_files[0].lower().endswith(".bed") or args.input_files[
        0
    ].lower().endswith(".bed.gz"):
        filetype = "bed"
    else:
        filetype = "unknown"
    return args, filetype


if __name__ == "__main__":
    args, ftype = run()
    if ftype == "spread":
        spreader(
            args.input_files,
            args.output_files,
            args.chunk_size,
            args.bgzip_output,
            args.threads,
        )
    elif ftype == "vcf":
        vcf(args.input_files, args.output_file, args.bgzip_output)
    elif ftype == "bed":
        bed(args.input_files, args.output_file, args.bgzip_output)
    elif ftype == "tsv":
        tsv(args.input_files, args.output_file, args.bgzip_output)
    elif ftype == "unknown":
        tsv(args.input_files, args.output_file, args.bgzip_output)
