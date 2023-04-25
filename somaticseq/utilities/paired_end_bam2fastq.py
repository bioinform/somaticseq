#!/usr/bin/env python3

import argparse
import gzip

import pysam

NT_PAIRS = {
    "G": "C",
    "T": "A",
    "C": "G",
    "A": "T",
    "N": "N",
    "g": "c",
    "t": "a",
    "c": "g",
    "a": "t",
    "n": "n",
}


def reverse_complement(seq):
    seq_j = "".join([NT_PAIRS[base_i] for base_i in seq[::-1]])
    return seq_j


def text_open_write(filename):
    if str(filename).endswith(".gz"):
        return gzip.open(filename, "wt")
    else:
        return open(filename, "w")


def bam2fq(bam_file, fastq1, fastq2):

    with pysam.AlignmentFile(bam_file) as bam, text_open_write(
        fastq1
    ) as fq1, text_open_write(fastq2) as fq2:

        reads1 = {}
        reads2 = {}

        reads = bam.fetch()
        for read_i in reads:

            if not read_i.is_secondary:

                seq_i = (
                    reverse_complement(read_i.query_sequence)
                    if read_i.is_reverse
                    else read_i.query_sequence
                )
                qual_i = read_i.qual[::-1] if read_i.is_reverse else read_i.qual

                if read_i.is_read1:

                    if read_i.query_name in reads2:

                        fq1.write("@{}/1\n".format(read_i.query_name))
                        fq1.write(seq_i + "\n")
                        fq1.write("+\n")
                        fq1.write(qual_i + "\n")

                        read_2 = reads2.pop(read_i.query_name)

                        fq2.write("@{}/2\n".format(read_2["qname"]))
                        fq2.write(read_2["seq"] + "\n")
                        fq2.write("+\n")
                        fq2.write(read_2["bq"] + "\n")

                    else:
                        reads1[read_i.query_name] = {}
                        reads1[read_i.query_name]["qname"] = read_i.query_name
                        reads1[read_i.query_name]["seq"] = seq_i
                        reads1[read_i.query_name]["bq"] = qual_i

                elif read_i.is_read2:

                    if read_i.query_name in reads1:

                        read_1 = reads1.pop(read_i.query_name)

                        fq1.write("@{}/1\n".format(read_1["qname"]))
                        fq1.write(read_1["seq"] + "\n")
                        fq1.write("+\n")
                        fq1.write(read_1["bq"] + "\n")

                        fq2.write("@{}/2\n".format(read_i.query_name))
                        fq2.write(seq_i + "\n")
                        fq2.write("+\n")
                        fq2.write(qual_i + "\n")

                    else:
                        reads2[read_i.query_name] = {}
                        reads2[read_i.query_name]["qname"] = read_i.query_name
                        reads2[read_i.query_name]["seq"] = seq_i
                        reads2[read_i.query_name]["bq"] = qual_i

    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Convert paired-end BAM to FASTQ1 and 2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-bam", "--bam", type=str, help="bam file in")
    parser.add_argument("-fq1", "--fastq1", type=str, help="fastq1 out")
    parser.add_argument("-fq2", "--fastq2", type=str, help="fastq2 out")

    args = parser.parse_args()

    bam2fq(args.bam, args.fastq1, args.fastq2)
