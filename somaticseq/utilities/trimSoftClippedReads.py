#!/usr/bin/env python3

import argparse
import re

import pysam

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "-bamin",
    "--bam-file-in",
    type=str,
    help="Input BAM file",
    required=True,
    default=None,
)
parser.add_argument(
    "-bamout",
    "--bam-file-out",
    type=str,
    help="Output BAM file",
    required=True,
    default=None,
)

args = parser.parse_args()
bam_file = args.bam_file_in
bam_out = args.bam_file_out

with pysam.AlignmentFile(bam_file) as bam, pysam.AlignmentFile(
    bam_out, "wb", template=bam
) as bamout:

    reads = bam.fetch()

    total_trimmed_bases = 0

    for read_i in reads:

        if read_i.cigarstring and "S" in read_i.cigarstring:

            front_clipped = re.search(r"^([0-9]+)S", read_i.cigarstring)
            back_clipped = re.search(r"([0-9]+)S$", read_i.cigarstring)

            if front_clipped and back_clipped:

                front_num = int(front_clipped.groups()[0])
                back_num = int(back_clipped.groups()[0])
                read_i.cigarstring = re.sub(
                    r"^[0-9]+S|[0-9]+S$", "", read_i.cigarstring
                )

                qual_i = read_i.qual[front_num::][:-back_num]

                read_i.seq = read_i.seq[front_num::][:-back_num]
                read_i.qual = qual_i

                if read_i.has_tag("BI"):
                    read_i.set_tag(
                        tag="BI",
                        value=read_i.get_tag("BI")[front_num::][:-back_num],
                        value_type="Z",
                        replace=True,
                    )

                if read_i.has_tag("BD"):
                    read_i.set_tag(
                        tag="BD",
                        value=read_i.get_tag("BD")[front_num::][:-back_num],
                        value_type="Z",
                        replace=True,
                    )

                total_trimmed_bases += front_num
                total_trimmed_bases += back_num

            elif front_clipped:

                num_bases = int(front_clipped.groups()[0])
                read_i.cigarstring = re.sub(r"^([0-9]+)S", "", read_i.cigarstring)

                qual_i = read_i.qual[num_bases::]

                read_i.seq = read_i.seq[num_bases::]
                read_i.qual = qual_i

                if read_i.has_tag("BI"):
                    read_i.set_tag(
                        tag="BI",
                        value=read_i.get_tag("BI")[num_bases::],
                        value_type="Z",
                        replace=True,
                    )

                if read_i.has_tag("BD"):
                    read_i.set_tag(
                        tag="BD",
                        value=read_i.get_tag("BD")[num_bases::],
                        value_type="Z",
                        replace=True,
                    )

                total_trimmed_bases += num_bases

            elif back_clipped:

                num_bases = int(back_clipped.groups()[0])
                read_i.cigarstring = re.sub("[0-9]+S$", "", read_i.cigarstring)

                qual_i = read_i.qual[:-num_bases]

                read_i.seq = read_i.seq[:-num_bases]
                read_i.qual = qual_i

                if read_i.has_tag("BI"):
                    read_i.set_tag(
                        tag="BI",
                        value=read_i.get_tag("BI")[:-num_bases],
                        value_type="Z",
                        replace=True,
                    )

                if read_i.has_tag("BD"):
                    read_i.set_tag(
                        tag="BD",
                        value=read_i.get_tag("BD")[:-num_bases],
                        value_type="Z",
                        replace=True,
                    )

                total_trimmed_bases += num_bases

        # Mate CIGAR
        if read_i.has_tag("MC"):
            mate_cigar = read_i.get_tag("MC")
            if "S" in mate_cigar:
                new_cigar = re.sub(r"^[0-9]+S|[0-9]+S$", "", mate_cigar)
                read_i.set_tag(tag="MC", value=new_cigar, value_type="Z", replace=True)

        bamout.write(read_i)

print("A total of {} bases trimmed.".format(total_trimmed_bases))
