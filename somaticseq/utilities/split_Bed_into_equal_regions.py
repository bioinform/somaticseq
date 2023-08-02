#!/usr/bin/env python3

import argparse
import math
import os
import re
import sys


def run():

    # argparse Stuff
    parser = argparse.ArgumentParser(
        description="""Given an input bed file, this program will output a number of bed files, each will have same number of total base pairs. 
        This routine is used to parallelize SomaticSeq tasks. 
        One limitation, however, is that some regions of the genome have much higher coverage than others. 
        This is the reason some regions run much slower than others.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Variant Call Type, i.e., snp or indel
    parser.add_argument(
        "-infile",
        "--input-file",
        type=str,
        help="Input merged BED file",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-num", "--num-of-files", type=int, help="1", required=False, default=1
    )
    parser.add_argument(
        "-outfiles",
        "--output-files",
        type=str,
        help="Output BED file",
        required=False,
        default=sys.stdout,
    )

    # Parse the arguments:
    args = parser.parse_args()

    infile = args.input_file
    outfiles = args.output_files
    num = args.num_of_files

    return infile, outfiles, num


def fai2bed(fai, bedout):

    with open(fai) as fai, open(bedout, "w") as bed:

        fai_i = fai.readline().rstrip()

        while fai_i:
            fai_item = fai_i.split("\t")
            bed.write("{}\t{}\t{}\n".format(fai_item[0], "0", fai_item[1]))
            fai_i = fai.readline().rstrip()

    return bedout


def split(infile, outfiles, num):

    outfilesWritten = []

    out_basename = os.path.basename(outfiles)
    out_directory = os.path.dirname(outfiles)
    os.makedirs(out_directory, exist_ok=True)

    if not out_directory:
        out_directory = os.curdir

    with open(infile) as bedin:

        line_i = bedin.readline().rstrip()

        while re.match(r"track|browser|#", line_i):
            line_i = bedin.readline().rstrip()

        total_region_size = 0
        original_regions = []
        while line_i:

            items = line_i.split("\t")

            chr_i = items[0]
            start_i = int(items[1])
            end_i = int(items[2])

            total_region_size = total_region_size + (end_i - start_i)
            original_regions.append((chr_i, start_i, end_i))

            line_i = bedin.readline().rstrip()

    # For each bed file, this is the total base pairs in that file
    size_per_file = math.ceil(total_region_size / num)

    # Go through every original region and split
    current_size = 0
    current_region = []
    ith_split = 1
    for region_i in original_regions:

        chr_i = region_i[0]
        start_i = region_i[1]
        end_i = region_i[2]

        # If the "now size" is still less than size/file requirement
        if current_size + (end_i - start_i) <= size_per_file:

            # Need to collect more to fulfill the size/file requirement, so append to current_region list
            current_region.append("{}\t{}\t{}\n".format(chr_i, start_i, end_i))
            current_size = current_size + (end_i - start_i)

        # If the "now size" exceeds the size/file requirement, need to start splitting:
        elif current_size + (end_i - start_i) > size_per_file:

            # Split a big region into a smaller regino, such that the size of "current_region" is equal to the size/file requirement:
            breakpoint_i = size_per_file + start_i - current_size

            # Write these regions out, , reset "current_region," then add 1 to ith_split afterward to keep track:
            outfilesWritten.append(
                "{}{}{}.{}".format(out_directory, os.sep, ith_split, out_basename)
            )
            with open(
                "{}{}{}.{}".format(out_directory, os.sep, ith_split, out_basename), "w"
            ) as ith_out:
                for line_i in current_region:
                    ith_out.write(line_i)

                # Make sure it doesn't write a 0-bp region:
                if breakpoint_i > start_i:
                    ith_out.write("{}\t{}\t{}\n".format(chr_i, start_i, breakpoint_i))
            ith_split += 1
            current_region = []

            # The remaining, is the end position of the original region and the previous breakpoint:
            remaining_length = end_i - breakpoint_i
            remaining_region = (chr_i, breakpoint_i, end_i)

            # If the remnant of the previous region is less than the size/file requirement, simply make it "current_region" and then move on:
            if remaining_length <= size_per_file:
                current_region.append("{}\t{}\t{}\n".format(chr_i, breakpoint_i, end_i))
                current_size = remaining_length

            # If the renmant of the previuos region exceed the size/file requirement, it needs to be properly split until it's small enough:
            elif remaining_length > size_per_file:

                # Each sub-region, if large enough, will have its own file output:
                while (end_i - breakpoint_i) > size_per_file:

                    end_j = breakpoint_i + size_per_file

                    outfilesWritten.append(
                        "{}{}{}.{}".format(
                            out_directory, os.sep, ith_split, out_basename
                        )
                    )
                    with open(
                        "{}{}{}.{}".format(
                            out_directory, os.sep, ith_split, out_basename
                        ),
                        "w",
                    ) as ith_out:

                        if end_j > breakpoint_i:
                            ith_out.write(
                                "{}\t{}\t{}\n".format(chr_i, breakpoint_i, end_j)
                            )
                    ith_split += 1

                    breakpoint_i = end_j

                # After every sub-region has its own bed file, the remnant is added to "current_region" to deal with the next line of the "original_regions"
                current_region.append("{}\t{}\t{}\n".format(chr_i, breakpoint_i, end_i))
                current_size = end_i - breakpoint_i

    # The final region to write out:
    ithOutName = "{}{}{}.{}".format(out_directory, os.sep, ith_split, out_basename)
    outfilesWritten.append(ithOutName)
    with open(ithOutName, "w") as ith_out:
        for line_i in current_region:
            ith_out.write(line_i)

    return outfilesWritten


def split_vcf_file(vcf_file, work_dir=os.curdir, num=1):

    num_lines = 0
    with open_textfile(vcf_file) as vcf:
        line_i = vcf.readline()
        header = []
        while line_i.startswith("#"):
            header.append(line_i)
            line_i = vcf.readline()
        while line_i:
            num_lines += 1
            line_i = vcf.readline()

    lines_per_file = math.ceil(float(num_lines) / num)

    with open_textfile(vcf_file) as vcf:

        outnames = [
            os.curdir
            + os.sep
            + str(i)
            + "_"
            + re.sub(r".vcf(.gz)?", "", os.path.basename(vcf_file))
            + ".vcf"
            for i in range(num)
        ]
        outhandles = [open(i, "w") for i in outnames]
        [write_header(header, i) for i in outhandles]

        line_i = vcf.readline()

        while line_i.startswith("#"):
            line_i = vcf.readline()

        while line_i:

            i = 0
            n = 0
            while line_i:

                outhandles[n].write(line_i)
                i += 1

                if i == lines_per_file:
                    i = 0
                    n += 1

                line_i = vcf.readline()

        [i.close() for i in outhandles]

    return outnames


if __name__ == "__main__":
    infile, outfiles, num = run()
    split(infile, outfiles, num)
