#!/usr/bin/env python3

import argparse
import math
import os
import re
import sys


def run() -> tuple[str, str, int]:
    # argparse Stuff
    parser = argparse.ArgumentParser(
        description=(
            "Given an input bed file, it will output a number of bed files, "
            "each with the same number of total base pairs. "
            "This routine is used to parallelize SomaticSeq tasks. "
            "One contiguous region may be split into multiple regions in "
            "multiple output files. The limitation is that some regions of the "
            "genome have much higher coverage than others. "
            "This is the reason some regions run much slower than others. "
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-infile",
        "--input-file",
        type=str,
        help="Input merged BED file",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-num",
        "--num-of-files",
        type=int,
        help="Number of bed files to split into",
        required=False,
        default=1,
    )
    parser.add_argument(
        "-outfiles",
        "--output-files",
        type=str,
        help="Output BED file",
        required=False,
        default=sys.stdout,
    )
    args = parser.parse_args()
    return args.input_file, args.output_files, args.num_of_files


def fai2bed(fai: str, bedout: str) -> str:
    """Create a whole genome bed file based on .fai file."""
    with open(fai) as fai_h, open(bedout, "w") as bed_h:
        fai_i = fai_h.readline().rstrip()
        while fai_i:
            fai_item = fai_i.split("\t")
            bed_h.write(f"{fai_item[0]}\t0\t{fai_item[1]}\n")
            fai_i = fai_h.readline().rstrip()
    return bedout


def split(infile: str, outfiles: str, num: int) -> list[str]:
    """Split a bed file into n bed files of equal-sized regions

    Args:
        infile: input bed file
        outfiles: output bed files to which "n." will be appended to its basename,
            e.g., if outfiles="/PATH/TO/important_regions.bed", the output bed files
            written will be /PATH/TO/1.important_regions.bed,
            /PATH/TO/2.important_regions.bed, etc.
        num: number of output bed files

    Returns:
        List of output bed files written
    """
    outfiles_written = []
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
    base_pairs_per_bed = math.ceil(total_region_size / num)

    # Go through every original region and split
    current_size = 0
    current_region = []
    ith_split = 1
    for region_i in original_regions:
        chr_i = region_i[0]
        start_i = region_i[1]
        end_i = region_i[2]

        # If the "now size" is still less than size/file requirement
        if current_size + (end_i - start_i) <= base_pairs_per_bed:
            # Need to collect more to fulfill the size/file requirement, so
            # append to current_region list
            current_region.append(f"{chr_i}\t{start_i}\t{end_i}\n")
            current_size = current_size + (end_i - start_i)

        # If the "now size" exceeds the size/file requirement, need to start
        # splitting:
        elif current_size + (end_i - start_i) > base_pairs_per_bed:
            # Split a big region into a smaller regino, such that the size of
            # "current_region" is equal to the size/file requirement:
            breakpoint_i = base_pairs_per_bed + start_i - current_size

            # Write these regions out, , reset "current_region," then add 1 to
            # ith_split afterward to keep track:
            outfiles_written.append(
                os.path.join(out_directory, f"{ith_split}.{out_basename}")
            )
            with open(
                os.path.join(out_directory, f"{ith_split}.{out_basename}"), "w"
            ) as ith_out:
                for line_i in current_region:
                    ith_out.write(line_i)

                # Make sure it doesn't write a 0-bp region:
                if breakpoint_i > start_i:
                    ith_out.write(f"{chr_i}\t{start_i}\t{breakpoint_i}\n")
            ith_split += 1
            current_region = []

            # The remaining, is the end position of the original region and the
            # previous breakpoint:
            remaining_length = end_i - breakpoint_i

            # If the remnant of the previous region is less than the size/file
            # requirement, simply make it "current_region" and then move on:
            if remaining_length <= base_pairs_per_bed:
                current_region.append(f"{chr_i}\t{breakpoint_i}\t{end_i}\n")
                current_size = remaining_length

            # If the renmant of the previuos region exceed the size/file
            # requirement, it needs to be properly split until it's small
            # enough:
            elif remaining_length > base_pairs_per_bed:
                # Each sub-region, if large enough, will have its own file output:
                while (end_i - breakpoint_i) > base_pairs_per_bed:
                    end_j = breakpoint_i + base_pairs_per_bed
                    outfiles_written.append(
                        os.path.join(out_directory, f"{ith_split}.{out_basename}")
                    )
                    with open(
                        os.path.join(out_directory, f"{ith_split}.{out_basename}"),
                        "w",
                    ) as ith_out:
                        if end_j > breakpoint_i:
                            ith_out.write(f"{chr_i}\t{breakpoint_i}\t{end_j}\n")
                    ith_split += 1
                    breakpoint_i = end_j

                # After every sub-region has its own bed file, the remnant is
                # added to "current_region" to deal with the next line of the
                # "original_regions"
                current_region.append(f"{chr_i}\t{breakpoint_i}\t{end_i}\n")
                current_size = end_i - breakpoint_i

    # The final region to write out:
    ith_basename = os.path.join(out_directory, f"{ith_split}.{out_basename}")
    outfiles_written.append(ith_basename)
    with open(ith_basename, "w") as ith_out:
        for line_i in current_region:
            ith_out.write(line_i)

    return outfiles_written


def main() -> None:
    infile, outfiles, num = run()
    split(infile, outfiles, num)


if __name__ == "__main__":
    main()
