#!/usr/bin/env python3

import argparse
import gzip
import re
import sys

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "-fai", "--fai-file", type=str, help=".fa.fai file", required=True, default=None
)
parser.add_argument(
    "-beds",
    "--bed-files",
    type=str,
    help="BED files",
    nargs="*",
    required=True,
    default=None,
)
parser.add_argument(
    "-out ",
    "--bed-out",
    type=str,
    help="BED file out",
    required=False,
    default=sys.stdout,
)

args = parser.parse_args()

fai_file = args.fai_file
bed_files = args.bed_files
bed_out = args.bed_out


def fai2bed(file_name):

    with open(file_name) as gfile:
        line_i = gfile.readline().rstrip("\n")

        callableLociBoundries = {}
        callableLociCounters = {}
        orderedContig = []
        while line_i:

            contig_match = re.match(r"([^\t]+)\t", line_i)

            if contig_match:

                contig_i = contig_match.groups()[0].split(" ")[
                    0
                ]  # some .fai files have space after the contig for descriptions.
                orderedContig.append(contig_i)

                contig_size = int(line_i.split("\t")[1])

                callableLociBoundries[contig_i] = [0, contig_size]
                callableLociCounters[contig_i] = [0]

            else:
                raise Exception(".fai file format not as expected.")

            line_i = gfile.readline().rstrip("\n")

    return callableLociBoundries, callableLociCounters, orderedContig


def bed2regions(bed_file):

    regions = {}

    with open(bed_file) as gfile:

        line_i = gfile.readline().rstrip("\n")

        while line_i:
            item = line_i.split("\t")

            chrom = item[0]
            startPos = int(item[1])
            endPos = int(item[2])

            try:
                regions[chrom].append((startPos, endPos))
            except KeyError:
                regions[chrom] = []
                regions[chrom].append((startPos, endPos))

            line_i = gfile.readline().rstrip("\n")

    return regions


def collapseIdenticalBoundries(boundries, counters):

    assert len(boundries) == len(counters) + 1

    outBoundries = []
    outCounters = []

    outBoundries.append(boundries[0])

    i = 0
    while i <= len(boundries) - 2:
        j = i + 1

        if boundries[i] != boundries[j]:
            outBoundries.append(boundries[j])
            outCounters.append(counters[j - 1])

        i += 1

    return outBoundries, outCounters


def countIntersectedRegions(original_boundry, original_counter, additional_regions):

    secondary_boundry = []
    secondary_counter = []
    for region_i in additional_regions:
        secondary_boundry.append(region_i[0])
        secondary_boundry.append(region_i[1])
        secondary_counter.append(1)
        secondary_counter.append(0)

    newBoundry = []
    newCounter = []

    secondaryIterator = iter(secondary_boundry)
    boundry_j = next(secondaryIterator)
    j = 0
    counter_j = secondary_counter[j]
    secondaryMore = True

    for i, boundry_i in enumerate(original_boundry):

        if secondaryMore:
            while boundry_j < boundry_i:

                newBoundry.append(boundry_j)

                if j % 2 == 0:
                    newCounter.append(counter_i + 1)
                elif j % 2 == 1:
                    newCounter.append(counter_i)

                try:
                    boundry_j = next(secondaryIterator)
                    j += 1
                except StopIteration:
                    secondaryMore = False
                    break

        # Move onto the next original boundry
        newBoundry.append(boundry_i)

        try:
            counter_i = original_counter[i]
        except IndexError:
            counter_i = None

        if not secondaryMore and counter_i != None:
            newCounter.append(counter_i)
        elif j % 2 == 0:
            newCounter.append(counter_i)
        elif j % 2 == 1:
            if counter_i != None:
                newCounter.append(counter_i + 1)

    consolidatedBoundries, consolidatedCounters = collapseIdenticalBoundries(
        newBoundry, newCounter
    )

    return consolidatedBoundries, consolidatedCounters


# Start routine:
contigBoundries, contigCounters, orderedContigs = fai2bed(fai_file)

# Look at BED files
for bed_file_i in bed_files:

    bedRegions = bed2regions(bed_file_i)

    for chrom in bedRegions:
        contigBoundries[chrom], contigCounters[chrom] = countIntersectedRegions(
            contigBoundries[chrom], contigCounters[chrom], bedRegions[chrom]
        )


for contig_i in orderedContigs:

    if contigCounters[contig_i] != [0]:

        for i, count_i in enumerate(contigCounters[contig_i]):

            out_string = "{}\t{}\t{}\t{}".format(
                contig_i,
                contigBoundries[contig_i][i],
                contigBoundries[contig_i][i + 1],
                count_i,
            )

            bed_out.write(out_string + "\n")
