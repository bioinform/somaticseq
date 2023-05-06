#!/usr/bin/env python3

import argparse
import gzip
import itertools
import math
import os
import re
import sys

import somaticseq.genomicFileHandler.genomic_file_handlers as genome

# argparse Stuff
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "-vcf",
    "--input-vcf",
    type=str,
    help="SomaticSeq VCF file",
    required=True,
    default=None,
)
parser.add_argument(
    "-combo",
    "--combo-code",
    type=str,
    help="E.g., MVJSDULK",
    required=True,
    default="MVJSDULK",
)

args = parser.parse_args()
vcf = args.input_vcf
combo = args.combo_code

tool_code = list(combo)

all_combos = {}
for i in range(1, len(tool_code) + 1):
    combo_gen = itertools.combinations(tool_code, i)
    for j in combo_gen:
        all_combos[j] = [0, 0]


with open(vcf) as vcf:

    line_i = vcf.readline().rstrip()

    while line_i.startswith("#"):
        line_i = vcf.readline().rstrip()

    print("#ToolCombo\tTruePositiveCalls\tAllCalls")

    while line_i:

        vcf_i = genome.VcfLine(line_i)
        combo_i = vcf_i.get_info_value(combo)
        tool_i = combo_i.split(",")
        tool_i = [int(i) for i in tool_i]

        current_call_set = set()
        for tool_code_j, tool_j in zip(tool_code, tool_i):
            if tool_j == 1:
                current_call_set.add(tool_code_j)

        for combo_j in all_combos:
            if set.intersection(set(combo_j), current_call_set):
                all_combos[combo_j][0] += 1

                if "TruePositive" in vcf_i.identifier:
                    all_combos[combo_j][1] += 1

        line_i = vcf.readline().rstrip()


for i in sorted(all_combos):
    print("".join(i) + "\t" + str(all_combos[i][1]) + "\t" + str(all_combos[i][0]))
