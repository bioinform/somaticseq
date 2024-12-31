#!/usr/bin/env python3

# Post-process GATK4's MuTect2 output. The main purpose is to split
# multi-allelic records into one variant record per line.

import argparse
import os
import re
import tempfile
import uuid

import somaticseq.genomic_file_parsers.genomic_file_handlers as genome
from somaticseq.vcf_modifier.bed_util import vcfsorter
from somaticseq.vcf_modifier.split_vcf import (
    split_complex_variants_into_snvs_and_indels,
)


def run() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-infile", "--input-vcf", type=str, help="Input VCF file", required=True
    )
    parser.add_argument(
        "-snv", "--snv-out", type=str, help="Output VCF file", required=True
    )
    parser.add_argument(
        "-indel", "--indel-out", type=str, help="Output VCF file", required=True
    )
    parser.add_argument(
        "-ref", "--genome-reference", type=str, help="genome reference", required=True
    )
    args = parser.parse_args()
    return args


def convert(infile, snv_out, indel_out, genome_reference):
    tempdir = tempfile.mkdtemp()
    tmp_snv_vcf = os.path.join(tempdir, uuid.uuid4().hex) + ".vcf"
    tmp_indel_vcf = os.path.join(tempdir, uuid.uuid4().hex) + ".vcf"
    info_to_split = "NLOD", "TLOD"
    info_to_keep = "STR", "ECNT"
    with (
        genome.open_textfile(infile) as vcf_in,
        open(tmp_snv_vcf, "w") as snvout,
        open(tmp_indel_vcf, "w") as indelout,
    ):
        line_i = vcf_in.readline().rstrip()
        while line_i.startswith("#"):
            snvout.write(line_i + "\n")
            indelout.write(line_i + "\n")
            line_i = vcf_in.readline().rstrip()

        while line_i:
            vcf_i = genome.VCFVariantRecord.from_vcf_line(line_i)

            # If "germlinerisk" is the only flag, then make it PASS since there
            # is no matched normal
            if vcf_i.filters == "germline_risk":
                vcf_i.filters = "PASS"

            if "," not in vcf_i.altbase:
                item = line_i.split("\t")
                if item[6] == "germline_risk":
                    item[6] = "PASS"

                new_line = "\t".join(item)
                if len(vcf_i.refbase) == 1 and len(vcf_i.altbase) == 1:
                    snvout.write(new_line + "\n")
                elif len(vcf_i.refbase) == 1 or len(vcf_i.altbase) == 1:
                    indelout.write(new_line + "\n")
                else:
                    snvs_and_indels = split_complex_variants_into_snvs_and_indels(vcf_i)
                    for snv_or_indel in snvs_and_indels:
                        if len(snv_or_indel.refbase) == len(snv_or_indel.altbase) == 1:
                            snvout.write(snv_or_indel.to_vcf_line() + "\n")
                        else:
                            indelout.write(snv_or_indel.to_vcf_line() + "\n")

            else:
                alt_bases = vcf_i.altbase.split(",")
                measures = []
                still_measures = []

                for measure_i in info_to_split:
                    try:
                        measures.append(vcf_i.get_info_value(measure_i).split(","))
                    except AttributeError:
                        measures.append(None)

                for measure_i in info_to_keep:
                    try:
                        still_measures.append(vcf_i.get_info_value(measure_i))
                    except AttributeError:
                        still_measures.append(None)

                for ith_base, altbase_i in enumerate(alt_bases):
                    split_infos = [
                        f"{info_variable}={info_value[ith_base]}"
                        for info_variable, info_value in zip(info_to_split, measures)
                        if info_value is not None
                    ]
                    still_infos = [
                        f"{info_variable}={info_value}"
                        for info_variable, info_value in zip(
                            info_to_keep, still_measures
                        )
                        if info_value is not False
                    ]
                    split_infos.extend(still_infos)
                    info_string = ";".join(split_infos)
                    GT0 = vcf_i.get_sample_value("GT", idx=0)
                    if GT0 != "0/0" and GT0 != "0/1":
                        sample_0 = re.sub(r"^[^:]+", "0/1", vcf_i.samples[0])
                    else:
                        sample_0 = vcf_i.samples[0]

                    new_line = "\t".join(
                        (
                            vcf_i.chromosome,
                            str(vcf_i.position),
                            vcf_i.identifier,
                            vcf_i.refbase,
                            altbase_i,
                            str(vcf_i.qual or "."),
                            vcf_i.filters,
                            info_string,
                            vcf_i.field,
                            sample_0,
                        )
                    )
                    if len(vcf_i.refbase) == 1 and len(altbase_i) == 1:
                        snvout.write(new_line + "\n")
                    elif len(vcf_i.refbase) == 1 or len(altbase_i) == 1:
                        indelout.write(new_line + "\n")
                    else:
                        complex_call = genome.VCFVariantRecord.from_vcf_line(new_line)
                        snvs_and_indels = split_complex_variants_into_snvs_and_indels(
                            complex_call
                        )
                        for snv_or_indel in snvs_and_indels:
                            if (
                                len(snv_or_indel.refbase)
                                == len(snv_or_indel.altbase)
                                == 1
                            ):
                                snvout.write(snv_or_indel.to_vcf_line() + "\n")
                            else:
                                indelout.write(snv_or_indel.to_vcf_line() + "\n")

            line_i = vcf_in.readline().rstrip()

    vcfsorter(genome_reference, tmp_snv_vcf, snv_out)
    vcfsorter(genome_reference, tmp_indel_vcf, indel_out)
    os.remove(tmp_snv_vcf)
    os.remove(tmp_indel_vcf)


if __name__ == "__main__":
    args = run()
    convert(args.input_vcf, args.snv_out, args.indel_out, args.genome_reference)
