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
    parser.add_argument(
        "-tnscope",
        "--is-tnscope",
        action="store_true",
        help="Actually TNscope VCF",
        required=False,
        default=False,
    )
    args = parser.parse_args()
    return args


def convert(infile, snv_out, indel_out, is_tnscope, genome_reference):
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
        while line_i.startswith("##"):
            if line_i.startswith("##normal_sample="):
                normal_name = line_i.split("=")[1]
            if line_i.startswith("##tumor_sample="):
                tumor_name = line_i.split("=")[1]
            if line_i.startswith("##INFO=<ID=SOR,"):
                line_i = re.sub(r"Float", "String", line_i)
            snvout.write(line_i + "\n")
            indelout.write(line_i + "\n")
            line_i = vcf_in.readline().rstrip()

        # This line will be #CHROM:
        snvout.write(line_i + "\n")
        indelout.write(line_i + "\n")
        header = line_i.split("\t")

        if is_tnscope:
            # Doesn't matter which one is normal/tumor. These information are not used.
            normal_index, tumor_index = 1, 0
        else:
            normal_index = header.index(normal_name) - 9
            tumor_index = header.index(tumor_name) - 9

        # This will be the first variant line:
        line_i = vcf_in.readline().rstrip()
        while line_i:
            vcf_i = genome.VCFVariantRecord.from_vcf_line(line_i)
            if "," not in vcf_i.altbase:
                if len(vcf_i.refbase) == 1 and len(vcf_i.altbase) == 1:
                    snvout.write(line_i + "\n")
                elif len(vcf_i.refbase) == 1 or len(vcf_i.altbase) == 1:
                    indelout.write(line_i + "\n")
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

                    gt_n = vcf_i.get_sample_value("GT", idx=normal_index)
                    if gt_n != "0/0" and gt_n != "0/1":
                        sample_n = re.sub(r"^[^:]+", "0/1", vcf_i.samples[normal_index])
                    else:
                        sample_n = vcf_i.samples[normal_index]

                    gt_1 = vcf_i.get_sample_value("GT", idx=tumor_index)
                    if gt_1 != "0/0" and gt_n != "0/1":
                        sample_t = re.sub(r"^[^:]+", "0/1", vcf_i.samples[tumor_index])
                    else:
                        sample_t = vcf_i.samples[tumor_index]

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
                            sample_n,
                            sample_t,
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
    convert(
        args.input_vcf,
        args.snv_out,
        args.indel_out,
        args.is_tnscope,
        args.genome_reference,
    )
