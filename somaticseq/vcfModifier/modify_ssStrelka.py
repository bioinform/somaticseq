#!/usr/bin/env python3

# Split into SNV and INDEL VCF files

import argparse
import re

import somaticseq.genomicFileHandler.genomic_file_handlers as genome


def run():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Variant Call Type, i.e., snp or indel
    parser.add_argument(
        "-infile", "--input-vcf", type=str, help="Input VCF file", required=True
    )
    parser.add_argument(
        "-snv", "--output-snv", type=str, help="Output SNV VCF file", required=True
    )
    parser.add_argument(
        "-indel",
        "--output-indel",
        type=str,
        help="Output INDEL VCF file",
        required=True,
    )

    # Parse the arguments:
    args = parser.parse_args()

    infile = args.input_vcf
    snv_out = args.output_snv
    indel_out = args.output_indel

    return infile, snv_out, indel_out


def convert(infile, snv_out, indel_out):

    info_to_split = "REFREP", "IDREP", "RU"
    info_to_keep = ("MQ",)

    with genome.open_textfile(infile) as vcf_in, open(snv_out, "w") as snv_out, open(
        indel_out, "w"
    ) as indel_out:

        line_i = vcf_in.readline().rstrip()

        while line_i.startswith("##"):

            snv_out.write(line_i + "\n")
            indel_out.write(line_i + "\n")
            line_i = vcf_in.readline().rstrip()

        # This is the #CHROM line:
        headers = line_i.split("\t")
        snv_out.write(line_i + "\n")
        indel_out.write(line_i + "\n")

        line_i = vcf_in.readline().rstrip()
        while line_i:

            items = line_i.split("\t")

            vcf_i = genome.VcfLine(line_i)

            if "," not in vcf_i.altbase:

                if len(vcf_i.refbase) == 1 and len(vcf_i.altbase) == 1:
                    snv_out.write(line_i + "\n")
                else:
                    indel_out.write(line_i + "\n")

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
                        "{}={}".format(info_variable, info_value[ith_base])
                        for info_variable, info_value in zip(info_to_split, measures)
                        if info_value != None
                    ]

                    still_infos = [
                        "{}={}".format(info_variable, info_value)
                        for info_variable, info_value in zip(
                            info_to_keep, still_measures
                        )
                        if info_value != False
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
                            vcf_i.qual,
                            vcf_i.filters,
                            info_string,
                            vcf_i.field,
                            sample_0,
                        )
                    )

                    if len(vcf_i.refbase) == 1 and len(altbase_i) == 1:
                        snv_out.write(new_line + "\n")
                    else:
                        indel_out.write(new_line + "\n")

            line_i = vcf_in.readline().rstrip()


if __name__ == "__main__":
    infile, snv_out, indel_out = run()
    convert(infile, snv_out, indel_out)
