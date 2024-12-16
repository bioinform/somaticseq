#!/usr/bin/env python3

# Post-process GATK4's MuTect2 output. The main purpose is to split multi-allelic records into one variant record per line.
# The objective is to seperate SNV and INDEL into two files.
# The primary objective is to extract variant candidate list, i.e., genomic coordinate, reference, and alternate alleles.
# It does NOT try to resolve different genotypes if there are multiple on a single VCF file. So keep in mind of this limitation.

import argparse
from copy import copy

import somaticseq.genomic_file_parsers.genomic_file_handlers as genome
from somaticseq.vcf_modifier.complex2indel import (
    resolve_complex_variants_into_snvs_and_indels,
)


def run():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Split a VCF file into SNVs and INDELs.",
    )
    # Variant Call Type, i.e., snp or indel
    parser.add_argument(
        "-infile", "--input-vcf", type=str, help="Input VCF file", required=True
    )
    parser.add_argument(
        "-indel", "--indel-out", type=str, help="Output VCF file", required=True
    )
    parser.add_argument(
        "-snv", "--snv-out", type=str, help="Output VCF file", required=True
    )
    # Parse the arguments:
    args = parser.parse_args()
    return args.input_vcf, args.snv_out, args.indel_out


def split_complex_variants_into_snvs_and_indels(
    vcf_record: genome.VCFVariantRecord,
) -> list[genome.VCFVariantRecord]:
    """
    Split a complex variant into snvs and indels
    """
    assert vcf_record.refbase
    assert vcf_record.altbase
    decomplexed_variants = resolve_complex_variants_into_snvs_and_indels(
        refbases=vcf_record.refbase, altbases=vcf_record.altbase
    )
    assert decomplexed_variants
    snv_and_indel_records = []
    for decomplex_dict in decomplexed_variants:
        offset = decomplex_dict["OFFSET"]
        ref_i = decomplex_dict["REF"]
        alt_i = decomplex_dict["ALT"]
        variant_i = copy(vcf_record)
        variant_i.position = vcf_record.position + offset
        variant_i.refbase = ref_i
        variant_i.altbase = alt_i
        snv_and_indel_records.append(variant_i)
    return snv_and_indel_records


def split_into_snv_and_indel(infile, snv_out, indel_out):
    with (
        genome.open_textfile(infile) as vcf_in,
        open(snv_out, "w") as snv_out,
        open(indel_out, "w") as indel_out,
    ):
        line_i = vcf_in.readline().rstrip()

        while line_i.startswith("#"):
            snv_out.write(line_i + "\n")
            indel_out.write(line_i + "\n")
            line_i = vcf_in.readline().rstrip()

        while line_i:
            vcf_i = genome.VCFVariantRecord.from_vcf_line(line_i)

            if ("," not in vcf_i.altbase) and ("/" not in vcf_i.altbase):
                if len(vcf_i.refbase) == 1 and len(vcf_i.altbase) == 1:
                    snv_out.write(line_i + "\n")
                elif len(vcf_i.refbase) == 1 or len(vcf_i.altbase) == 1:
                    indel_out.write(line_i + "\n")
                else:
                    snvs_and_indels = split_complex_variants_into_snvs_and_indels(vcf_i)
                    for snv_or_indel in snvs_and_indels:
                        if len(snv_or_indel.refbase) == len(snv_or_indel.altbase) == 1:
                            snv_out.write(snv_or_indel.to_vcf_line() + "\n")
                        else:
                            indel_out.write(snv_or_indel.to_vcf_line() + "\n")
            else:
                item = line_i.split("\t")
                if "," in vcf_i.altbase:
                    alt_bases = vcf_i.altbase.split(",")
                elif "/" in vcf_i.altbase:
                    alt_bases = vcf_i.altbase.split("/")
                else:
                    raise ValueError(f"Check the line: {line_i}")

                for ith_base, altbase_i in enumerate(alt_bases):
                    # snv
                    if len(vcf_i.refbase) == 1 and len(altbase_i) == 1:
                        item_j = copy(item)
                        item_j[4] = altbase_i
                        new_line = "\t".join(item_j)
                        snv_out.write(new_line + "\n")
                    # indel
                    elif len(vcf_i.refbase) == 1 or len(altbase_i) == 1:
                        item_j = copy(item)
                        item_j[4] = altbase_i
                        new_line = "\t".join(item_j)
                        indel_out.write(new_line + "\n")
                    # complex variants
                    else:
                        vcf_j = copy(vcf_i)
                        vcf_j.altbase = altbase_i
                        snvs_and_indels = split_complex_variants_into_snvs_and_indels(
                            vcf_j
                        )
                        for snv_or_indel in snvs_and_indels:
                            if (
                                len(snv_or_indel.refbase)
                                == len(snv_or_indel.altbase)
                                == 1
                            ):
                                snv_out.write(snv_or_indel.to_vcf_line() + "\n")
                            else:
                                indel_out.write(snv_or_indel.to_vcf_line() + "\n")

            line_i = vcf_in.readline().rstrip()


def main() -> None:
    infile, snv_out, indel_out = run()
    split_into_snv_and_indel(infile, snv_out, indel_out)


if __name__ == "__main__":
    main()
