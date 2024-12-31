#!/usr/bin/env python3
# flake8: noqa: E501

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
        "-snv", "--output-snv", type=str, help="Output SNV VCF file", required=True
    )
    parser.add_argument(
        "-indel",
        "--output-indel",
        type=str,
        help="Output INDEL VCF file",
        required=True,
    )
    parser.add_argument(
        "-ref", "--genome-reference", type=str, help="genome reference", required=True
    )
    args = parser.parse_args()
    return args


def _make_new_vcf_line(
    vcfcall: genome.VCFVariantRecord,
    new_format_string: str,
    tumor_sample: str,
    normal_sample: str | None = None,
) -> str:
    assert vcfcall.chromosome is not None
    assert vcfcall.position is not None
    assert vcfcall.identifier is not None
    assert vcfcall.refbase is not None
    assert vcfcall.altbase is not None
    assert vcfcall.filters is not None
    assert vcfcall.info is not None
    if normal_sample:
        return "\t".join(
            (
                vcfcall.chromosome,
                str(vcfcall.position),
                vcfcall.identifier,
                vcfcall.refbase,
                vcfcall.altbase,
                str(vcfcall.qual or "."),
                vcfcall.filters,
                vcfcall.info,
                new_format_string,
                normal_sample,
                tumor_sample,
            )
        )
    return "\t".join(
        (
            vcfcall.chromosome,
            str(vcfcall.position),
            vcfcall.identifier,
            vcfcall.refbase,
            vcfcall.altbase,
            str(vcfcall.qual or "."),
            vcfcall.filters,
            vcfcall.info,
            new_format_string,
            tumor_sample,
        )
    )


def convert(infile, snv_out, indel_out, genome_reference):
    tempdir = tempfile.mkdtemp()
    tmp_snv_vcf = os.path.join(tempdir, uuid.uuid4().hex) + ".vcf"
    tmp_indel_vcf = os.path.join(tempdir, uuid.uuid4().hex) + ".vcf"
    with (
        genome.open_textfile(infile) as vcf,
        open(tmp_snv_vcf, "w") as snpout,
        open(tmp_indel_vcf, "w") as indelout,
    ):
        line_i = vcf.readline().rstrip()
        while line_i.startswith("##"):
            if re.match(r"^##INFO=<ID=(LSEQ|RSEQ),", line_i):
                line_i = line_i.replace("Number=G", "Number=1")

            elif line_i.startswith("##FORMAT=<ID=BIAS,"):
                line_i = line_i.replace("Number=1", "Number=.")

            elif (
                line_i.startswith("##FORMAT=<ID=PSTD,")
                or line_i.startswith("##FORMAT=<ID=QSTD,")
                or line_i.startswith("##INFO=<ID=SOR,")
            ):
                line_i = line_i.replace("Type=Float", "Type=String")

            snpout.write(line_i + "\n")
            indelout.write(line_i + "\n")
            line_i = vcf.readline().rstrip()

        addition_header = []
        addition_header.append(
            '##INFO=<ID=Germline,Number=0,Type=Flag,Description="VarDict Germline">'
        )
        addition_header.append(
            '##INFO=<ID=StrongSomatic,Number=0,Type=Flag,Description="VarDict Strong Somatic">'
        )
        addition_header.append(
            '##INFO=<ID=LikelySomatic,Number=0,Type=Flag,Description="VarDict Likely Somatic">'
        )
        addition_header.append(
            '##INFO=<ID=LikelyLOH,Number=0,Type=Flag,Description="VarDict Likely LOH">'
        )
        addition_header.append(
            '##INFO=<ID=StrongLOH,Number=0,Type=Flag,Description="VarDict Strong LOH">'
        )
        addition_header.append(
            '##INFO=<ID=AFDiff,Number=0,Type=Flag,Description="VarDict AF Diff">'
        )
        addition_header.append(
            '##INFO=<ID=Deletion,Number=0,Type=Flag,Description="VarDict Deletion">'
        )
        addition_header.append(
            '##INFO=<ID=SampleSpecific,Number=0,Type=Flag,Description="VarDict SampleSpecific">'
        )
        addition_header.append(
            '##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">'
        )
        for item_i in addition_header:
            snpout.write(item_i + "\n")
            indelout.write(item_i + "\n")

        # This is the #CHROM line
        header_main_item = line_i.split("\t")
        num_header = len(header_main_item)

        if num_header == 10:
            paired = False
        elif num_header == 11:
            paired = True
        else:
            raise ValueError("New VarDict VCF file SomaticSeq did not handle properly.")

        snpout.write(line_i + "\n")
        indelout.write(line_i + "\n")

        line_i = vcf.readline().rstrip()
        while line_i:
            vcfcall = genome.VCFVariantRecord.from_vcf_line(line_i)

            # Fix the occasional error where ALT and REF are the same:
            if vcfcall.refbase == vcfcall.altbase:
                continue

            # In the REF/ALT field, non-GCTA characters should be changed to N
            # to fit the VCF standard:
            vcfcall.refbase = re.sub(r"[^GCTA]", "N", vcfcall.refbase, flags=re.I)
            vcfcall.altbase = re.sub(r"[^GCTA]", "N", vcfcall.altbase, flags=re.I)

            ## To be consistent with other tools, Combine AD:RD or ALD:RD into DP4.
            # VarDict puts Tumor first and Normal next Also, the old version has
            # no ALD (somatic.pl). The new version has ALD (paired.pl).
            format_field = vcfcall.field.split(":")
            idx_rd = format_field.index("RD")
            tumor_sample = vcfcall.samples[0].split(":")
            tumor_dp4 = tumor_sample.pop(idx_rd)
            if paired:
                normal_sample = vcfcall.samples[1].split(":")
                normal_dp4 = normal_sample.pop(idx_rd)

            format_field.pop(idx_rd)

            # As right now, the old version has no ALD. The new version has ALD.
            # If the VCF has no ALD, then the AD means the same thing ALD is
            # supposed to mean.
            try:
                idx_ad = format_field.index("ALD")
            except ValueError:
                idx_ad = format_field.index("AD")

            if paired:
                normal_dp4 = normal_dp4 + "," + normal_sample.pop(idx_ad)

            tumor_dp4 = tumor_dp4 + "," + tumor_sample.pop(idx_ad)
            format_field.pop(idx_ad)
            # Re-format the strings:
            format_field.append("DP4")

            if paired:
                normal_sample.append(normal_dp4)
            tumor_sample.append(tumor_dp4)

            if paired:
                normal_sample = ":".join(normal_sample)
            tumor_sample = ":".join(tumor_sample)
            new_format_string = ":".join(format_field)

            # VarDict's END tag has caused problem with GATK CombineVariants.
            # Simply get rid of it.
            vcfcall.info = re.sub(r"END=[0-9]+;", "", vcfcall.info)

            if paired:
                line_i = _make_new_vcf_line(
                    vcfcall=vcfcall,
                    new_format_string=new_format_string,
                    tumor_sample=tumor_sample,
                    normal_sample=normal_sample,
                )
            else:
                line_i = _make_new_vcf_line(
                    vcfcall=vcfcall,
                    new_format_string=new_format_string,
                    tumor_sample=tumor_sample,
                )
            # Write to snp and indel into different files:
            if "TYPE=SNV" in vcfcall.info:
                snpout.write(line_i + "\n")

            elif "TYPE=Deletion" in vcfcall.info or "TYPE=Insertion" in vcfcall.info:
                indelout.write(line_i + "\n")

            else:  # "TYPE=Complex"
                complex_call = genome.VCFVariantRecord.from_vcf_line(line_i)
                snvs_and_indels = split_complex_variants_into_snvs_and_indels(
                    complex_call
                )
                for snv_or_indel in snvs_and_indels:
                    if len(snv_or_indel.refbase) == len(snv_or_indel.altbase) == 1:
                        snpout.write(snv_or_indel.to_vcf_line() + "\n")
                    else:
                        indelout.write(snv_or_indel.to_vcf_line() + "\n")

            # Continue:
            line_i = vcf.readline().rstrip()

    vcfsorter(genome_reference, tmp_snv_vcf, snv_out)
    vcfsorter(genome_reference, tmp_indel_vcf, indel_out)
    os.remove(tmp_snv_vcf)
    os.remove(tmp_indel_vcf)


if __name__ == "__main__":
    args = run()
    convert(args.input_vcf, args.output_snv, args.output_indel, args.genome_reference)
