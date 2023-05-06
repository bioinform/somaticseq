#!/usr/bin/env python3

import argparse
import re

import somaticseq.genomicFileHandler.genomic_file_handlers as genome


def run():
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

    args = parser.parse_args()
    infile = args.input_vcf
    snv_out = args.output_snv
    indel_out = args.output_indel

    return infile, snv_out, indel_out


def convert(infile, snv_out, indel_out):

    with genome.open_textfile(infile) as vcf, open(snv_out, "w") as snpout, open(
        indel_out, "w"
    ) as indelout:

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

        snpout.write(line_i + "\n")
        indelout.write(line_i + "\n")

        line_i = vcf.readline().rstrip()
        while line_i:

            vcfcall = genome.VcfLine(line_i)

            # Fix the occasional error where ALT and REF are the same:
            if vcfcall.refbase != vcfcall.altbase:

                # In the REF/ALT field, non-GCTA characters should be changed to N to fit the VCF standard:
                vcfcall.refbase = re.sub(r"[^GCTA]", "N", vcfcall.refbase, flags=re.I)
                vcfcall.altbase = re.sub(r"[^GCTA]", "N", vcfcall.altbase, flags=re.I)

                ## To be consistent with other tools, Combine AD:RD or ALD:RD into DP4.
                # VarDict puts Tumor first and Normal next
                # Also, the old version has no ALD (somatic.pl). The new version has ALD (paired.pl).
                format_field = vcfcall.field.split(":")
                idx_rd = format_field.index("RD")

                tumor_sample = vcfcall.samples[0].split(":")
                tumor_dp4 = tumor_sample.pop(idx_rd)

                if paired:
                    normal_sample = vcfcall.samples[1].split(":")
                    normal_dp4 = normal_sample.pop(idx_rd)

                format_field.pop(idx_rd)

                # As right now, the old version has no ALD. The new version has ALD.
                # If the VCF has no ALD, then the AD means the same thing ALD is supposed to mean.
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

                # VarDict's END tag has caused problem with GATK CombineVariants. Simply get rid of it.
                vcfcall.info = re.sub(r"END=[0-9]+;", "", vcfcall.info)

                if paired:
                    line_i = "\t".join(
                        (
                            vcfcall.chromosome,
                            str(vcfcall.position),
                            vcfcall.identifier,
                            vcfcall.refbase,
                            vcfcall.altbase,
                            vcfcall.qual,
                            vcfcall.filters,
                            vcfcall.info,
                            new_format_string,
                            normal_sample,
                            tumor_sample,
                        )
                    )
                else:
                    line_i = "\t".join(
                        (
                            vcfcall.chromosome,
                            str(vcfcall.position),
                            vcfcall.identifier,
                            vcfcall.refbase,
                            vcfcall.altbase,
                            vcfcall.qual,
                            vcfcall.filters,
                            vcfcall.info,
                            new_format_string,
                            tumor_sample,
                        )
                    )

                # Write to snp and indel into different files:
                if "TYPE=SNV" in vcfcall.info:
                    snpout.write(line_i + "\n")

                elif (
                    "TYPE=Deletion" in vcfcall.info or "TYPE=Insertion" in vcfcall.info
                ):
                    indelout.write(line_i + "\n")

                elif "TYPE=Complex" in vcfcall.info and (
                    len(vcfcall.refbase) == len(vcfcall.altbase)
                ):
                    i = 0

                    for ref_i, alt_i in zip(vcfcall.refbase, vcfcall.altbase):

                        if ref_i != alt_i:
                            if paired:
                                line_i = "\t".join(
                                    (
                                        vcfcall.chromosome,
                                        str(vcfcall.position + i),
                                        vcfcall.identifier,
                                        ref_i,
                                        alt_i,
                                        vcfcall.qual,
                                        vcfcall.filters,
                                        vcfcall.info,
                                        new_format_string,
                                        normal_sample,
                                        tumor_sample,
                                    )
                                )
                            else:
                                line_i = "\t".join(
                                    (
                                        vcfcall.chromosome,
                                        str(vcfcall.position + i),
                                        vcfcall.identifier,
                                        ref_i,
                                        alt_i,
                                        vcfcall.qual,
                                        vcfcall.filters,
                                        vcfcall.info,
                                        new_format_string,
                                        tumor_sample,
                                    )
                                )

                            snpout.write(line_i + "\n")

                        i += 1

            # Continue:
            line_i = vcf.readline().rstrip()


if __name__ == "__main__":
    infile, snv_out, indel_out = run()
    convert(infile, snv_out, indel_out)
