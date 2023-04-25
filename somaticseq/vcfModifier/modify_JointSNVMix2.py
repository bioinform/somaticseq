#!/usr/bin/env python3

import argparse

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
        "-outfile", "--output-vcf", type=str, help="Output VCF file", required=True
    )

    # Parse the arguments:
    args = parser.parse_args()
    infile = args.input_vcf
    outfile = args.output_vcf

    return infile, outfile


def convert(infile, outfile):

    (
        idx_chrom,
        idx_pos,
        idx_id,
        idx_ref,
        idx_alt,
        idx_qual,
        idx_filter,
        idx_info,
        idx_format,
        idx_SM1,
        idx_SM2,
    ) = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    with genome.open_textfile(infile) as vcf, open(outfile, "w") as vcfout:

        line_i = vcf.readline().rstrip()

        # VCF header
        while line_i.startswith("#"):

            if line_i.startswith("##FORMAT=<ID=AD,"):
                line_i = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'

            vcfout.write(line_i + "\n")

            line_i = vcf.readline().rstrip()

        while line_i:

            item = line_i.split("\t")

            format_items = item[idx_format].split(":")
            if "AD" in format_items and "RD" in format_items:

                # NORMAL
                idx_ad = format_items.index("AD")
                idx_rd = format_items.index("RD")
                format_items.pop(idx_rd)

                item_normal = item[idx_SM1].split(":")
                normal_ad = int(item_normal[idx_ad])
                normal_rd = int(item_normal[idx_rd])

                try:
                    vaf = normal_ad / (normal_ad + normal_rd)
                except ZeroDivisionError:
                    vaf = 0

                if vaf > 0.8:
                    normal_gt = "1/1"
                elif vaf > 0.25:
                    normal_gt = "0/1"
                else:
                    normal_gt = "0/0"

                item_normal[idx_ad] = "{},{}".format(
                    item_normal[idx_rd], item_normal[idx_ad]
                )
                item_normal.pop(idx_rd)
                item_normal = [normal_gt] + item_normal

                # TUMOR
                item_tumor = item[idx_SM2].split(":")
                tumor_ad = int(item_tumor[idx_ad])
                tumor_rd = int(item_tumor[idx_rd])

                try:
                    vaf = tumor_ad / (tumor_ad + tumor_rd)
                except ZeroDivisionError:
                    vaf = 0

                if vaf > 0.8:
                    tumor_gt = "1/1"
                else:
                    tumor_gt = "0/1"

                item_tumor[idx_ad] = "{},{}".format(
                    item_tumor[idx_rd], item_tumor[idx_ad]
                )
                item_tumor.pop(idx_rd)
                item_tumor = [tumor_gt] + item_tumor

                # Rewrite
                item[idx_format] = "GT:" + ":".join(format_items)
                item[idx_SM1] = ":".join(item_normal)
                item[idx_SM2] = ":".join(item_tumor)

            line_i = "\t".join(item)

            vcfout.write(line_i + "\n")

            line_i = vcf.readline().rstrip()


if __name__ == "__main__":
    infile, outfile = run()
    convert(infile, outfile)
