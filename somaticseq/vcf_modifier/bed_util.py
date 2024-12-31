import os
import re
import subprocess
import uuid

import somaticseq.genomic_file_parsers.genomic_file_handlers as genome


def remove_vcf_illegal_lines(invcf, outvcf):
    """
    In VarDict v1.7, there are lines with <XXX> in ALT without END in info,
    which will cause bedtools to fail. This program will check if these things
    exist, and if they do, remove them. If the input VCF has illegal lines, it
    will return the modified output VCF file excluding those lines. If the input
    VCF file does not have such illegal lines, it will return False.
    """

    hasIllegalLine = False
    with genome.open_textfile(invcf) as vcf:
        line_i = vcf.readline().rstrip()
        while line_i.startswith("#"):
            line_i = vcf.readline().rstrip()

        while line_i:
            vcf_i = genome.VCFVariantRecord.from_vcf_line(line_i)

            if re.match(r"<\w+>", vcf_i.altbase) and (not vcf_i.get_info_value("END")):
                hasIllegalLine = True
                break

            line_i = vcf.readline().rstrip()

    if hasIllegalLine:
        with genome.open_textfile(invcf) as vcf, open(outvcf, "w") as out:
            line_i = vcf.readline().rstrip()
            while line_i.startswith("#"):
                out.write(line_i + "\n")
                line_i = vcf.readline().rstrip()

            while line_i:
                vcf_i = genome.VCFVariantRecord.from_vcf_line(line_i)

                if not (
                    re.match(r"<\w+>", vcf_i.altbase)
                    and (not vcf_i.get_info_value("END"))
                ):
                    out.write(line_i + "\n")

                line_i = vcf.readline().rstrip()

        return outvcf

    else:
        return hasIllegalLine


def bed_include(infile, inclusion_region, outfile):
    assert infile != outfile

    if inclusion_region:
        cmd_line = "bedtools intersect -header -a {} -b {} | uniq > {}".format(
            infile, inclusion_region, outfile
        )
        subprocess.check_call(cmd_line, shell=True)

    else:
        outfile = None

    return outfile


def bed_exclude(infile, exclusion_region, outfile):
    assert infile != outfile

    if exclusion_region:
        cmd_line = "bedtools intersect -header -a {} -b {} -v | uniq > {}".format(
            infile, exclusion_region, outfile
        )
        subprocess.check_call(cmd_line, shell=True)

    else:
        outfile = None

    return outfile


def bed_intersector(infile, outfile, inclusion_region=None, exclusion_region=None):
    assert infile != outfile
    from shutil import copyfile

    # Get the input file name minus the extention, and also get the extension
    infile_noext = re.sub(r"\.\w+$", "", infile)
    file_ext = re.search(r"\.\w+$", infile).group()

    temp_files = []

    if inclusion_region:
        included_temp_file = infile_noext + uuid.uuid4().hex + file_ext

        cmd_line = "bedtools intersect -header -a {} -b {} | uniq > {}".format(
            infile, inclusion_region, included_temp_file
        )
        subprocess.check_call(cmd_line, shell=True)

        infile = included_temp_file
        temp_files.append(included_temp_file)

    if exclusion_region:
        cmd_line = "bedtools intersect -header -a {} -b {} -v | uniq > {}".format(
            infile, exclusion_region, outfile
        )
        subprocess.check_call(cmd_line, shell=True)

    if inclusion_region and not exclusion_region:
        copyfile(included_temp_file, outfile)

    elif not (inclusion_region or exclusion_region):
        if infile.endswith(".gz"):
            cmd_line = f"gunzip -c {infile} > {outfile}"
            subprocess.check_call(cmd_line, shell=True)

        else:
            copyfile(infile, outfile)

    for file_i in temp_files:
        os.remove(file_i)

    return outfile


def bed_sort_and_merge(infile: str, outfile: str, fai: str) -> None:
    cmd_line = (
        f"bedtools sort -faidx {fai} -header -i {infile} | "
        f"bedtools merge -i stdin > {outfile}"
    )
    subprocess.check_call(cmd_line, shell=True)


def vcfsorter(ref: str, vcfin: str, vcfout: str):
    """
    Uses bedtools sort and then make sure there is no duplicate line.

    Args:
        ref: Fasta file. Code assumes .fai will exist.
        vcfin: vcf file input
        vcfout: vcf file output

    Returns:
        Exit code for bedtools sort.
    """
    fai = ref + ".fai"
    cmd_line = f"bedtools sort -faidx {fai} -header -i {vcfin} | uniq > {vcfout}"
    subprocess.check_call(cmd_line, shell=True)
