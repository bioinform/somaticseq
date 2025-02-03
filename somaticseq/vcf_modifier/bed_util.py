import os
import re
import tempfile

from pybedtools import BedTool

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


def bed_include(infile: str, inclusion_region: str | None, outfile: str) -> str:
    assert infile != outfile
    if not inclusion_region:
        return infile

    # Load the BED files
    a = BedTool(infile)
    b = BedTool(inclusion_region)

    # Perform the intersection with header
    intersected = a.intersect(b, header=True)
    _, temp_file = tempfile.mkstemp(text=True)
    out = intersected.saveas(temp_file)
    genome.uniq(out.fn, outfile)
    os.remove(out.fn)
    return outfile


def bed_exclude(infile: str, exclusion_region: str | None, outfile: str) -> str:
    assert infile != outfile
    if not exclusion_region:
        return infile

    # Load the BED files
    a = BedTool(infile)
    b = BedTool(exclusion_region)

    # Perform the intersection with header. 'v' is exclude
    excluded = a.intersect(b, v=True, header=True)
    _, temp_file = tempfile.mkstemp(text=True)
    out = excluded.saveas(temp_file)
    genome.uniq(out.fn, outfile)
    os.remove(out.fn)
    return outfile


def bed_intersector(
    infile: str,
    outfile: str,
    inclusion_region: str | None = None,
    exclusion_region: str | None = None,
) -> str:
    assert infile != outfile

    # Get the input file name minus the extention, and also get the extension
    infile_noext, file_ext = os.path.splitext(infile)

    _, intermediate_file = tempfile.mkstemp(text=True)
    included = bed_include(infile, inclusion_region, intermediate_file)
    included_then_exclude = bed_exclude(included, exclusion_region, outfile)
    return included_then_exclude


def bed_sort_and_merge(infile: str, outfile: str, fai: str) -> None:
    sorted_bed = BedTool(infile).sort(faidx=fai, header=True)
    merged_bed = sorted_bed.merge()
    merged_bed.saveas(outfile)


def vcfsorter(ref: str, vcfin: str, vcfout: str) -> None:
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

    _, intermediate_file = tempfile.mkstemp(text=True)
    sorted_bed = BedTool(vcfin).sort(faidx=fai, header=True)
    intermediate_bed = sorted_bed.saveas(intermediate_file)
    genome.uniq(intermediate_bed.fn, vcfout)
    os.remove(intermediate_bed.fn)
