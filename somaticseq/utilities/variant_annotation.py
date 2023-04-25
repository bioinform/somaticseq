#!/usr/bin/env python3

import argparse
import gzip
import itertools
import logging
import multiprocessing
import os
import re
import subprocess
import tempfile
import uuid

import pysam

COSMIC_STRING = "GENE,CDS,AA,CNT"
DBSNP_STRING = (
    "RSPOS,GENEINFO,dbSNPBuildID,SAO,SSR,VC,PM,MUT,KGPhase1,KGPhase3,OM,CDA,CAF,COMMON"
)


def snpsift_snp(snpsift_jar, input_vcf, dbsnp_vcf, output_vcf, info_string):

    logger = logging.getLogger(snpsift_snp.__name__)
    sift_command = "java -Xmx8g -jar {} annotate -info {} {} {} > {}".format(
        snpsift_jar, info_string, dbsnp_vcf, input_vcf, output_vcf
    )
    logger.info(sift_command)
    subprocess.check_call(sift_command, shell=True)

    return output_vcf


def snpsift_cosmic(snpsift_jar, input_vcf, cosmic_vcf, output_vcf, info_string):

    logger = logging.getLogger(snpsift_cosmic.__name__)
    sift_command = "java -Xmx8g -jar {} annotate -info {} {} {} > {}".format(
        snpsift_jar, info_string, cosmic_vcf, input_vcf, output_vcf
    )
    logger.info(sift_command)
    subprocess.check_call(sift_command, shell=True)

    return output_vcf


def snpeff_annotate(snpeff_jar, input_vcf, output_vcf, db):

    logger = logging.getLogger(snpeff_annotate.__name__)
    eff_command = "java -Xmx8g -jar {} -noStats {} {} > {}".format(
        snpeff_jar, db, input_vcf, output_vcf
    )
    logger.info(eff_command)
    subprocess.check_call(eff_command, shell=True)

    return output_vcf


def annotate_small_variants(
    snpsift_jar,
    snpeff_jar,
    input_vcf,
    dbsnp_vcf,
    cosmic_vcf,
    output_vcf,
    snp_string,
    cosmic_string,
    eff_db,
):

    dirname = tempfile.gettempdir()

    dbsnp_annotated = snpsift_snp(
        snpsift_jar,
        input_vcf,
        dbsnp_vcf,
        os.path.join(dirname, uuid.uuid4().hex + ".vcf"),
        snp_string,
    )
    cosmic_annotated = snpsift_cosmic(
        snpsift_jar,
        dbsnp_annotated,
        cosmic_vcf,
        os.path.join(dirname, uuid.uuid4().hex + ".vcf"),
        cosmic_string,
    )
    output_vcf = snpeff_annotate(snpeff_jar, cosmic_annotated, output_vcf, eff_db)

    os.remove(dbsnp_annotated)
    os.remove(cosmic_annotated)

    pysam.tabix_index(output_vcf, force=True, preset="vcf")

    return output_vcf + ".gz"


if __name__ == "__main__":

    FORMAT = "%(levelname)s %(asctime)-15s %(name)-20s %(message)s"
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    parser = argparse.ArgumentParser(
        description="Annotate with snpSift and snpEff with dbSNP and COSMIC",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-infile", "--infile", help="input vcf file")
    parser.add_argument("-outfile", "--outfile", help="output vcf file")
    parser.add_argument(
        "-dbsnp", "--dbsnp", help="dbsnp vcf file to feed into GATK4 HaplotypeCaller"
    )
    parser.add_argument(
        "-cosmic", "--cosmic", help="cosmic vcf file to feed into GATK4 HaplotypeCaller"
    )
    parser.add_argument("-snpsift", "--snpsift", help="SnpSift JAR")
    parser.add_argument("-snpeff", "--snpeff", help="snpEff JAR")
    parser.add_argument("-db", "--snpeff-db", help="snpEff db", default="GRCh38.86")

    args = parser.parse_args()

    annotate_small_variants(
        args.snpsift,
        args.snpeff,
        args.infile,
        args.dbsnp,
        args.cosmic,
        args.outfile,
        DBSNP_STRING,
        COSMIC_STRING,
        args.snpeff_db,
    )
