#!/usr/bin/env python3
# mypy: ignore-errors
import os

import setuptools_scm  # noqa
from setuptools import find_packages, setup

with open("README.md") as fn:
    long_description = fn.read()

# Read the version from the _version.py file
version_file = os.path.join("somaticseq", "_version.py")
with open(version_file) as f:
    exec(f.read())  # This will define __version__

setup(
    name="somaticseq",
    description=(
        "SomaticSeq: "
        "An ensemble approach to accurately detect somatic mutations using SomaticSeq"
    ),
    version=__version__,  # noqa
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Li Tai Fang",
    author_email="ltfang@gmail.com",
    url="https://github.com/bioinform/somaticseq",
    packages=find_packages(),
    package_data={"": ["*.R"]},
    python_requires=">=3.11.0",
    setup_requires=["setuptools>=42", "setuptools_scm"],
    install_requires=[  # overridden by pyproject.toml
        "pysam",
        "numpy",
        "scipy",
        "pandas",
        "xgboost>=1.4",
        "pydantic>=2.0.0,<3.0",
    ],
    scripts=[
        "somaticseq/somaticseq_parallel.py",
        "somaticseq/run_somaticseq.py",
        "somaticseq/single_sample_vcf2tsv.py",
        "somaticseq/somatic_vcf2tsv.py",
        "somaticseq/somatic_xgboost.py",
        "somaticseq/somatic_tsv2vcf.py",
        "somaticseq/genomic_file_parsers/concat.py",
        "somaticseq/utilities/linguistic_sequence_complexity.py",
        "somaticseq/utilities/lociCounterWithLabels.py",
        "somaticseq/utilities/paired_end_bam2fastq.py",
        "somaticseq/utilities/remove_callers_from_somaticseq_tsv.py",
        "somaticseq/utilities/split_bed_into_equal_regions.py",
        "somaticseq/utilities/tally_variants_from_multiple_vcfs.py",
        "somaticseq/utilities/variant_annotation.py",
        "somaticseq/utilities/vcfsorter.pl",
        "somaticseq/utilities/dockered_pipelines/makeAlignmentScripts.py",
        "somaticseq/utilities/dockered_pipelines/makeSomaticScripts.py",
        "somaticseq/utilities/dockered_pipelines/run_workflows.py",
        "somaticseq/vcf_modifier/split_vcf.py",
        "r_scripts/ada_model_builder_ntChange.R",
        "r_scripts/ada_model_predictor.R",
    ],
)
