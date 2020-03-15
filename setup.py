#!/usr/bin/env python3

from setuptools import setup, find_packages
from somaticseq._version import  __version__ as version

print(version)

setup(
    name='SomaticSeq',
    version=version,
    description='SomaticSeq: An ensemble approach to accurately detect somatic mutations using SomaticSeq',
    author='Li Tai Fang',
    author_email='li_tai.fang@roche.com',
    url='https://github.com/bioinform/somaticseq',
    packages=find_packages(),
    package_data={'': ['*.R']},
    install_requires=['pysam', 'numpy', 'scipy'],
    scripts=['somaticseq_parallel.py',
             'somaticseq/run_somaticseq.py',
             'utilities/attach_pileupVAF.py',
             'utilities/bamQC.py',
             'utilities/lociCounterWithLabels.py',
             'utilities/remove_callers_from_somaticseq_tsv.py',
             'utilities/split_Bed_into_equal_regions.py',
             'utilities/tally_variants_from_multiple_vcfs.py',
             'utilities/variant_annotation.py',
             'utilities/dockered_pipelines/makeSomaticScripts.py',
             'somaticseq/single_sample_vcf2tsv.py',
             'somaticseq/somatic_vcf2tsv.py',
             'somaticseq/SSeq_tsv2vcf.py',]
)
