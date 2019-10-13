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
    scripts=['somaticseq/run_somaticseq.py',
             'somaticseq_parallel.py',
             'utilities/dockered_pipelines/makeSomaticScripts.py',
             'utilities/lociCounterWithLabels.py',
             'utilities/split_Bed_into_equal_regions.py']
)
