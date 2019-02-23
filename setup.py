#!/usr/bin/env python3

from setuptools import setup, find_packages

version = 'Unknown'
for line in open('VERSION'):
    if line.startswith('##SomaticSeq='):
        version = line.strip().split("=")[1].strip()

print(version)

setup(
    name='SomaticSeq',
    version=version,
    description='SomaticSeq: An ensemble approach to accurately detect somatic mutations using SomaticSeq',
    author='Li Tai Fang',
    author_email='li_tai.fang@roche.com',
    url='https://github.com/bioinform/somaticseq',
    packages=find_packages(),
    install_requires=['pysam', 'numpy', 'scipy'],
    scripts=['somaticseq/run_somaticseq.py',
             'somaticseq_parallel.py',
             'utilities/dockered_pipelines/makeSomaticScripts.py',],
)
