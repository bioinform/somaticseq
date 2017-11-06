FROM ubuntu:14.04

RUN apt-get update && apt-get install -y python python-dev wget build-essential samtools zlib1g-dev cython && apt-get clean
RUN cd /opt && wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/pysam/pysam-0.5.tar.gz && tar -xvf pysam-0.5.tar.gz
RUN cd /opt/pysam-0.5 && wget https://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11-py2.7.egg && python setup.py build && python setup.py install

RUN cd /opt && wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/joint-snv-mix/JointSNVMix-0.7.5.tar.gz && tar -xvf JointSNVMix-0.7.5.tar.gz && cd JointSNVMix-0.7.5 && python setup.py install
RUN cd /opt && wget https://www.dropbox.com/s/rbegan3opz2fc4k/vcfsorter.pl && chmod a+x vcfsorter.pl
