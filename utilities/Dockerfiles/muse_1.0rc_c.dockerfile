FROM ubuntu:20.04

RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install wget
RUN cd /usr/local/bin/ && wget http://bioinformatics.mdanderson.org/Software/MuSE/MuSEv1.0rc_submission_c039ffa && chmod a+x MuSEv1.0rc_submission_c039ffa
