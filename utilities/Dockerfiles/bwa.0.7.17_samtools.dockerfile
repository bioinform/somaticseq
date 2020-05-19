FROM ubuntu:20.04

MAINTAINER Li Tai Fang <li_tai.fang@roche.com>

RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install bwa samtools && apt-get clean
