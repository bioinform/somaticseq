FROM ubuntu:18.04

MAINTAINER Li Tai Fang <li_tai.fang@roche.com>

RUN apt-get update && apt-get install -y bwa samtools && apt-get clean
