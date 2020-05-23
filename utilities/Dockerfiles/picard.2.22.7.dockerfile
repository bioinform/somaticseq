FROM ubuntu:18.04

RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install openjdk-8-jdk wget && apt clean
RUN cd /opt && wget https://github.com/broadinstitute/picard/releases/download/2.22.7/picard.jar
RUN apt -y autoremove wget
