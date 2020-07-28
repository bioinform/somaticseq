FROM ubuntu:20.04

RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install wget unzip default-jdk && apt clean
RUN cd /opt && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip && ln -s Trimmomatic-0.39 Trimmomatic && cd Trimmomatic && ln -s trimmomatic-0.39.jar trimmomatic.jar
