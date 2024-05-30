FROM ubuntu:24.04
RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install samtools && apt clean
