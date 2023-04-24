FROM ubuntu:23.04
RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install sambamba
