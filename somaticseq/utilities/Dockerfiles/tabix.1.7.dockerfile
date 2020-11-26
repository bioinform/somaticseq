FROM ubuntu:18.04
RUN apt update && apt install -y tabix && apt clean
