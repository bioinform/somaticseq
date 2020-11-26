FROM ubuntu:16.04
 
RUN apt-get update && apt-get install -y build-essential git-core cmake zlib1g-dev libncurses-dev
RUN cd /opt/ && git clone https://github.com/genome/somatic-sniper.git && mkdir -p /opt/somatic-sniper/build && cd /opt/somatic-sniper/build && cmake ../ && make deps && make -j
