FROM ubuntu:16.04

RUN apt-get update && apt-get -y install wget tar mono-runtime mono-complete libunwind-dev libcurl3 libssl1.0.0 libssl-dev && apt-get clean
RUN cd /opt && wget https://download.microsoft.com/download/2/E/C/2EC018A0-A0FC-40A2-849D-AA692F68349E/dotnet-sdk-2.1.105-linux-x64.tar.gz && tar -xvf dotnet-sdk-2.1.105-linux-x64.tar.gz && ln -s /opt/dotnet /usr/local/bin/dotnet
RUN cd /opt && wget https://github.com/Illumina/Nirvana/archive/v2.0.9.tar.gz && tar -xvf v2.0.9.tar.gz && ln -s Nirvana-2.0.9 Nirvana && cd Nirvana-2.0.9 && /opt/dotnet build -c Release
RUN chmod a+rx /root/
