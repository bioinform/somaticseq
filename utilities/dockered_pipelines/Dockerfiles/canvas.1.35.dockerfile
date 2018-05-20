FROM ubuntu:16.04

RUN apt-get update && apt-get install -y wget mono-runtime mono-complete tar libunwind-dev && apt-get clean
RUN cd /opt && wget https://download.microsoft.com/download/D/7/A/D7A9E4E9-5D25-4F0C-B071-210CB8267943/dotnet-ubuntu.16.04-x64.1.1.2.tar.gz && tar -xvf dotnet-ubuntu.16.04-x64.1.1.2.tar.gz && ln -s /opt/shared/Microsoft.NETCore.App/1.1.2/dotnet /usr/bin/dotnet
RUN cd /opt && wget https://github.com/Illumina/canvas/releases/download/1.35.1.1316%2Bmaster/Canvas-1.35.1.1316.master_x64.tar.gz && tar -xvf Canvas-1.35.1.1316.master_x64.tar.gz && ln -s 'Canvas-1.35.1.1316+master_x64/' Canvas && chmod a+x Canvas/tabix
