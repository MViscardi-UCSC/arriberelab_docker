FROM ubuntu:focal

RUN apt update -y
RUN apt upgrade -y

ENV DEBIAN_FRONTEND=noninteractive
# This installs 3.5 with xenial, but 3.7 in bionic
RUN apt install python3-pip -y
RUN pip3 install --upgrade pip

# Install STAR for alignment
RUN apt install wget -y
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
RUN tar -xzf 2.7.10a.tar.gz
WORKDIR "STAR-2.7.10a/source"
RUN ls -a
RUN make STAR
WORKDIR "/"

# Install Samtools:
RUN apt install samtools -y

# Bash starts up when using 'docker run -it <name>'
RUN apt install nano
RUN mkdir /usr/src/working
WORKDIR /usr/src/working

RUN apt install tree

RUN wget 

# COPY docker_script.sh /bin/startup
CMD /bin/bash
