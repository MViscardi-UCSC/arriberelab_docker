FROM ubuntu:focal

RUN apt update -y
RUN apt upgrade -y

ENV DEBIAN_FRONTEND=noninteractive
# This installs 3.5 with xenial, but 3.7 in bionic
RUN apt install python3-pip -y
RUN pip3 install --upgrade pip

# Install STAR for alignments
RUN apt install wget -y
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
RUN tar -xzf 2.7.10a.tar.gz
WORKDIR "STAR-2.7.10a/source"
RUN ls -a
RUN make STAR
WORKDIR "/"
ENV PATH="/STAR-2.7.10a/bin/Linux_x86_64/:${PATH}"

# Install Samtools:
RUN apt install samtools -y

RUN apt install nano -y
RUN apt install tree -y
RUN apt install git -y

RUN pip3 install pandas

RUN mkdir /usr/src/working
RUN mkdir /usr/src/working/scripts_dir
WORKDIR /usr/src/working/scripts_dir
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" /tmp/skipcache.txt
RUN git clone https://github.com/MViscardi-UCSC/arriberelab_docker
RUN chmod -R +x arriberelab_docker
ENV PATH="/usr/src/working/scripts_dir/arriberelab_docker:${PATH}"

WORKDIR /usr/src/working
