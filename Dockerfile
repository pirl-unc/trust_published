# trust_published:1
# this build is on a platform that uses aquasecurity: the image is built as root but intended to be run as the user
# to run this image as root get rid of the --user flags

FROM debian:stretch

USER root

# install tools
RUN \
  apt-get update && \
  apt-get install -yq \
    less \
    nano \
    sudo \
    wget \ 
    unzip

# add python
# pip 10 has a bug and shouldn't be used as of 20180614.
RUN \
  apt-get -yq install \
    python-dev \
    build-essential \
    python-pip && \
 pip install --upgrade pip==9.0.3 && \
 pip install \
  Biopython \
  iterutils \
  numpy \
  pysam \
  argparse && \
 apt-get update && \
 apt-get install python-gmpy

RUN \
  pip install parasail

# now add trust
RUN \
  cd /opt && \
  wget https://media.nature.com/original/nature-assets/ng/journal/v49/n4/extref/ng.3820-s2.zip && \
  unzip -o ng.3820-s2.zip && \
  rm -r __MACOSX && \
  rm ng.3820-s2.zip && \
  mv SupplementarySoftware TRUST

RUN \ 
  apt-get clean  

# Need to set the home so trust can write and avoid this error :
#  Permission denied: '/root/.cache'
ENV \ 
  HOME=/tmp

CMD \
  bash -c "trust --help"
