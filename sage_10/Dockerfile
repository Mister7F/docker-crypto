FROM sagemath/sagemath:10.2
ARG timezone=Europe/Brussels

ENV TERM xterm-256color

# Bases
RUN sudo apt-get -y update && sudo apt-get -y upgrade \
    && sudo apt-get -y install wget git curl build-essential libboost-all-dev

RUN cd /tmp && git clone https://github.com/Mister7F/docker-crypto

RUN sudo apt-get -y install python3-pip
RUN sage --pip install pycryptodome

# docker build . -t "crypto-sage_10"
# docker run -it --rm crypto-sage_10:latest
ENTRYPOINT bash
