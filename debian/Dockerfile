FROM debian:12
ARG timezone=Europe/Brussels

ENV TERM xterm-256color

# Bases
RUN apt -y update && apt -y upgrade \
    && apt -y install sagemath wget git curl build-essential libboost-all-dev

RUN cd /tmp && git clone https://github.com/Mister7F/docker-crypto

###################
# BreakingRainbow #
###################

# https://github.com/WardBeullens/BreakingRainbow
RUN cd /tmp/docker-crypto/sources/BreakingRainbow/xl-20160426 && make

# docker build . -t "crypto-debian"
# docker run -it --rm crypto-debian:latest
ENTRYPOINT bash
