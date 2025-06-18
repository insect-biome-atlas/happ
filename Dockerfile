FROM ubuntu:latest

LABEL description="This is a Dockerfile recipe for the HAPP workflow."
LABEL maintainer="John Sundh (NBIS)"
ARG TARGETPLATFORM
ARG TARGETARCH
ARG TARGETVARIANT
ARG TARGETOS

WORKDIR /work

SHELL ["/bin/bash", "-c"]

ENV TZ=Europe/London

RUN echo "I'm building for $TARGETARCH on $TARGETOS"

RUN apt update && apt install -y wget software-properties-common tzdata && \
    dpkg-reconfigure -f noninteractive tzdata && \
    add-apt-repository -y ppa:apptainer/ppa && apt update && apt install -y apptainer && \
    apt-get clean

RUN wget https://pixi.sh/install.sh && bash install.sh && rm install.sh

RUN mkdir -p data

COPY pixi.toml environment.yml ./

RUN $HOME/.pixi/bin/pixi install -e default && $HOME/.pixi/bin/pixi clean cache -y && echo 'eval "$($HOME/.pixi/bin/pixi shell-hook)"' >> $HOME/.bashrc

RUN eval "$($HOME/.pixi/bin/pixi shell-hook)" && conda config --add channels defaults && conda config --set channel_priority strict

COPY config config
COPY data/test data/test
COPY dardel dardel
COPY test test
COPY local local
COPY slurm slurm
COPY workflow workflow

CMD ["/bin/bash"]