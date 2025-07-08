FROM ghcr.io/prefix-dev/pixi:latest AS build

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

COPY pixi.lock pixi.toml environment.yml ./

RUN pixi install -e default && pixi clean cache -y
# Create the shell-hook bash script to activate the environment
RUN pixi shell-hook -e default > /shell-hook.sh
# extend the shell-hook script to run the command passed to the container
RUN echo 'exec "$@"' >> /shell-hook.sh
# configure conda
RUN eval "$(pixi shell-hook)" && conda config --add channels defaults && conda config --set channel_priority strict

FROM ubuntu:24.04 AS production

SHELL ["/bin/bash", "-c"]

ENV TZ=Europe/London

RUN echo "I'm building for $TARGETARCH on $TARGETOS"

RUN apt update && apt install -y wget software-properties-common tzdata && \
    dpkg-reconfigure -f noninteractive tzdata && \
    add-apt-repository -y ppa:apptainer/ppa && apt update && apt install -y apptainer && \
    apt-get clean

COPY --from=build /work/.pixi/envs/default /work/.pixi/envs/default
COPY --from=build /shell-hook.sh /shell-hook.sh
WORKDIR /work

RUN mkdir -p data
RUN mkdir -p /tmp/qiime2

COPY config config
COPY data/test data/test
COPY data/neeat_test data/neeat_test
COPY dardel dardel
COPY test test
COPY local local
COPY slurm slurm
COPY workflow workflow

ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]

CMD ["/bin/bash"]