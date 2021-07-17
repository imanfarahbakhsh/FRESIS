# hash:sha256:afebed87ce16b6f26225c9eb2121d10ac044597edc4739773ef00cbe19c8d6ff
FROM registry.codeocean.com/codeocean/ubuntu:20.04.2

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gfortran \
    && rm -rf /var/lib/apt/lists/*
