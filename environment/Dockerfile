# hash:sha256:b511f80f2e74eb2312f551f3c3bf5ba5ec2c62df75233eef8fe5ae1c0c603f9b
FROM registry.codeocean.com/codeocean/ubuntu:20.04.2

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gfortran=4:9.3.0-1ubuntu2 \
    && rm -rf /var/lib/apt/lists/*
