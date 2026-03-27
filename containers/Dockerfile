FROM mambaorg/micromamba:1.5 as base

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

WORKDIR /workspace

COPY fusilli_env.yaml /tmp/fusilli_env.yaml

RUN micromamba install -y -n base -f /tmp/fusilli_env.yaml && \
    micromamba clean --all --yes && \
    mkdir -p /workspace

CMD ["/bin/bash"]
