FROM condaforge/mambaforge:latest

ARG TZ=America/Chicago
ENV TZ=${TZ}

# System packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    zsh \
    sudo \
    tzdata \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=${USER_UID}
RUN groupadd --gid ${USER_GID} ${USERNAME} \
    && useradd --uid ${USER_UID} --gid ${USER_GID} -m ${USERNAME} -s /bin/zsh \
    && echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Copy environment spec and create conda environment
COPY fusilli_env.yaml /tmp/fusilli_env.yaml
RUN mamba env create -f /tmp/fusilli_env.yaml && mamba clean -afy

# Make the conda environment activate by default for all shells
RUN echo "conda activate FUSILLI" >> /home/${USERNAME}/.bashrc \
    && echo "conda activate FUSILLI" >> /home/${USERNAME}/.zshrc

# Install pytest in the environment (for running tests)
RUN mamba run -n FUSILLI pip install pytest --no-cache-dir

WORKDIR /workspace
USER ${USERNAME}
