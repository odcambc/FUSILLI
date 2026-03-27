# Container Definitions for FUSILLI

This directory contains container definitions for running FUSILLI in isolated environments.

## Supported Container Runtimes

- **Docker**: Local development, CI/CD
- **Singularity/Apptainer**: HPC clusters (SLURM, PBS, etc.)

## Quick Start

### Docker

```bash
cd containers

# Build image
./build_docker.sh

# Run pipeline with Docker
snakemake -s workflow/Snakefile \
  --software-deployment-method conda \
  --conda-frontend docker \
  --docker-image fusilli:latest
```

### Singularity

```bash
cd containers

# Build SIF image
./build_singularity.sh

# Run pipeline with Singularity
snakemake -s workflow/Snakefile \
  --software-deployment-method apptainer \
  --apptainer-args "--bind /data:/data" \
  --apptainer-prefix ./containers
```

## Image Details

### Base Image

- **Docker**: `mambaorg/micromamba:1.5`
- **Singularity**: Same base, converted via definition file

### Installed Software

See `fusilli_env.yaml` for the full list of dependencies installed via micromamba/conda.

## Customization

### Adding Custom Resources

To include custom adapter files or contaminant sequences in the container:

1. Edit the Dockerfile or fusilli.def
2. Add `%files` section:

```singularity
%files
    /path/to/custom/adapters.fa /workspace/adapters.fa
```

### Environment Variables

- `MAMBA_ROOT_PREFIX=/opt/conda` - Conda package directory
- `PATH=$MAMBA_ROOT_PREFIX/bin:$PATH` - Executable search path

## Testing Containers

```bash
# Test Docker container
docker run --rm fusilli:latest snakemake --version

# Test Singularity container
singularity exec fusilli_dev.sif snakemake --version
```
