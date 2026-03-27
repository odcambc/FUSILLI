# Snakemake Profiles for FUSILLI

This directory contains Snakemake profiles for different execution environments.

## Available Profiles

### `slurm/`

SLURM cluster profile for High Performance Computing (HPC) environments.

**Features:**
- Job submission via `sbatch`
- Status monitoring via `sacct`/`squeue`
- Configurable partitions and resources
- Automatic job cancellation
- Retry on failure (3 attempts)
- Filesystem sync latency handling (60s)

**Configuration:**
Edit `slurm/config.yaml` to customize:
- Default partition
- SLURM account/QoS
- Resource limits

**Usage:**
```bash
snakemake -s workflow/Snakefile --profile profiles/slurm
```

**Required SLURM modules:**
- `sbatch`, `squeue`, `sacct`, `scancel`

### `generic/`

Generic profile for local execution with resource management.

**Usage:**
```bash
snakemake -s workflow/Snakefile --profile profiles/generic
```

## Quick Reference

| Profile | Use Case | Job Submission |
|---------|----------|----------------|
| `slurm` | HPC clusters | sbatch |
| `generic` | Local workstation | local |

## Log Output

When using SLURM, logs are written to:
- `logs/slurm/{rule}-%j.out` (stdout)
- `logs/slurm/{rule}-%j.err` (stderr)

Where `%j` is the SLURM job ID.

## Customization

Create your own profile by copying an existing one:

```bash
cp -r profiles/slurm profiles/my_cluster
# Edit profiles/my_cluster/config.yaml
snakemake -s workflow/Snakefile --profile profiles/my_cluster
```
