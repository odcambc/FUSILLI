# FUSILLI: Fusion Utility for Scanning and Identification of Library Linked Interactions

**A Snakemake pipeline for fusion mutational scanning**

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for processing short-read sequencing data from libraries of
fusion constructs. The pipeline is designed to take in raw reads
and produce counts of detected fusion breakpoints.

## Quick start

```bash
git clone https://github.com/odcambc/FUSILLI
cd FUSILLI
conda env create --file fusilli_env.yaml
conda activate FUSILLI
```

Note: for ARM64 Macs, try the following (assuming [Rosetta](https://support.apple.com/en-us/102527) is installed):

```bash
git clone https://github.com/odcambc/FUSILLI
cd FUSILLI
CONDA_SUBDIR=osx-64 conda env create --file fusilli_env.yaml
conda activate FUSILLI
```

If the environment installed and activated properly,
edit the configuration files in the `config` directory as needed. Then run the pipeline with:

```bash
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 16
```

## License

This is licensed under the MIT license. See the LICENSE file for details.

## Contributing

Contributions and feedback are welcome. Please submit an issue or pull request.

## Getting help

For any issues, please open an issue on the GitHub repository. For
questions or feedback, [email Chris](https://www.wcoyotelab.com/members/).
