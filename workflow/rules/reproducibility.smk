"""
Reproducibility capture rules.

This module captures environment and invocation details for each pipeline run,
storing artifacts in the pipeline output directory.

Outputs:
- results/{experiment}/repro/metadata.json - Machine-readable metadata
- results/{experiment}/repro/metadata.txt - Human-readable metadata
- results/{experiment}/repro/conda-env.yaml - Conda environment export (if available)
- results/{experiment}/repro/pip-freeze.txt - Pip freeze output (if available)
"""


rule capture_reproducibility:
    """
    Capture reproducibility metadata for this pipeline run.

    This rule runs automatically as part of the pipeline and records:
    - Command-line invocation
    - Snakemake version
    - Python version
    - OS information
    - Dependency snapshots (conda env export and pip freeze)

    The outputs are stored in results/{experiment}/repro/ and are included
    in the default pipeline targets.
    """
    output:
        metadata_json="results/{experiment}/repro/metadata.json",
        metadata_txt="results/{experiment}/repro/metadata.txt",
        conda_env="results/{experiment}/repro/conda-env.yaml",
        pip_freeze="results/{experiment}/repro/pip-freeze.txt"
    params:
        output_dir=lambda wildcards: f"results/{wildcards.experiment}/repro"
    log:
        "logs/{experiment}/reproducibility.log"
    shell:
        """
        python workflow/scripts/capture_reproducibility.py {params.output_dir} 2>&1 | tee {log}
        """

