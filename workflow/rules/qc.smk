def get_multiqc_inputs(wildcards):
    inputs = []
    inputs.extend(
        expand(
            "stats/{{experiment}}/fastqc/{sample}_{read}.fastqc.zip",
            sample=SAMPLES,
            read=["R1", "R2"],
        )
    )
    inputs.append("stats/{experiment}/fastqc")
    inputs.extend(expand("stats/{{experiment}}/merge/{sample}.ihist", sample=SAMPLES))
    inputs.extend(expand("stats/{{experiment}}/trim/{sample}.trim.stats.txt", sample=SAMPLES))
    inputs.extend(expand("stats/{{experiment}}/contam/{sample}.contam.stats.txt", sample=SAMPLES))
    inputs.extend(
        expand(
            "stats/{{experiment}}/fastqc/{sample}_{read}.fastqc.html",
            sample=SAMPLES,
            read=["R1", "R2"],
        )
    )
    inputs.extend(expand("logs/{{experiment}}/bbmerge/{sample}.log", sample=SAMPLES))
    inputs.extend(expand("logs/{{experiment}}/bbduk/{sample}.trim.log", sample=SAMPLES))
    inputs.extend(expand("logs/{{experiment}}/bbduk/{sample}.clean.log", sample=SAMPLES))
    inputs.extend(expand("logs/{{experiment}}/bbduk/{sample}.quality.log", sample=SAMPLES))
    inputs.extend([
        "results/{experiment}/fusion_qc_metrics.csv",
        "results/{experiment}/sensitivity_metrics.csv",
        "results/{experiment}/decay_metrics.csv",
        "results/{experiment}/trim_metrics.csv",
        "results/{experiment}/contam_metrics.csv",
        "results/{experiment}/quality_metrics.csv",
        "results/{experiment}/partner_counts_summary.csv",
        "results/{experiment}/fusion_counts_summary.csv",
    ])
    inputs.extend(
        expand(
            "results/{{experiment}}/counts/{sample}.fusion_metrics.json",
            sample=SAMPLES
        )
    )
    inputs.extend(
        expand(
            "results/{{experiment}}/counts/{sample}.partner_counts.csv",
            sample=SAMPLES
        )
    )
    if UNMERGED_DETECTION:
        inputs.extend([
            "results/{experiment}/unmerged_qc_metrics.csv",
            "results/{experiment}/unmerged_counts_summary.csv",
            "results/{experiment}/unmerged_partner_counts_summary.csv",
        ])
        inputs.extend(
            expand(
                "results/{{experiment}}/counts/{sample}.{mate}.unmerged_fusion_metrics.json",
                sample=SAMPLES,
                mate=["R1", "R2"]
            )
        )
        inputs.extend(
            expand(
                "results/{{experiment}}/counts/{sample}.{mate}.unmerged_partner_counts.csv",
                sample=SAMPLES,
                mate=["R1", "R2"]
            )
        )
    return inputs


rule multiqc_dir:
    """Final QC: aggregate FastQC and intermediate log files into a final report with MultiQC."""
    input:
        get_multiqc_inputs
    output:
        "stats/{experiment}/{experiment}_multiqc.html",
    benchmark:
        "benchmarks/{experiment}/multiqc.benchmark.txt"
    log:
        "logs/{experiment}/multiqc.log",
    conda:
        "../envs/qc.yaml",
    params:
        config="config/multiqc_config.yaml"
    wrapper:
        "v3.1.0/bio/multiqc"


rule fastqc:
    """Initial QC: run FastQC on all input reads (raw or subsampled)."""
    input:
        fq=lambda wildcards: get_input_fastq_r1(wildcards) if wildcards.read == "R1" else get_input_fastq_r2(wildcards),
    output:
        html="stats/{experiment}/fastqc/{sample}_{read}.fastqc.html",
        zip="stats/{experiment}/fastqc/{sample}_{read}.fastqc.zip",
    params:
        "--quiet",
    benchmark:
        "benchmarks/{experiment}/{sample}_{read}.fastqc.benchmark.txt"
    log:
        "logs/{experiment}/fastqc/{sample}_{read}.log",
    threads: 8
    resources:
        mem_mb=config["mem_fastqc"],
    conda:
        "../envs/qc.yaml",
    wrapper:
        "v3.1.0/bio/fastqc"
