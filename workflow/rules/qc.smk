rule multiqc_dir:
    """Final QC: aggregate FastQC and intermediate log files into a final report with MultiQC."""
    input:
        expand("stats/{{experiment}}/fastqc/{sample}_{read}.fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("stats/{{experiment}}/merge/{sample}.ihist", sample=SAMPLES),
    output:
        "stats/{experiment}/{experiment}_multiqc.html",
    benchmark:
        "benchmarks/{experiment}/multiqc.benchmark.txt"
    log:
        "logs/{experiment}/multiqc.log",
    conda:
        "../envs/qc.yaml",
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
