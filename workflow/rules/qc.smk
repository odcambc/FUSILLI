rule multiqc_dir:
    """Final QC: aggregate FastQC and intermediate log files into a final report with MultiQC.
  """
    input:
        expand(
            "stats/{{experiment}}/{sample_prefix}_map.stats.txt",
            sample_prefix=samples,
        ),
        expand(
            "stats/{{experiment}}/{sample_prefix}_map.idxstats.txt",
            sample_prefix=samples,
        ),
        expand(
            "stats/{{experiment}}/{sample_prefix}_map.flagstats.txt",
            sample_prefix=samples,
        ),
        expand("stats/{{experiment}}/fastqc/{file}_R1_001_fastqc.html", file=files),
        expand("stats/{{experiment}}/fastqc/{file}_R2_001_fastqc.html", file=files),
    output:
        "stats/{experiment}/{experiment}_multiqc.html",
    benchmark:
        "benchmarks/{experiment}/multiqc.benchmark.txt"
    log:
        "logs/{experiment}/multiqc.log",
    wrapper:
        "v3.1.0/bio/multiqc"


rule fastqc:
    """Initial QC: run FastQC on all input reads."""
    input:
        expand("{data_dir}/{{sample_read}}.fastq.gz", data_dir=config["data_dir"]),
    output:
        html="stats/{experiment}/fastqc/{sample_read}_fastqc.html",
        zip="stats/{experiment}/fastqc/{sample_read}_fastqc.zip",
    params:
        "--quiet",
    benchmark:
        "benchmarks/{experiment}/{sample_read}.fastqc.benchmark.txt"
    log:
        "logs/{experiment}/fastqc/{sample_read}.log",
    threads: 8
    resources:
        mem_mb=config["mem_fastqc"],
    wrapper:
        "v3.1.0/bio/fastqc"
