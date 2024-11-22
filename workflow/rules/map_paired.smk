
rule map_to_reference_minimap2:
    """Map reads to reference sequence using minimap2."""
    input:
        R1_ec="results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz",
        R2_ec="results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz",
    output:
        "results/{experiment}/{sample_prefix}.minimap2_mapped.sam",
    params:
        ref=config["reference"],
        ref_dir=config["ref_dir"],
    log:
        "logs/{experiment}/minimap2/{sample_prefix}.minimap2_map.log",
    threads: 32
    shell:
        "minimap2 -ax sr {params.ref_dir}/{params.ref} {input.R1_ec} {input.R2_ec} -o {output} 2> {log}"


rule sort_samtools:
    """Sort mapped reads using samtools."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_mapped.sam",
    output:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam",
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.samtools_sort.benchmark.txt"
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools.sort.log",
    wrapper:
        "v4.6.0/bio/samtools/sort"


rule index_samtools:
    """Index mapped reads using samtools."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam",
    output:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam.bai",
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.samtools_index.benchmark.txt"
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools.index.log",
    wrapper:
        "v4.6.0/bio/samtools/index"


rule stats_samtools:
    """Generate alignment stats with samtools."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam",
    output:
        "stats/{experiment}/{sample_prefix}_map.stats.txt",
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools_stats.log",
    wrapper:
        "v4.6.0/bio/samtools/stats"


rule flagstat_samtools:
    """Generate alignment stats with samtools."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam",
    output:
        "stats/{experiment}/{sample_prefix}_map.flagstats.txt",
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools_flagstats.log",
    wrapper:
        "v4.6.0/bio/samtools/flagstat"


rule idxstats_samtools:
    """Generate alignment stats with samtools."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam",
    output:
        "stats/{experiment}/{sample_prefix}_map.idxstats.txt",
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools_idxstats.log",
    wrapper:
        "v4.6.0/bio/samtools/idxstats"


rule count_sam_matches:
    """Count the number of mapped reads."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam",
    output:
        expand(
            "stats/{{experiment}}/{{sample_prefix}}_map.{reference_names}_count.txt",
            reference_names=reference_names,
        ),
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools_count.log",
    shell:
        "samtools view -F 0x800 {input} > {output}"
