rule samtools_filter_chimeric:
    """Filter output bams by presence of supplemtary alignment."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam",
    output:
        temp("results/{experiment}/{sample_prefix}.minimap2_chimeric.bam"),
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools_count.log",
    shell:
        "samtools view -d SA {input} > {output}"


rule process_chimeric:
    """Process chimeric reads."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_sorted.bam",
    output:
        "results/{experiment}/counts/{sample_prefix}.fusion_counts.csv",
    params:
        domain_list=domain_list,
        domain_dict=domain_dict,
    log:
        "logs/{experiment}/samtools/{sample_prefix}.fusion_counts.log",
    script:
        "scripts/parse_fusion_bam.py"
