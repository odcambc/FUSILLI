rule filter_awk:
    """Filter output of cutadapt to find matching sequences.
    This rule finds reads where the barcode was detected and outputs
    the upstream sequence along with the identified barcode."""
    input:
        "results/cutadapt/{sample}.barcodes.info.tsv",
    output:
        "results/cutadapt/{sample}.barcodes.matches.txt",
    params:
        max_barcode_length=40,
    shell:
        "awk '($2 !~ /-1/ && $8 ~ /1;2/ && length($5) < {params.max_barcode_length}) {{print $5}}' {input} > {output}"


rule make_bam:
    """Convert sam to bam using samtools."""
    input:
        "results/clusters/{sample}/mapped/{barcode}.sam",
    output:
        temp("results/clusters/{sample}/mapped/{barcode}.bam"),
    log:
        "logging/consensus/{sample}/{barcode}_bam.log",
    shell:
        "samtools sort {input} > {output} 2> {log}"


rule index_bam:
    """Index bam files using samtools."""
    input:
        "results/clusters/{sample}/mapped/{barcode}.bam",
    output:
        temp("results/clusters/{sample}/mapped/{barcode}.bam.bai"),
    log:
        "logging/consensus/{sample}/{barcode}_bam.log",
    shell:
        "samtools index {input} 2> {log}"


rule mpileup:
    """Call variants using bcftools mpileup."""
    input:
        bam="results/consensus/{sample}/{barcode}_consensus.bam",
        reference=reference_file,
        indexed="results/consensus/{sample}/{barcode}_consensus.bam.bai",
    output:
        temp("results/consensus/{sample}/{barcode}_consensus.bcf"),
    log:
        "logging/consensus/{sample}/{barcode}_mpileup.log",
    shell:
        "bcftools mpileup -d 5000 -Ou -f {input.reference} {input.bam} | bcftools call -vm -Ob --ploidy 1 -o {output} 2> {log}"
