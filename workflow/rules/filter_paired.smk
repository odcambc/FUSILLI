"""
Paired-end read preprocessing rules.

This module contains rules for:
1. Adapter trimming
2. Contaminant removal
3. Error correction and merging
"""


rule trim_adapters:
    """
    Remove sequencing adapters using BBDuk.

    Uses the standard Illumina adapter sequences from BBTools.
    Also generates QC histograms for downstream analysis.
    """
    input:
        R1=get_raw_fastq_r1,
        R2=get_raw_fastq_r2
    output:
        R1_trim=temp("results/{experiment}/trimmed/{sample}_R1.trim.fastq.gz"),
        R2_trim=temp("results/{experiment}/trimmed/{sample}_R2.trim.fastq.gz"),
        qhist="stats/{experiment}/trim/{sample}.qhist",
        bhist="stats/{experiment}/trim/{sample}.bhist",
        gchist="stats/{experiment}/trim/{sample}.gchist",
        aqhist="stats/{experiment}/trim/{sample}.aqhist",
        lhist="stats/{experiment}/trim/{sample}.lhist",
        stats="stats/{experiment}/trim/{sample}.stats.txt"
    params:
        adapters=ADAPTERS_REF
    log:
        "logs/{experiment}/bbduk/{sample}.trim.log"
    benchmark:
        "benchmarks/{experiment}/{sample}.bbduk_trim.txt"
    threads: DEFAULT_THREADS
    resources:
        mem_mb=DEFAULT_MEMORY
    shell:
        """
        bbduk.sh \
            in1={input.R1} \
            in2={input.R2} \
            out1={output.R1_trim} \
            out2={output.R2_trim} \
            ref={params.adapters} \
            ktrim=r k=23 mink=10 hdist=1 tpe tbo \
            bhist={output.bhist} \
            qhist={output.qhist} \
            gchist={output.gchist} \
            aqhist={output.aqhist} \
            lhist={output.lhist} \
            stats={output.stats} \
            overwrite=true \
            t={threads} \
            gcbins=auto \
            2> {log}
        """


rule remove_contaminants:
    """
    Remove contaminants (PhiX, sequencing artifacts) using BBDuk.
    """
    input:
        R1_trim="results/{experiment}/trimmed/{sample}_R1.trim.fastq.gz",
        R2_trim="results/{experiment}/trimmed/{sample}_R2.trim.fastq.gz"
    output:
        R1_clean=temp("results/{experiment}/cleaned/{sample}_R1.clean.fastq.gz"),
        R2_clean=temp("results/{experiment}/cleaned/{sample}_R2.clean.fastq.gz"),
        stats="stats/{experiment}/contam/{sample}.stats.txt"
    params:
        contaminants=CONTAMINANTS_REF
    log:
        "logs/{experiment}/bbduk/{sample}.clean.log"
    benchmark:
        "benchmarks/{experiment}/{sample}.bbduk_clean.txt"
    threads: DEFAULT_THREADS
    resources:
        mem_mb=DEFAULT_MEMORY
    shell:
        """
        bbduk.sh \
            in={input.R1_trim} \
            in2={input.R2_trim} \
            out={output.R1_clean} \
            out2={output.R2_clean} \
            stats={output.stats} \
            ref={params.contaminants} \
            k=31 \
            overwrite=true \
            t={threads} \
            2> {log}
        """


rule merge_reads:
    """
    Merge paired-end reads with error correction using BBMerge.

    This step:
    1. Performs error correction using paired-end overlap
    2. Merges overlapping reads into single sequences

    The merged reads are used for fusion detection as they provide
    higher confidence sequence spanning the breakpoint.
    """
    input:
        R1_clean="results/{experiment}/cleaned/{sample}_R1.clean.fastq.gz",
        R2_clean="results/{experiment}/cleaned/{sample}_R2.clean.fastq.gz"
    output:
        merged="results/{experiment}/merged/{sample}_merged.fastq.gz",
        ihist="stats/{experiment}/merge/{sample}.ihist"
    log:
        "logs/{experiment}/bbmerge/{sample}.log"
    benchmark:
        "benchmarks/{experiment}/{sample}.bbmerge.txt"
    threads: DEFAULT_THREADS
    resources:
        mem_mb=DEFAULT_MEMORY
    shell:
        """
        bbmerge.sh \
            in={input.R1_clean} \
            in2={input.R2_clean} \
            out={output.merged} \
            ihist={output.ihist} \
            ecco mix \
            showhiststats=t \
            overwrite=true \
            t={threads} \
            2> {log}
        """
