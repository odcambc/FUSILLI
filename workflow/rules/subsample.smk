"""
Quick-mode subsampling for fast sanity check runs.
"""

rule subsample_reads:
    """
    Subsample paired FASTQs for quick mode using bbmap reformat.sh.
    """
    input:
        R1=get_raw_fastq_r1,
        R2=get_raw_fastq_r2,
    output:
        R1_sub=temp("results/{experiment}/quick/{sample}_R1.sub.fastq.gz"),
        R2_sub=temp("results/{experiment}/quick/{sample}_R2.sub.fastq.gz"),
    params:
        max_reads=lambda wildcards: QUICK_MAX_READS,
        fraction=lambda wildcards: "" if QUICK_FRACTION in (None, "None") else QUICK_FRACTION,
        seed=lambda wildcards: QUICK_SEED,
    log:
        "logs/{experiment}/subsample/{sample}.log",
    benchmark:
        "benchmarks/{experiment}/{sample}.subsample.txt"
    threads: DEFAULT_THREADS
    resources:
        mem_mb=DEFAULT_MEMORY
    shell:
        r"""
        samplerate_flag=""
        limit_flag=""
        if [ -n "{params.fraction}" ]; then
            samplerate_flag="samplerate={params.fraction} sampleseed={params.seed}"
        else
            # Take the first N reads for speed (no full-file sampling)
            limit_flag="reads={params.max_reads}"
        fi

        reformat.sh \
            in={input.R1:q} \
            in2={input.R2:q} \
            out={output.R1_sub:q} \
            out2={output.R2_sub:q} \
            $samplerate_flag \
            $limit_flag \
            overwrite=true \
            t={threads} \
            2> {log:q}
        """
