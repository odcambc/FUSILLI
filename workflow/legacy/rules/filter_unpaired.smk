rule bbduk_trim_adapters:
    """Remove adapters with BBDuk. This uses all of the adapters in the adapters.fa file, which is included with BBTools.
    This also performs the first sample renaming step, using the mapping file -> sample from the provided experiment CSV."""
    input:
        get_file_from_sample,
    output:
        trim=temp("results/{experiment}/{sample_prefix}.trim.fastq.gz"),
        qhist="stats/{experiment}/{sample_prefix}_trim.qhist",
        bhist="stats/{experiment}/{sample_prefix}_trim.bhist",
        gchist="stats/{experiment}/{sample_prefix}_trim.gchist",
        aqhist="stats/{experiment}/{sample_prefix}_trim.aqhist",
        lhist="stats/{experiment}/{sample_prefix}_trim.lhist",
        stats="stats/{experiment}/{sample_prefix}_trim.stats.txt",
    params:
        adapters=adapters_ref,
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.bbduk_trim.benchmark.txt"
    log:
        "logs/{experiment}/bbduk/{sample_prefix}.trim.bbduk.log",
    threads: 16
    shell:
        "bbduk.sh in={input} "
        "ref={params.adapters} ktrim=r k=23 mink=10 hdist=1 tpe tbo "
        "out={output.trim} "
        "bhist={output.bhist} "
        "qhist={output.qhist} "
        "gchist={output.gchist} "
        "aqhist={output.aqhist} "
        "lhist={output.lhist} "
        "stats={output.stats} "
        "overwrite=true "
        "t={threads} "
        "gcbins=auto 2> {log}"


rule remove_contaminants:
    """Remove typical contaminants with BBDuk.
    The contaminants files used, by default, are included with BBTools, and consist of PhiX and some (possibly JGI-specific) other sequencing artifacts."""
    input:
        "stats/{experiment}/{sample_prefix}_trim.stats.txt",
        trim="results/{experiment}/{sample_prefix}.trim.fastq.gz",
    output:
        clean=temp("results/{experiment}/{sample_prefix}.clean.fastq.gz"),
        stats="stats/{experiment}/{sample_prefix}_trim_contam.stats.txt",
    params:
        contaminants=contaminants_ref,
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.bbduk_clean.benchmark.txt"
    log:
        "logs/{experiment}/bbduk/{sample_prefix}.clean.bbduk.log",
    threads: 16
    shell:
        "bbduk.sh "
        "in={input.trim} "
        "out={output.clean} "
        "stats={output.stats} "
        "overwrite=true "
        "k=31 "
        "t={threads} "
        "ref={params.contaminants} 2> {log}"
