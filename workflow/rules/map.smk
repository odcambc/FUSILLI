rule map_to_reference_bbmap:
    """Map reads to reference sequence using BBMap. covhist output is currently disabled as it causes MultiQC bloat."""
    input:
        "ref/genome/1/chr1.chrom.gz",
        R1_ec="results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz",
        R2_ec="results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz",
    output:
        mapped=temp("results/{experiment}/{sample_prefix}.mapped.sam"),
        covstats="stats/{experiment}/{sample_prefix}_map.covstats",
        basecov="stats/{experiment}/{sample_prefix}_map.basecov",
        bincov="stats/{experiment}/{sample_prefix}_map.bincov",
        ehist="stats/{experiment}/{sample_prefix}_map.ehist",
        indelhist="stats/{experiment}/{sample_prefix}_map.indelhist",
        mhist="stats/{experiment}/{sample_prefix}_map.mhist",
        idhist="stats/{experiment}/{sample_prefix}_map.idhist",
    params:
        ref=config["reference"],
        ref_dir=config["ref_dir"],
        sam=config["sam"],
        kmers=config["kmers"],
        mem=config["mem"],
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.bbmap_map.benchmark.txt"
    log:
        "logs/{experiment}/bbmap/{sample_prefix}.bbmap_map.log",
    threads: 16
    shell:
        "bbmap.sh "
        "-Xmx{params.mem}g "
        "in1={input.R1_ec} "
        "in2={input.R2_ec} "
        "sam={params.sam} 32bit=t "
        "ref={params.ref_dir}/{params.ref} "
        "outm={output.mapped} "
        "k={params.kmers} "
        "t={threads} "
        "nzo=true "
        "covstats={output.covstats} "
        "basecov={output.basecov} "
        "bincov={output.bincov} "
        "ehist={output.ehist} "
        "indelhist={output.indelhist} "
        "mhist={output.mhist} "
        "idhist={output.idhist} 2> {log}"


rule map_to_reference_bbmap_minimap2:
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


rule sort_index_samtools:
    """Sort and index mapped reads using samtools."""
    input:
        "results/{experiment}/{sample_prefix}.minimap2_mapped.sam",
    output:
        bam="results/{experiment}/{sample_prefix}.minimap2_mapped.bam",
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.samtools_sort.benchmark.txt"
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools.log",
    params:
        extra="-O bam --write-index -o {output.bam}",
    wrapper:
        "v4.6.0/bio/samtools/sort"


rule stats_samtools:
    """Generate alignment stats with samtools."""
    input:
        bam="results/{experiment}/{sample_prefix}.minimap2_mapped.bam",
    output:
        "stats/{experiment}/{sample_prefix}_map.stats.txt",
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools_stats.log",
    wrapper:
        "v4.6.0/bio/samtools/stats"


rule flagstat_samtools:
    """Generate alignment stats with samtools."""
    input:
        bam="results/{experiment}/{sample_prefix}.minimap2_mapped.bam",
    output:
        "stats/{experiment}/{sample_prefix}_map.flagstats.txt",
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools_flagstats.log",
    wrapper:
        "v4.6.0/bio/samtools/flagstat"


rule idxstats_samtools:
    """Generate alignment stats with samtools."""
    input:
        bam="results/{experiment}/{sample_prefix}.minimap2_mapped.bam",
    output:
        "stats/{experiment}/{sample_prefix}_map.idxstats.txt",
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools_idxstats.log",
    wrapper:
        "v4.6.0/bio/samtools/idxstats"
