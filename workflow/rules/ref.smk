rule prepare_bbmap_index:
    """Generate the index for mapping with bbmap. This must be run once before mapping."""
    input:
        expand(
            "{ref_dir}/{reference}",
            ref_dir=config["ref_dir"],
            reference=config["reference"],
        ),
    output:
        "ref/genome/1/chr1.chrom.gz",
    threads: 16
    log:
        expand(
            "logs/{experiment}/bbmap/{reference}.bbmap_index.log",
            experiment=experiment,
            reference=config["reference"],
        ),
    shell:
        "bbmap.sh "
        "ref={input}"
