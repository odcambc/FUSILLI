def get_file_from_sample(wildcards):
    """Maps from the sequencing output file names to the sample names defined in the experiment CSV.
    This is used in rule bbduk_trim_adapters"""
    read_number = wildcards.sample_prefix[-2:]
    sample = wildcards.sample_prefix[:-3]

    if read_number == "R1" or read_number == "R2":
        filename = experiments.loc[experiments["sample"] == sample, "file"].squeeze()
        prefix = config["data_dir"] + "/" + filename
        full_file = prefix + "_" + read_number + "_001.fastq.gz"
        return full_file
    else:
        warnings.warn(
            "The read number is not R1 or R2. Please check the naming convention."
        )


def get_ref(wildcards):
    """Removes file suffix from reference fasta file. This is used in rule bwa_index"""
    prefix = config["reference"].split(".fasta")[0]
    return prefix


def get_ref_names(ref_file):
    """Returns a list of reference names from a fasta file. This is used to extract
    the potentially fused sequences."""
    ref_names = []
    for record in SeqIO.parse(ref_file, "fasta"):
        ref_names.append(record.id)
    return ref_names


def pass_names(names):
    """Takes either an array of strings (such as files) or a single string and returns a comma-separated string.
    This is used while passing references to bbduk."""
    if isinstance(names, str):
        return names
    else:
        return ",".join(names)


def get_baseline_samples(experiments, samples):
    """Returns a list of baseline sample and files"""
    baseline_samples = []
    baseline_files = []
    for sample in samples:
        if experiments.loc[sample, "condition"] == config["baseline_condition"]:
            baseline_samples.append(sample)
            baseline_files.append(experiments.loc[sample, "file"])
    return baseline_samples, baseline_files


def get_experiment_samples(experiments, samples):
    """Returns a list of experiment (i.e., non-baseline) sample and files"""
    experiment_samples = []
    experiment_files = []
    for sample in samples:
        if experiments.loc[sample, "condition"] != config["baseline_condition"]:
            experiment_samples.append(sample)
            experiment_files.append(experiments.loc[sample, "file"])
    return experiment_samples, experiment_files


def get_input(wildcards):
    """Generate the input files for the dummy rule all.
    This is necessary to allow optional pipeline outputs."""

    input_list = []

    if experiment_samples:
        if paired:
            input_list.extend(
                expand(
                    "results/{experiment_name}/{sample_name}.fusion_counts.csv",
                    experiment_name=config["experiment"],
                    sample_name=samples,
                )
            )
        else:
            input_list.extend(
                expand(
                    "results/{experiment_name}/counts/{sample_name}.R1.fusion_counts.csv",
                    experiment_name=config["experiment"],
                    sample_name=samples,
                )
            )
            input_list.extend(
                expand(
                    "results/{experiment_name}/counts/{sample_name}.R2.fusion_counts.csv",
                    experiment_name=config["experiment"],
                    sample_name=samples,
                )
            )
        if config["run_qc"]:
            input_list.extend(
                expand(
                    "stats/{experiment_name}/{experiment_name}_multiqc.html",
                    experiment_name=config["experiment"],
                )
            )
    return input_list


# Validate config and experiment files
validate(config, "../schemas/config.schema.yaml")
# validate(config["experiment_file"], "../schemas/experiments.schema.yaml")
# validate(config["fusion_list"], "../schemas/domains.schema.yaml")


# Set variables from config and experiment files
experiments = (
    pd.read_csv(config["experiment_file"], header=0)
    .dropna(how="all")
    .set_index("sample", drop=False, verify_integrity=True)
)

# Set up sample variables
experiment = config["experiment"]
samples = experiments["sample"]
files = experiments["file"]

domain_dict = {}
domain_list = []

with open(config["fusion_list"], "r") as f:
    lines = f.readlines()[1:]

domain_dict = {line.split(",")[0]: int(line.split(",")[1]) for line in lines}
domain_list = list(domain_dict.keys())

experiment_samples, experiment_files = get_experiment_samples(experiments, samples)
baseline_samples, baseline_files = get_baseline_samples(experiments, samples)

conditions = set(experiments["condition"])
experimental_conditions = conditions - set([config["baseline_condition"]])

# Load experiment config and domain list files


fusion_list = (
    pd.read_csv(config["fusion_list"], header=0).dropna(how="all")
).set_index("domain", drop=False, verify_integrity=True)

# Set reference file variables
reference_names = get_ref_names(config["ref_dir"] + "/" + config["reference"])
reference_name = get_ref(config["reference"])

# Set filtering and mapping variables
adapters_ref = pass_names(config["adapters"])
contaminants_ref = pass_names(config["contaminants"])
samtools_local = config["samtools_local"]

# Set experimental design variables
paired = config["paired_reads"]

tiles = experiments["tile"].unique()
if len(tiles) > 1:
    config["tiled"] = True
else:
    config["tiled"] = False
