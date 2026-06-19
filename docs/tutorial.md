# Tutorial: your first FUSILLI run

This walkthrough takes you from a fresh clone to interpreting fusion counts, using the
small synthetic dataset bundled in `tests/`. No real sequencing data is required. It
takes a few minutes and doubles as a check that your installation works.

You will:

1. Install and activate the environment
2. Run the pipeline on the bundled test dataset
3. Find and read the output
4. Confirm you got the expected result

## The test library

The fixture mimics a tiny kinase-fusion library:

- **Retained domain:** `Met_WT` (240 nt) — the constant domain present in every fusion
- **Fusion partners:** `TPR` (426 nt) and `CCDC6` (303 nt), each fused to the retained domain through a `GGGAGC` linker
- **Reads:** 170 synthetic paired-end read pairs in `tests/data/`

The reads were built to contain specific fusion junctions, so we know what the pipeline
should find — the expected counts are recorded in `tests/expected_counts.csv`.

## 1. Install and activate

```bash
git clone https://github.com/odcambc/FUSILLI
cd FUSILLI
conda env create --file fusilli_env.yaml
conda activate FUSILLI
```

This environment provides Snakemake, BBTools (`bbduk.sh`, `bbmerge.sh`), and the
fast-matching libraries the detector uses.

## 2. Dry run, then run

Everything is preconfigured in `config/test.yaml`, which points at the `tests/` fixtures.
First preview the plan without executing anything:

```bash
snakemake -s workflow/Snakefile counts_only --configfile config/test.yaml -n
```

Then run the counting pipeline:

```bash
snakemake -s workflow/Snakefile counts_only --configfile config/test.yaml --cores 4
```

(The target — `counts_only` — comes **before** `--configfile`, because `--configfile`
accepts multiple files and would otherwise swallow the target name.)

The `counts_only` target runs the full read-processing chain (adapter trim →
decontaminate → quality filter → merge) followed by the two-stage fusion matcher, then
aggregates the counts. It skips the optional FastQC/MultiQC reports. (Omitting
`counts_only` runs the default target, which additionally writes reproducibility
metadata under `results/test_experiment/repro/`.)

## 3. Find the output

Results land under `results/test_experiment/`:

```
results/test_experiment/
├── references/
│   ├── junction_sequences.csv     # candidate junction k-mers generated from the references
│   └── domain_ends.csv            # the fusion partners' 3′ ends used for the fast pre-filter
├── counts/
│   └── test_sample.fusion_counts.csv   # per-sample counts  ← the main result
├── fusion_counts_summary.csv      # counts as a matrix across samples
├── fusion_qc_metrics.csv          # detection efficiency, diversity, coverage
├── sensitivity_metrics.csv
└── decay_metrics.csv              # read retention through each processing step
```

Look at the per-sample counts:

```bash
column -s, -t results/test_experiment/counts/test_sample.fusion_counts.csv | head
```

## 4. Confirm the expected result

Five fusions are detected (all other candidates are reported with a count of 0):

| fusion_id        | count |
| ---------------- | ----- |
| TPR_426_Met_WT   | 30    |
| TPR_423_Met_WT   | 25    |
| TPR_420_Met_WT   | 20    |
| CCDC6_303_Met_WT | 15    |
| CCDC6_300_Met_WT | 10    |

These are `TPR` and `CCDC6` fusion partners, truncated at several breakpoints, fused to the
retained `Met_WT` domain — exactly the counts in `tests/expected_counts.csv`. If you see
them, the pipeline is working.

## What a fusion_id means

`TPR_426_Met_WT` reads as **fusion partner `TPR`, breakpoint coordinate 426, retained
domain `Met_WT`**. The test config sets `truncated_component: 'partner'`, so the
breakpoint coordinate is measured from the fusion partner's start — `426` is the full length of
TPR, i.e. a fusion that retains all of TPR. See
[Breakpoint Conventions](../README.md#breakpoint-conventions) for the full explanation,
including how each junction k-mer is built (here, `2 × breakpoint_window` = 24 nt).

## Next steps

- **Use your own data.** Copy `config/examples/config.yaml`, point `data_dir`,
  `samples_file`, and the reference/partner files at your own inputs, and run without
  `--configfile config/test.yaml`.
- **Look up any option.** `config/examples/config.yaml` is the annotated template;
  `workflow/schemas/config.schema.yaml` is the authoritative schema (types, defaults,
  allowed values) the pipeline validates against at startup.
