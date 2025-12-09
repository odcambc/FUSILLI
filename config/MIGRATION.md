# Configuration Migration Guide

This guide helps migrate from the old configuration format to the new, cleaner system.

## Overview of Changes

| Old Format | New Format |
|------------|------------|
| `kinases_second.yaml` | `config.yaml` |
| `kinase_fusions.csv` (domain list) | `fusion_partners.csv` |
| `kinase_samples.csv` | `samples.csv` |
| `fusion_target` + `fusion_list` | `fusion_library` section |

## Step-by-Step Migration

### 1. Main Config File

**Old format (`kinases_second.yaml`):**
```yaml
experiment: 'kinase_fusions_second'
data_dir: '/path/to/data'
ref_dir: 'references'
experiment_file: 'config/kinase_samples_hov.csv'
fusion_list: 'config/kinase_fusions_hov.csv'
reference: 'kinase_sequences_hov.fasta'
paired_reads: true
fusion_target: 'Met_WT'
linker_sequence: GGGAGC
baseline_condition: 'baseline'
run_qc: false
kmers: 15
```

**New format (`config.yaml`):**
```yaml
experiment: 'kinase_fusions_second'

data_dir: '/path/to/data'
ref_dir: 'references'
samples_file: 'config/samples.csv'

fusion_library:
  anchor:
    name: 'Met_WT'
    position: 'downstream'
  linker_sequence: 'GGGAGC'
  partners_file: 'config/fusion_partners.csv'
  sequences_file: 'kinase_sequences_hov.fasta'

detection:
  method: 'string'
  breakpoint_window: 12
  maintain_frame: true
  kmer_size: 15

sequencing:
  paired: true
  min_quality: 30

qc:
  run_qc: false
  baseline_condition: 'baseline'
```

### 2. Fusion Partners File

**Old format (`kinase_fusions.csv`):**
```csv
domain,end_pos,full_length,upstream
Met_WT,1377,false,0
TPR,426,true,1
CCDC6,303,true,1
```

**New format (`fusion_partners.csv`):**
```csv
partner_name,sequence_length,include,description
TPR,426,true,TPR-MET fusion partner
CCDC6,303,true,RET fusion partner
```

**Key changes:**
- Renamed `domain` → `partner_name`
- Renamed `end_pos` → `sequence_length` (clearer meaning)
- Renamed `full_length` → `include`
- Removed `upstream` column (now in main config as `anchor.position`)
- Anchor domain (Met_WT) is no longer listed here - it's defined in config
- Added optional `description` column

### 3. Samples File

The samples file format is largely unchanged, but comments are now supported:

**New format supports comments:**
```csv
# Sample definitions for kinase fusion experiment
# Lines starting with # are ignored
sample,condition,replicate,time,tile,file
gDNA_fusions,baseline,1,0,1,gDNA_S241_L004
pDNA_fusions,baseline,2,0,1,pDNA_S242_L004
```

## Conversion Script

If you have many domains to convert, here's a quick Python script:

```python
import csv

# Read old format
with open('old_kinase_fusions.csv', 'r') as f:
    reader = csv.DictReader(f)
    old_rows = list(reader)

# Write new format
with open('fusion_partners.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f,
        fieldnames=['partner_name', 'sequence_length', 'include', 'description'])
    writer.writeheader()

    for row in old_rows:
        # Skip the anchor domain
        if row['upstream'] == '0':
            continue

        writer.writerow({
            'partner_name': row['domain'],
            'sequence_length': row['end_pos'],
            'include': 'true' if row['full_length'].lower() == 'true' else 'false',
            'description': ''
        })
```

## Validating Your Migration

After migrating, run a dry-run to check for errors:

```bash
cd fusilli
snakemake -s workflow/Snakefile -n
```

If successful, you'll see the pipeline plan. If there are config errors, you'll get
helpful error messages indicating what needs to be fixed.

