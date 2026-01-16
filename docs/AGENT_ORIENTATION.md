# Agent Orientation Guide

This document provides comprehensive orientation for AI agents and human editors working on the FUSILLI project. It explains both the biological context and technical implementation to ensure effective collaboration.

---

## Quick Reference Summary

**Technical Purpose:** FUSILLI is a Snakemake pipeline that takes short-read sequencing data and provides the per-sample counts of particular variants that are identified in the library. The variants are expected to be fusions of different genes occuring at different positions (truncation lengths) of one gene. The library is expected to be constructed by pooling variants together and sequencing the pooled library.

The pipeline performs read pre-processing with bbtools, preparation of the sequences of expected variants, variant detection by string matching, and counting. It also performs quality control and reproducibility capture.

Because sequencing reads are short (Illumina, 50-300 bp), not all reads will cover the breakpoint and allow identification of the variant. Reads that only cover either the partner domain or the anchor domain cannot be identified. Therefore we filter reads based on which can plausibly contain a breakpoint by first searching for domain end k-mers and then matching against the full set of breakpoints.

**Biological Purpose:** FUSILLI quantifies variant-specific fusion breakpoints in deep mutational scanning (DMS) experiments on fusion protein libraries (e.g., kinase domain fusions like MET), enabling researchers to measure the abundance of different fusion variants in their experimental samples. In these experiments, taking timepoints of a pooled library and measuring the change in abundance of different variants is used to identify the fitness effects of different variants under some condition.

---

## Biological Context

### Experimental Design

FUSILLI is designed for analyzing mutational scanning data from fusion protein libraries. The tool was created to support analysis of libraries consisting of various domains (e.g., TPR, CCDC6) fused to a variably truncated anchor domain (e.g., MET kinase).

**Key Biological Concepts:**

- **Fusion Proteins:** Constructs where a partner domain is fused to an anchor domain (typically a kinase domain)
- **Deep Mutational Scanning (DMS):** Experimental approach that systematically tests many variants to understand sequence-function relationships
- **Breakpoints:** The junction points between partner and anchor domains, representing different truncation positions
- **Variant Abundance:** The count of each fusion variant in a sample, which reflects its representation in the library through relative abundance.

### Fusion Library Structure

A typical fusion construct has this structure:

```(markdown)
[Partner N-term]────[Linker]────[Kinase Domain (truncated)]
      ↑                              ↑
Full-length partner            Variable truncation
(breakpoint)
```

**Example:** TPR-MET fusion

- Partner: TPR domain (amino acids 1-142, full-length)
- Linker: GGGAGC (or other linker sequence)
- Anchor: MET kinase domain (variably N-terminally truncated)

The "variable truncation" means different fusion variants have the anchor domain truncated at different positions, creating a library of fusion variants where each variant contains a different length of an N-terminally truncated anchor domain with some full-length partner domain.

In addition to the fusion variants, the library may also contain variants that are not fusions, such as control fusions with the anchor domain or variants that are not fusions at all. These will also be counted and reported as non-fusion variants with a distinct identification logic.

### Research Context

This tool builds on work studying MET kinase domain fusions and exon skipping in disrupting signaling pathways using DMS (see references in [README.md](README.md)). The pipeline processes sequencing data from libraries where:

- Multiple fusion variants are pooled together
- Each variant has a unique breakpoint position
- Sequencing reads contain information about which variant is present
- The goal is to count how many times each variant appears in each sample

---

## Technical Overview

### Pipeline Architecture

FUSILLI is a Snakemake-based workflow pipeline. See [ARCHITECTURE.md](ARCHITECTURE.md) for detailed system design.

**Main Components:**

1. **Preprocessing** (`filter_paired.smk`): Adapter trimming, quality filtering, read merging
2. **Reference Generation** (`process_strings.smk`): Creates breakpoint k-mer sequences
3. **Detection** (`process_strings.smk`): String matching to identify fusions in reads
4. **Counting** (`process_strings.smk`): Aggregates matches per fusion variant across all samples
5. **Quality Control** (`qc.smk`): FastQC, MultiQC reports
6. **Reproducibility** (`reproducibility.smk`): Captures run metadata

### Data Flow

```(markdown)
FASTQ files (paired-end)
    ↓
Preprocessing (trimming, merging)
    ↓
Reference generation (breakpoint k-mers)
    ↓
Detection (string matching)
    ↓
Counting (per-sample fusion counts)
    ↓
Summary (aggregated counts matrix)
```

See [ARCHITECTURE.md](ARCHITECTURE.md) for the complete data flow diagram.

### Key Technical Concepts

**Breakpoint K-mers:**
- Short sequences spanning the fusion junction
- Uniquely identifies each fusion variant
- Include nucleotides from both partner and anchor domains
- Include linker sequence if present
- Window size configurable (default: 12 nt each side of breakpoint)

**Pre-filtering (string matching):**
- Uses domain end k-mers to quickly filter potentially informative reads
- Only reads that contain a domain end k-mer can contain a breakpoint for identification
- Reduces search space before exact breakpoint matching
- Improves performance for large datasets

**Fusion ID Format:**
- `{partner}_{breakpoint_nt}_{anchor}`
- Example: `TPR_126_Met_WT`
  - Partner: TPR
  - Breakpoint: 126 nucleotides from anchor start ( .e.g, amino acids 1-42 of MET missing)
  - Anchor: Met_WT

---

## Key Terminology

| Term | Definition |
|------|------------|
| **Breakpoint** | The junction point between partner and anchor domains in a fusion protein, represented as a nucleotide position |
| **Fusion Variant** | A specific fusion construct with a particular breakpoint position |
| **Partner Domain** | The full-length domain in the fusion (e.g., TPR, CCDC6) |
| **Anchor Domain** | The domain that can be N-terminally truncated at different positions during fusion (e.g., MET kinase domain) |
| **Linker** | Short sequence connecting partner and anchor domains |
| **Breakpoint K-mer** | Short sequence spanning the fusion junction, used for detection |
| **Domain End** | The 3' end of the partner domain, used for pre-filtering reads |
| **Fusion ID** | Unique identifier encoding partner, breakpoint position, and anchor |
| **DMS** | Deep Mutational Scanning - experimental approach for systematic variant testing |

---

## Common Questions & Clarifications

**Note:** The following questions are placeholders for the project maintainer to fill in with specific details about their experimental setup and use cases. Edit this section to provide clarity for agents and editors.

### Biological Context Questions

**Q1: What is the exact experimental setup for the fusion library?**
- [ ] **Answer:** [The library is constructed by pooling variants together and sequencing the pooled library. Experiments are performed by taking timepoints of the pooled library under some experimental condition and using sequencing to measure the change in abundance of different variants over time.]

**Q2: How are the fusion variants created experimentally?**
- [ ] **Answer:** [The library is constructed through a Golden-Gate assembly with a designed oligo pool containing designed variants.]

**Q3: What does "variable truncation" mean in the experimental context?**
- [ ] **Answer:** [Variable truncation means that the partner domain is present fused with the anchor domain at different positions within the anchor domain. The truncation length is the number of nucleotides from the start of the anchor domain to the breakpoint that are deleted.]

**Q4: What is the relationship between nucleotide breakpoints and protein-level effects?**
- [ ] **Answer:** [This is dependent on the experimental system and context. In general, whatever variant effects occur are expected to result in changes of relative abundance over time which is interpreted as a fitness effect.]

**Q5: What specific research question does counting fusion variants answer?**
- [ ] **Answer:** [This varies by experiment and context and can not be answered in general.]

### Technical Questions

**Q6: Why is pre-filtering using domain ends necessary?**
- [ ] **Answer:** [This is necessary for performance and to avoid unnecessary computation. The breakpoint region which allows identification of the fusion variant is generally a small portion of the overall fusion sequence, and thereform random sequencing (as with tagmentation) will result in most reads not containing the breakpoint region. Tiled amplicon sequencing may not have this problem but is not yet explicitly supported.]

**Q8: How does the breakpoint window size affect detection sensitivity and specificity?**
- [ ] **Answer:** [Breakpoint window size changes the length of the breakpoint k-mer which is used for matching against reads. Larger windows are more specific but may miss reads that don't span the full window. Smaller windows are more sensitive but may have more false positives. Above a certain length, however, the differences are expected to be negligible.]

**Q9: What are the limitations of string matching compared to alignment-based methods?**
- [ ] **Answer:** [String matching as implemented requires a perfect match of the breakpoint k-mer to the read. This is a stringent requirement and may miss reads that are not a perfect match. Alignment or mismatch-tolerant matching will be slightly more sensitive but will be vastly slower. The accuracy of short-read sequencing is generally high enough that string matching is sufficient for the purpose of the tool.]

---

## Orientation Checklist

Before starting work on FUSILLI, an agent or editor should understand:

### Biological Understanding
- [ ] What fusion proteins are and how they're structured
- [ ] What DMS experiments are and why variant counting matters
- [ ] What breakpoints represent biologically
- [ ] How the fusion library is constructed (see Q1-Q2 above)
- [ ] What the research goals are (see Q5-Q6 above)

### Technical Understanding
- [ ] Snakemake workflow structure (see [ARCHITECTURE.md](ARCHITECTURE.md))
- [ ] Data flow through the pipeline (see [ARCHITECTURE.md](ARCHITECTURE.md))
- [ ] Breakpoint k-mer generation logic (see [README.md](README.md))
- [ ] String matching detection method (see [ARCHITECTURE.md](ARCHITECTURE.md))
- [ ] Output file formats (see [README.md](README.md))

### Code Understanding
- [ ] Project structure and file organization (see [README.md](README.md))
- [ ] Code style conventions (see [TECHNICAL.md](TECHNICAL.md))
- [ ] Testing patterns (see [TECHNICAL.md](TECHNICAL.md))
- [ ] Configuration system (see [README.md](README.md) and [ARCHITECTURE.md](ARCHITECTURE.md))

### Documentation Understanding
- [ ] Where to find different types of information:
  - [ ] Project overview: [README.md](README.md)
  - [ ] System design: [ARCHITECTURE.md](ARCHITECTURE.md)
  - [ ] Code patterns: [TECHNICAL.md](TECHNICAL.md)
  - [ ] AI collaboration: [AI_ASSISTED_CODING.md](AI_ASSISTED_CODING.md)
  - [ ] Project management: [PROJECT_MANAGEMENT.md](PROJECT_MANAGEMENT.md)
  - [ ] Current status: [status.md](status.md)
  - [ ] Tasks: [tasks/tasks.md](tasks/tasks.md)

---

## Quick Reference for Common Tasks

### Adding a New Feature
1. Review [tasks/tasks.md](tasks/tasks.md) for task requirements
2. Check [ARCHITECTURE.md](ARCHITECTURE.md) for extension points
3. Follow patterns in [TECHNICAL.md](TECHNICAL.md)
4. See [AI_ASSISTED_CODING.md](AI_ASSISTED_CODING.md) for collaboration patterns

### Fixing a Bug
1. Check error logs and messages
2. Review [TECHNICAL.md](TECHNICAL.md) error handling patterns
3. Check [status.md](status.md) for known issues
4. Review similar code sections for patterns

### Understanding the Pipeline
1. Start with [README.md](README.md) for overview
2. Review [ARCHITECTURE.md](ARCHITECTURE.md) for system design
3. Check [TECHNICAL.md](TECHNICAL.md) for implementation details
4. Review this document for biological context

### Working with Configuration
1. See [README.md](README.md) Configuration Reference section
2. Check [ARCHITECTURE.md](ARCHITECTURE.md) Configuration Layer section
3. Review config schema in `workflow/schemas/`

---

## Additional Resources

- **Main Entry Point:** [README.md](README.md)
- **System Architecture:** [ARCHITECTURE.md](ARCHITECTURE.md)
- **Code Patterns:** [TECHNICAL.md](TECHNICAL.md)
- **AI Collaboration:** [AI_ASSISTED_CODING.md](AI_ASSISTED_CODING.md)
- **Project Management:** [PROJECT_MANAGEMENT.md](PROJECT_MANAGEMENT.md)
- **Current Status:** [status.md](status.md)
- **Tasks:** [tasks/tasks.md](tasks/tasks.md)