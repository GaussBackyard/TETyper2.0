# 🧬💻 TETyper 2.0: Containerized Transposable Element Typing Pipeline

[![Version](https://img.shields.io/badge/version-2.0.0-blue.svg)]()
[![Status](https://img.shields.io/badge/status-Active%20Development-success.svg)]()
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)

---

## Overview

TETyper 2.0 is a containerized bioinformatics pipeline for characterizing transposable element (TE) from paired-end whole-genome sequencing (WGS) data. Building on the original TETyper 1.0 (Sheppard, 2018), version 2.0 delivers a fully containerized, multi-TE, multi-sample workflow with integrated variant calling.

Given paired FASTQ reads and reference sequences for your TE of interest, the pipeline:

- **Maps reads** to the TE's inverted repeat sequences (IRL/IRR) to identify insertion sites
- **Extracts** flanking sequences and **links** them to host contigs obtained from de novo assembly
- **Calls SNPs, indels, and structural variants** by mapping reads to the full transposon reference sequence (Case 2 only)
- **Scales to batch processing** of multiple samples and multiple TEs in a single run


---

## What's New in Version 2.0

- **Sensitive TE Detection**: IR-based read mapping and flanking sequence extraction at insertion sites
- **Flank-Contig Linkage**: [`BLAST+`](https://blast.ncbi.nlm.nih.gov/Blast.cgi)-based linking of extracted flanking sequences to [`SPAdes`](https://github.com/ablab/spades) assembled contigs
- **Updated Variant Calling**: SNP/indel detection and structural variant analysis via [`BCFtools`](https://samtools.github.io/bcftools/)
- **Improved Scalability**: Simultaneous analysis of multi-samples and multi-TEs in one run
- **Comprehensive Reporting**: JSON/TSV output formats with structured intermediate files for downstream analysis
- **Containerized Deployment**: All dependencies bundled in a single Apptainer/Singularity images


---

## Quick Start

All tools and dependencies (BWA, SAMtools, SPAdes, BLAST+, BCFtools, and supporting Python packages) are bundled into a single Apptainer image (`tetyper2.sif`). No separate installation of bioinformatics tools is required.

### Prerequisites

- **Linux** (local workstation or HPC cluster)
- **Apptainer ≥ 1.0** (formerly Singularity) — [installation guide](https://apptainer.org/docs/admin/main/installation.html)
- **`tetyper2.sif`** — the TETyper 2.0 container image (see [Packages](https://github.com/GaussBackyard?tab=packages&repo_name=TETyper2.0))
- Paired FASTQ files (plain or gzipped, `.fq`, `.fq.gz`, or `.fastq.gz`)
- Reference FASTA files for your TE (IRL, IRR, and optionally the full TE sequence)


### Smoke Test

Confirm the container is functional:

```bash
apptainer run --cleanenv tetyper2.sif --help
```

### Obtaining the Container Image

Download `tetyper2.sif` from the [Packages page](https://github.com/GaussBackyard?tab=packages&repo_name=TETyper2.0). Confirm the file is present before running:

```bash
ls -lh /path/to/tetyper2.sif
```

### Alternative: Source Installation (Under Construction)

For developers who want to run without the container


---

## Two Analysis Modes

TETyper 2.0 supports two analysis cases/modes, determined by the reference files you provide:

| Mode | Input | Pipeline Steps | Variant Output |
|------|-------|---------------|----------------|
| **Case 1** | IRL FASTA + IRR FASTA only | Map to IR boundaries, extract flanks and link to assembled contigs | None |
| **Case 2** | IRL FASTA + IRR FASTA + whole-TE FASTA | All Case 1 steps + map to full TE, assemble TE from BAM, call SNPs/indels/SVs | VCF files in `tn_variants/` |

> **IRL / IRR:** FASTA files containing the left and right inverted repeat sequences of your TE.

> **Whole-TE FASTA:** the complete reference sequence of the transposon.

> [!TIP]
> Use Case 1 if you don't have a good reference sequence and only have inverted repeat sequences. Use Case 2 if you have a complete transposon reference (e.g., from GenBank) and need full variant characterization.

---

## Basic Usage

### Container Command Pattern

```bash
apptainer run --cleanenv \
    -B /path/to/your/data:/path/to/your/data \
    tetyper2.sif \
    [TETyper arguments]
```

The `-B` flag **binds** (mounts) a host directory into the container. Use the same path on both sides for simplicity. The container must be able to see all input files and write to the output directory via this bind path.

---

## Example Commands

### Example 1 — Single Sample, Case 1 (IRL/IRR Only)

Identify TE insertion flanks for `SAMPLE01` using Tn2 inverted repeat references only.

```bash
apptainer run --cleanenv \
    -B /data/myproject:/data/myproject \
    tetyper2.sif \
    --reads1     /data/myproject/reads/SAMPLE01_R1.fq.gz \
    --reads2     /data/myproject/reads/SAMPLE01_R2.fq.gz \
    --sample-id  SAMPLE01 \
    --te-irl     /data/myproject/refs/Tn2_IRL.fasta \
    --te-irr     /data/myproject/refs/Tn2_IRR.fasta \
    --te-id      Tn2 \
    --flank-length 5 \
    --threads    8 \
    --memory-limit 16 \
    --output-dir /data/myproject/results/SAMPLE01
```

### Example 2 — Single Sample, Case 2 (IRL/IRR + Whole-TE FASTA)

Run the full pipeline including SNP/indel/SV calling against the complete TE reference. Add `--te-fasta` to trigger Case 2.

```bash
apptainer run --cleanenv \
    -B /data/myproject:/data/myproject \
    tetyper2.sif \
    --reads1     /data/myproject/reads/SAMPLE01_R1.fq.gz \
    --reads2     /data/myproject/reads/SAMPLE01_R2.fq.gz \
    --sample-id  SAMPLE01 \
    --te-irl     /data/myproject/refs/Tn2_IRL.fasta \
    --te-irr     /data/myproject/refs/Tn2_IRR.fasta \
    --te-id      Tn2 \
    --te-fasta   /data/myproject/refs/Tn2.fasta \
    --flank-length 5 \
    --threads    8 \
    --memory-limit 16 \
    --output-dir /data/myproject/results/SAMPLE01
```

### Example 3 — Single Sample, Multiple TEs (Case 2)

Type Tn2, Tn4401, and Tn2e simultaneously. Repeat `--te-irl`, `--te-irr`, `--te-id`, and `--te-fasta` for each TE in positional order.

> [!TIP]
> **Positional pairing rule:** The N-th `--te-irl`, `--te-irr`, `--te-id`, and `--te-fasta` all refer to the same TE. Tn2e shares its IRL/IRR sequences with Tn2, but requires its own whole-TE reference via `--te-fasta`.

```bash
apptainer run --cleanenv \
    -B /data/myproject:/data/myproject \
    tetyper2.sif \
    --reads1     /data/myproject/reads/SAMPLE01_R1.fq.gz \
    --reads2     /data/myproject/reads/SAMPLE01_R2.fq.gz \
    --sample-id  SAMPLE01 \
    --te-irl     /data/myproject/refs/Tn2_IRL.fasta \
    --te-irl     /data/myproject/refs/Tn4401_IRL.fasta \
    --te-irl     /data/myproject/refs/Tn2_IRL.fasta \
    --te-irr     /data/myproject/refs/Tn2_IRR.fasta \
    --te-irr     /data/myproject/refs/Tn4401_IRR.fasta \
    --te-irr     /data/myproject/refs/Tn2_IRR.fasta \
    --te-id      Tn2 \
    --te-id      Tn4401 \
    --te-id      Tn2e \
    --te-fasta   /data/myproject/refs/Tn2.fasta \
    --te-fasta   /data/myproject/refs/Tn4401.fasta \
    --te-fasta   /data/myproject/refs/Tn2e.fasta \
    --flank-length 5 \
    --threads    8 \
    --memory-limit 16 \
    --output-dir /data/myproject/results/SAMPLE01_multiTE
```

### Example 4 — JSON Config File (Reusable TE Definition)

Define a TE once in a JSON config and reuse it across runs. Including the `"tn_fasta"` key automatically triggers Case 2.

**`my_Tn2.json`:**

```json
{
    "name": "Tn2",
    "irl_fasta": "/data/myproject/refs/Tn2_IRL.fasta",
    "irr_fasta": "/data/myproject/refs/Tn2_IRR.fasta",
    "tn_fasta":  "/data/myproject/refs/Tn2.fasta"
}
```

```bash
apptainer run --cleanenv \
    -B /data/myproject:/data/myproject \
    tetyper2.sif \
    --reads1     /data/myproject/reads/SAMPLE01_R1.fq.gz \
    --reads2     /data/myproject/reads/SAMPLE01_R2.fq.gz \
    --sample-id  SAMPLE01 \
    --te-config  /data/myproject/my_Tn2.json \
    --flank-length 5 \
    --threads    8 \
    --memory-limit 16 \
    --output-dir /data/myproject/results/SAMPLE01
```

---

## Command-Line Reference

### Input Options

| Flag | Type | Description |
|------|------|-------------|
| `--reads1` / `-1` | str | Forward reads FASTQ (`.fq`, `.fq.gz`, `.fastq.gz`) |
| `--reads2` / `-2` | str | Reverse reads FASTQ |
| `--sample-id` | str | Sample identifier (default: stem of `--reads1`) |
| `--sample-sheet` | str | CSV batch file with columns `sample_id,reads1,reads2` |

### TE Configuration Options

Three mutually additive modes can be combined in a single run:

| Flag | Type | Description |
|------|------|-------------|
| `--te-name` | str, repeatable | Built-in TE name (`Tn2` or `Tn4401`); uses hardcoded IRL/IRR sequences; Case 1 only |
| `--te-config` | str, repeatable | Path to a JSON config file; supports `irl_fasta`/`irr_fasta` or `irl_sequence`/`irr_sequence`; include `"tn_fasta"` to trigger Case 2 |
| `--te-irl` | str, repeatable | IRL reference FASTA; must be paired positionally with `--te-irr` and `--te-id` |
| `--te-irr` | str, repeatable | IRR reference FASTA; positionally paired with `--te-irl` |
| `--te-id` | str, repeatable | Name for the TE defined by `--te-irl`/`--te-irr`; defaults to `custom_TE_1`, `custom_TE_2`, ... |
| `--te-fasta` | str, repeatable | **Whole-transposon reference FASTA**; positionally paired with `--te-irl`/`--te-irr`; triggers Case 2 |

### Processing Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--flank-length` | 10 | Base pairs of flanking sequence to extract at each IR boundary |
| `--min-mapped-length` | 40 | Minimum bp a read must align to be counted |
| `--min-reads` | 10 | Minimum read support required for a flank call |
| `--min-each-strand` | 3 | Minimum supporting reads per strand |
| `--min-quality` | 20 | Minimum Phred base quality score |
| `--threads` | 4 | CPU threads for BWA, SAMtools, SPAdes, BLAST, BCFtools |
| `--memory-limit` | 8 | SPAdes memory cap in GB |

### Output Options

| Flag | Default | Description |
|------|---------|-------------|
| `--output-dir` / `-o` | `.` | Root directory for all output files |
| `--output-format` | `both` | Output format: `json`, `tsv`, or `both` |

### Pipeline Options

| Flag | Description |
|------|-------------|
| `--assembly` | Path to pre-built contigs FASTA; skips isolate SPAdes assembly |

### Recommended Settings (Illumina 150 bp paired-end, for detecting Tn4401/Tn2)

```
--flank-length     5
--min-mapped-length 40
--min-reads        10
--min-each-strand   3
--min-quality      20
--threads           8
--memory-limit     32
```

---

## Built-in TE Reference Sequences

| TE | Family | IRL | IRR |
|----|--------|-----|-----|
| **Tn2** | Tn3 | `GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGG` | `GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGCA` |
| **Tn4401** | Tn3 | `GGGGTTCTAAGCGGGAATCCCAGAAAATTCCGTCATTCCG` | `GGGGGGGTAAGCGGGAECCCCAGAAAATTCCGCCATTCCG` |

Custom TEs can be provided via `--te-irl`/`--te-irr` FASTA files or via a JSON config with `"irl_sequence"`/`"irr_sequence"` keys.

---

## Output Structure

After a run completes, the output directory contains:

```
results/SAMPLE01/
├── SAMPLE01_<timestamp>.log          # Full run log (all pipeline steps)
├── SAMPLE01_Tn2.json                 # Results in JSON format
├── SAMPLE01_Tn2.tsv                  # Results in TSV format
├── references/                       # IRL/IRR FASTA copies (sequence-based input only)
└── SAMPLE01/                         # Working directory
    ├── assembly/                     # Whole-genome SPAdes assembly
    └── Tn2/
        ├── IR_mapping/               # BWA alignments to IRL/IRR 
        ├── flank_extract/            # Extracted flank sequences 
        ├── flank_contig_link/        # BLAST-based flank-to-contig links 
        ├── tn_mapping/               # Reads mapped to whole-TE (Case 2)
        ├── tn_assembly/              # TE-specific SPAdes assembly (Case 2)
        └── tn_variants/              # VCF files, SNP/SV calls (Case 2)
```

When running multiple TEs, each TE receives its own subdirectory (`Tn2/`, `Tn4401/`, `Tn2e/`) with fully isolated intermediate files.

### Result File Contents

**JSON output** (`SAMPLE01_Tn2.json`) contains:

- `sample_id`, `te_name`
- `flank_results` — left and right flank sequences with supporting read counts
- `contig_links` — linkages between flank sequences and assembled contigs
- `tn_assembly_stats` — SPAdes assembly statistics *(Case 2 only)*
- `variant_summary` — total SNPs, indels, and structural variants *(Case 2 only)*
- `summary` — top-level pass/fail status and aggregate counts

**TSV output** (`SAMPLE01_Tn2.tsv`) — tab-delimited flat representation of flank results, suitable for spreadsheet import or downstream scripting.

> [!NOTE]
> **Skipped samples:** If a sample has fewer than 50 properly paired reads mapping to the TE reference, the `tn_assembly_stats` field will contain `"skipped_reason": "too_few_paired_pairs"`. This indicates the TE is absent or at very low coverage — it is expected behavior, not a pipeline error. Check the `.log` file for the exact paired read count.

---

## Running on an HPC Cluster (SLURM)

The container runs without modification on any HPC cluster with Apptainer or Singularity available. Below is a minimal SLURM batch script template:

```bash
#!/bin/bash
#SBATCH --job-name=tetyper2
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=tetyper2_%j.log

module load apptainer   # or: module load singularity

CONTAINER=/path/to/tetyper2.sif
DATA=/path/to/your/data

apptainer run --cleanenv \
    -B ${DATA}:${DATA} \
    ${CONTAINER} \
    --reads1       ${DATA}/reads/SAMPLE01_R1.fq.gz \
    --reads2       ${DATA}/reads/SAMPLE01_R2.fq.gz \
    --sample-id    SAMPLE01 \
    --te-irl       ${DATA}/refs/Tn2_IRL.fasta \
    --te-irr       ${DATA}/refs/Tn2_IRR.fasta \
    --te-id        Tn2 \
    --te-fasta     ${DATA}/refs/Tn2.fasta \
    --threads      8 \
    --memory-limit 32 \
    --output-dir   ${DATA}/results/SAMPLE01
```

Submit with:

```bash
sbatch run_tetyper2.sh
```

---



## Troubleshooting

### Container File Not Found

```
FATAL: ... lstat /resolved/path/tetyper2.sif: no such file or directory
```

Confirm the `.sif` file exists and always use its absolute path:

```bash
ls -lh /path/to/tetyper2.sif
apptainer run --cleanenv /path/to/tetyper2.sif --help
```

### Input Files Not Found Inside the Container

If the pipeline reports missing reads or reference files, the bind path is likely the issue. Always bind the common parent directory containing all input files:

```bash
# If your data is at /data/myproject/reads/ and /data/myproject/refs/
-B /data/myproject:/data/myproject
```

### Python Import Errors at Runtime

Always use `--cleanenv` to prevent the host's Python environment from interfering with the container's bundled conda environment:

```bash
apptainer run --cleanenv ...
```

### Too Few Reads / Skipped Assembly

If results show `"skipped_reason": "too_few_paired_pairs"`, the TE is likely absent or present at very low coverage in that sample. This is a valid result, not a crash. Check the run log (`.log` file) for the exact paired read count.


---

## FAQ

**Q: Can I use single-end reads?**  
A: Not currently. The pipeline is optimized for paired-end Illumina data.

**Q: What TE references are supported out of the box?**  
A: Tn2 and Tn4401 (both Tn3-family) are bundled as built-in sequences via `--te-name`. Any other TE can be supplied via `--te-irl`/`--te-irr` FASTA files or a JSON config.

**Q: Can I type multiple TEs from the same sample simultaneously?**  
A: Yes. Repeat `--te-irl`, `--te-irr`, `--te-id`, and (for Case 2) `--te-fasta` for each TE, following the positional pairing rule.

**Q: What if Tn2e shares IRL/IRR with Tn2?**  
A: Reuse `Tn2_IRL.fasta` and `Tn2_IRR.fasta` for the Tn2e entry, but supply the distinct `Tn2e.fasta` as `--te-fasta` for that positional slot.

**Q: How do I speed up runs?**  
A: Pass a pre-built contigs FASTA via `--assembly <path/to/contigs.fasta>` to skip the SPAdes assembly step.

---

## Issues / Feedback

Feedback is essential for improving TETyper 2.0. 

Please submit feedback or suggestions via [GitHub Issues](https://github.com/GaussBackyard/TETyper2.0/issues) or email [Thomas Li](mailto:gauss.backyard@gmail.com).

---


## Citation

If you use TETyper 2.0 in your research, please cite:

> Li, T.J. *et al.* TETyper 2.0: a containerized pipeline for transposable element typing from whole-genome sequencing data. *(in preparation, 2026)*

The original TETyper 1.0 was developed by Anna Sheppard (2018).

---



## License

This repository is licensed under the GPL-3.0 License. 

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](LICENSE) file for details.



---

**Last Updated:** 2026-05-26  
**Maintainer:** Thomas Jiaxian Li  
**Repository:** [https://github.com/GaussBackyard/TETyper2.0](https://github.com/GaussBackyard/TETyper2.0) 
