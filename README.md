# Coxsackievirus A10 Nextstrain Analysis

This repository provides a comprehensive Nextstrain analysis of Coxsackievirus A10 (CVA10), suitable for both **VP1 run (≥600 bp)** and **whole genome run (≥6400 bp)** analyses.

For those unfamiliar with Nextstrain or needing installation guidance, see the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/).

---

## Repository Organization

- **ingest/**: Python scripts and Snakefile for automatic downloading of CVA10 sequences and metadata (GenBank), formatting, and preparation.
- **scripts/**: Custom Python scripts called by the main Snakefile.
- **snakefile**: The main computational pipeline, managed using Snakemake. [Snakemake docs](https://snakemake.readthedocs.io/en/stable/)
- **vp1/**: Sequences and configuration files for the **VP1 run**.
- **whole_genome/**: Sequences and configuration files for the **whole genome run**.
- **data/**: Static data files used by the workflow.

### Configuration Files

The following files may be found under `config/`, `vp1/config/`, or `whole_genome/config/`:
- `colors.tsv`: Color scheme
- `geo_regions.tsv`: Geographical regions
- `lat_longs.tsv`: Latitude/longitude locations
- `dropped_strains.txt`: List of excluded strains
- `clades_genome.tsv`: Virus clade definitions/assignments
- `reference_sequence.gb`: Reference sequence (GenBank format)
- `auspice_config.json`: Auspice display config

> ⚠️ **Note:** The reference sequence used is [Kowalik, accession number AY421767](https://www.ncbi.nlm.nih.gov/nuccore/AY421767), sampled in 1950.

---

## Quickstart

### 1. Install Nextstrain

Follow instructions in the [Nextstrain installation guide](https://docs.nextstrain.org/en/latest/guides/install/local-installation.html).

### 2. Data Preparation & Ingest Pipeline

CVA10 sequence and metadata ingestion is automated via the **`ingest/`** workflow.

**Prepare reference files:**

- Check `config/config.yaml` and confirm taxid (NCBI) is correct.
- Run the reference extraction script:
  ```bash
  python3 ingest/bin/generate_from_genbank.py --reference "AY421767.1" --output-dir whole_genome/config/
  ```
  - The script may prompt for reference protein/CDS selection. Typical codes: `[0]` `[product]` `[2]`.

- Ensure attributes in `data/references/pathogen.json` are up-to-date.
  - Reference: [Nextclade pathogen config docs](https://docs.nextstrain.org/projects/nextclade/en/stable/how-to/custom-reference.html)

**Run the ingest workflow:**
```bash
cd ingest
chmod +x ./vendored/*; chmod +x ./bin/*
snakemake --cores 4 all
```
This gathers/fetches and processes the latest public (and optionally, private) CVA10 sequence and metadata for downstream analysis.

### 3. Running a Build

In the project root:
```bash
snakemake --cores 9 all
```

**Or, for a specific result:**
- VP1:
  ```bash
  snakemake  auspice/coxsackievirus_A10_vp1.json --cores 9
  ```
- Whole genome:
  ```bash
  snakemake auspice/coxsackievirus_A10_genome.json --cores 9
  ```

### 4. Visualizing

To view build results in Auspice:
```bash
auspice view --datasetDir auspice
```

---

## Data Input

- **Automatic:** The `ingest/` pipeline fetches up-to-date public data. Private sequences or metadata (e.g. new submissions) can be incorporated by adding to `data/meta_manually_added.tsv`.
- **Manual:** Download from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), searching for “CVA10” or Taxid.

---

## Updating Vendored Scripts

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage ingest scripts in `ingest/vendored`. To update:
```bash
git subrepo pull ingest/vendored
```
See details in `ingest/vendored/README.md`.

---

## Troubleshooting & Notes

- *Reference collection date/inclusion*: Make sure the true collection date for the reference sequence is in your metadata to avoid downstream tree artifacts.
- *Colors*: If not found in Auspice, make sure your `auspice_config.json` does include them.

---

## Feedback

For questions or comments, please contact us per e-mail (eve-group[at]swisstph.ch) or open an [issue](https://github.com/hodcroftlab/coxsackievirus_a10/issues).

---

*For advanced users and further workflow options, see the documentation in each subdirectory (especially `ingest/README.md`).*
