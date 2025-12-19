# Coxsackievirus A10 Nextstrain Analysis

This repository provides a comprehensive Nextstrain analysis of Coxsackievirus A10. You can choose to perform either a VP1 run (>=600 base pairs) or a whole genome run (>=6400 base pairs).

For those unfamiliar with Nextstrain or needing installation guidance, please refer to the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/).

### Enhancing the Analysis

If there are some metadata that should be added, they can be included in the `data/meta_manually_added.tsv`, such as for example, published clades.

## Repository Organization

This repository includes the following directories and files:

  - `scripts`: Custom Python scripts called by the snakefile.
  - `snakefile`: The entire computational pipeline, managed using Snakemake. Snakemake documentation can be found here.
  - `ingest`: Contains Python scripts and the snakefile for automatic downloading of EV-A71 sequences and metadata.
  - `vp1`: Sequences and configuration files for the VP1 run.
  - `whole_genome`: Sequences and configuration files for the whole genome run.
  - `data`: Files needed for running the snakefile.

### Configuration Files

The `config`, `vp1/config`, and `whole_genome/config` directories contain necessary configuration files:

  - `colors.tsv`: Color scheme
  - `config.yaml`: Configuration file
  - `geo_regions.tsv`: Geographical locations
  - `lat_longs.tsv`: Latitude data
  - `dropped_strains.txt`: Dropped strains
  - `include.txt`: Samples to include
  - `reference_sequence.gb`: Reference sequence
  - `clades_genome.tsv`: Clade definitions
  - `auspice_config.json`: Auspice configuration file

> ⚠️ **Note:** The reference sequence used is [Kowalik, accession number AY421767](https://www.ncbi.nlm.nih.gov/nuccore/AY421767), from 1950. You need to add the year in the `data/metadata.tsv`, in the `date` column as 1950-XX-XX.

## TODOs

- [ ] Automate the inclusion of the reference collection date. If its inferred, is set to 2022 and affects all the tree structure.
- [X] Remove the `date` 'Color by'. Not found in the `auspice_config.json` files. → This field was included in the snakefile in the `rule export`, in the `--color-by-metadata {params.coloring_fields}`
