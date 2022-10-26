# <img src="https://github.com/hydra-genetics/biomarker/blob/develop/images/hydragenetics.png" width=40 /> hydra-genetics/filtering

Collection of variant filters

![Lint](https://github.com/hydra-genetics/filtering/actions/workflows/lint.yaml/badge.svg?branch=develop)
![Snakefmt](https://github.com/hydra-genetics/filtering/actions/workflows/snakefmt.yaml/badge.svg?branch=develop)
![snakemake dry run](https://github.com/hydra-genetics/filtering/actions/workflows/snakemake-dry-run.yaml/badge.svg?branch=develop)
![integration test](https://github.com/hydra-genetics/filtering/actions/workflows/integration.yaml/badge.svg?branch=develop)
![pycodestyle](https://github.com/hydra-genetics/filtering/actions/workflows/pycodestyl.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

The module consists of rules used to filter VCF files. Currently the VCFs can be filtered based on regions and annotation.

### Region based filtering

Regions filtering is based on bcftools and inculdes or excludes variants based on a bed file.

### Annotation based filtering

Filter criteria are defined in a .yaml file and can be [hard filters](https://github.com/hydra-genetics/filtering/blob/develop/.tests/integration/config_hard_filter.yaml) or [soft filters](https://github.com/hydra-genetics/filtering/blob/develop/.tests/integration/config_soft_filter.yaml) where a filter flag is added to the filter column in the VCF file. Example of annotations that can be used are those found in the format column or added by VEP. Filter criteria can be combined by logical operators. For matching annotations there is the possibility to use regular expressions in combination with for example exist. Also, NA-behavior can be specified.

## :heavy_exclamation_mark: Dependencies

In order to use this module, the following dependencies are required:

[![hydra-genetics](https://img.shields.io/badge/hydragenetics-0.15.0-blue)](https://github.com/hydra-genetics/)
[![pandas](https://img.shields.io/badge/pandas-1.3.1-blue)](https://pandas.pydata.org/)
[![python](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/)
[![snakemake](https://img.shields.io/badge/snakemake-7.13.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.0.0-blue)](https://sylabs.io/docs/)
[![drmaa](https://img.shields.io/badge/drmaa-0.7.9-blue)](https://pypi.org/project/drmaa/)
[![tabulate](https://img.shields.io/badge/tabulate-0.8.10-blue)](https://pypi.org/project/tabulate/)

## :school_satchel: Preparations

### Sample and unit data

Input data should be added to [`samples.tsv`](https://github.com/hydra-genetics/prealignment/blob/develop/config/samples.tsv)
and [`units.tsv`](https://github.com/hydra-genetics/prealignment/blob/develop/config/units.tsv).
The following information need to be added to these files:

| Column Id | Description |
| --- | --- |
| **`samples.tsv`** |
| sample | unique sample/patient id, one per row |
| tumor_content | ratio of tumor cells to total cells |
| **`units.tsv`** |
| sample | same sample/patient id as in `samples.tsv` |
| type | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA |
| platform | type of sequencing platform, e.g. `NovaSeq` |
| machine | specific machine id, e.g. NovaSeq instruments have `@Axxxxx` |
| flowcell | identifer of flowcell used |
| lane | flowcell lane number |
| barcode | sequence library barcode/index, connect forward and reverse indices by `+`, e.g. `ATGC+ATGC` |
| fastq1/2 | absolute path to forward and reverse reads |
| adapter | adapter sequences to be trimmed, separated by comma |

## :white_check_mark: Testing

The workflow repository contains a small test dataset `.tests/integration` which can be run like so:

```bash
cd .tests/integration
$ snakemake -s ../../Snakefile -j1 --configfile config.yaml --use-singularity
```

## :rocket: Usage

To use this module in your workflow, follow the description in the
[snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
Add the module to your `Snakefile` like so:

```bash
module biomarker:
    snakefile:
        github(
            "hydra-genetics/biomarker",
            path="workflow/Snakefile",
            tag="v0.1.0",
        )
    config:
        config


use rule * from biomarker as biomarker_*
```

### Compatibility

Latest:
 - annotation:v0.1.0
 - snv_indel:v0.2.0

See [COMPATIBLITY.md](../master/COMPATIBLITY.md) file for a complete list of module compatibility.

### Input files

| File | Description |
|---|---|
| ***`hydra-genetics/annotation`*** |
| `annotation/{file}.vcf.gz` | annotated vcf |
| ***`hydra-genetics/snv_indel data`*** |
| `cnv_sv/{file}.vcf.gz` | non-annotated vcf |

### Output files

The following output files should be targeted via another rule:

| File | Description |
|---|---|
| `{file}.include.{tag}.vcf.gz` | vcf filtered by bcftools include |
| `{file}.exclude.{tag}.vcf.gz` | vcf filtered by bcftools exclude |
| `{file}.filter.{tag}.vcf.gz` | vcf filtered based on annotations |

## :judge: Rule Graph

### Biomarker

![rule_graph](images/biomarker.svg)
