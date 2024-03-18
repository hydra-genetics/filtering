# Softwares used in the filtering module

## [bcftools filter (include region)](https://samtools.github.io/bcftools/bcftools.html#filter)
Filter a VCF file by using a BED file to define the regions where variants should be kept.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bcftools__bcftools_filter_include_region#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bcftools__bcftools_filter_include_region#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bcftools_filter_include_region#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bcftools_filter_include_region#

## [bcftools filter (exclude regions)](https://samtools.github.io/bcftools/bcftools.html£filter)
Filter a VCF file based on a provided BED file that specifies regions where variants will be removed.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bcftools__bcftools_filter_exclude_region#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bcftools__bcftools_filter_exclude_region#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bcftools_filter_exclude_region#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bcftools_filter_exclude_region#


## [bcftools view](https://samtools.github.io/bcftools/bcftools.html#view)
Filter VCF files based on their content, including criteria such as flags (e.g. PASS) and variant types.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bcftools__bcftools_view#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bcftools__bcftools_view#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bcftools_view#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bcftools_view#

## [filter_vcf](https://github.com/hydra-genetics/filtering/blob/develop/workflow/scripts/filter_vcf.py)
In-house developed python script used to do more complex filtering on the fields of a variant.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__filter_vcf__filter_vcf#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__filter_vcf__filter_vcf#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__filter_vcf#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__filter_vcf#