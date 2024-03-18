# Config

The configs for all rules in this module follows the same format for settings: 

```yaml

rule_name:
   tag1: settings1
   tag2: settings2
```

# Example

all rule in Snakefile
```python

rule all:
   input:
       variants/sample1.include.BRCA1_EXONS.vcf.gz,
       variants/sample1.include.BRCA2_EXONS.vcf.gz,
       variants/sample1.exclude.INTRONS.vcf.gz,
       variants/sample1.bcftools_view.PASS.vcf.gz,
       variants/sample1.filter.GERMLINE.vcf.gz,
       
```

**config.yaml**
```yaml

# Triggered for the following files:
# - my_vcf.include.BRCA1_EXONS.vcf
# - my_vcf.include.BRCA2_EXONS.vcf
bcftools_filter_include_region:
   BRCA1_EXONS: "/data/design/brca1.bed"
   BRCA2_EXONS: "/data/design/brca2.bed"

# Triggered for my_vcf.exclude.INTRONS.vcf
bcftools_filter_exclude_region:
   INTRONS: "/data/design/introns.bed"


# Trigger for my_vcf.bcftools_view.PASS.vcf
bcftools_view:
   PASS: " --apply-filter PASS "

# Triggered for my_vcf.filter.GERMLINE.vcf
filter_vcf:
   GERMLINE: "config/filter.yaml"

```