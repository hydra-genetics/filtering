# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule filter_vcf_on_format:
    input:
        vcf="snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vep_annotated.vcf.gz",
    output:
        vcf=temp("filter/filter_vcf_on_format/{sample}_{type}.format_filt.vcf"),
    params:
        filter=config.get("filter_vcf_on_format", {}).get("filter", "-i 'FORMAT/DP > 100 & FORMAT/AD > 20 & FORMAT/AF > 0.05'"),
        extra=config.get("filter_vcf_on_format", {}).get("extra", "--mode '+' --soft-filter 'DP_AD_AF'"),
    log:
        "filter/filter_vcf_on_format/{sample}_{type}.log",
    benchmark:
        repeat(
            "filter/filter_vcf_on_format/{sample}_{type}.benchmark.tsv",
            config.get("filter_vcf_on_format", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("filter_vcf_on_format", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("filter_vcf_on_format", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("filter_vcf_on_format", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("filter_vcf_on_format", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("filter_vcf_on_format", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("filter_vcf_on_format", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("filter_vcf_on_format", {}).get("container", config["default_container"])
    conda:
        "../envs/filter_vcf_on_format.yaml"
    message:
        "{rule}: Filter vcf filter/filter_vcf_on_format/{wildcards.sample}_{wildcards.type} based on format"
    wrapper:
        "0.72.0/bio/bcftools/filter"
