# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


def get_include_bed_file(wildcards):
    return "-R %s" % config["bcftools_region_include_filter"][wildcards.tag]

def get_exclude_bed_file(wildcards):
    return "-T ^%s" % config["bcftools_region_exclude_filter"][wildcards.tag]

rule bcftools_region_include_filter:
    input:
        vcf="{file}.vcf.gz",
        tabix="{file}.vcf.gz.tbi",
    output:
        vcf=temp("{file}.include.{tag}.vcf.gz"),
    params:
        filter=get_include_bed_file,
        extra=config.get("bcftools_region_include_filter", {}).get("extra", ""),
    log:
        "{file}.include.{tag}.log",
    benchmark:
        repeat(
            "{file}.include.{tag}.benchmark.tsv",
            config.get("bcftools_region_include_filter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_region_include_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("bcftools_region_include_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_region_include_filter", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("bcftools_region_include_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_region_include_filter", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("bcftools_region_include_filter", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("bcftools_region_include_filter", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools_region_filter.yaml"
    message:
        "{rule}: Use bedtools to include variants in vcf overlapping bed: {wildcards.file}.include.{wildcards.tag}.vcf"
    wrapper:
        "v1.0.0/bio/bcftools/filter"


rule bcftools_region_exclude_filter:
    input:
        vcf="{file}.vcf.gz",
        tabix="{file}.vcf.gz.tbi",
    output:
        vcf=temp("{file}.exclude.{tag}.vcf.gz"),
    params:
        filter=get_exclude_bed_file,
        extra=config.get("bcftools_region_exclude_filter", {}).get("extra", ""),
    log:
        "{file}.exclude.{tag}.log",
    benchmark:
        repeat(
            "{file}.exclude.{tag}.benchmark.tsv",
            config.get("bcftools_region_exclude_filter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_region_exclude_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("bcftools_region_exclude_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_region_exclude_filter", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("bcftools_region_exclude_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_region_exclude_filter", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("bcftools_region_exclude_filter", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("bcftools_region_exclude_filter", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools_region_filter.yaml"
    message:
        "{rule}: Use bedtools to exclude variants in vcf overlapping bed: {wildcards.file}.exclude.{wildcards.tag}.vcf"
    wrapper:
        "v1.0.0/bio/bcftools/filter"
