# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2022, Jonas Almlöf, Patrik Smeds"
__email__ = "jonas.almlof@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule soft_filter_vcf:
    input:
        vcf="{file}.vcf.gz",
        vcf_index="{file}.vcf.gz.tbi",
    output:
        vcf=temp("{file}.soft_filter.vcf"),
    params:
        filter_config=config.get("soft_filter_vcf", {})["filter_config"],
    log:
        "{file}.soft_filter.log",
    benchmark:
        repeat("{file}.soft_filter.benchmark.tsv", config.get("soft_filter_vcf", {}).get("benchmark_repeats", 1))
    threads: config.get("soft_filter_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("soft_filter_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("soft_filter_vcf", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("soft_filter_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("soft_filter_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("soft_filter_vcf", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("soft_filter_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/filter_vcf.yaml"
    message:
        "{rule}: soft filter vcf {wildcards.file} based on annotations"
    script:
        "../scripts/filter_vcf.py"


rule hard_filter_vcf:
    input:
        vcf="{file}.vcf.gz",
        vcf_index="{file}.vcf.gz.tbi",
    output:
        vcf=temp("{file}.hard_filter.vcf"),
    params:
        filter_config=config.get("hard_filter_vcf", {})["filter_config"],
    log:
        "{file}.hard_filter.log",
    benchmark:
        repeat("{file}.hard_filter.benchmark.tsv", config.get("hard_filter_vcf", {}).get("benchmark_repeats", 1))
    threads: config.get("hard_filter_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("hard_filter_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("hard_filter_vcf", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("hard_filter_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("hard_filter_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("hard_filter_vcf", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("hard_filter_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/filter_vcf.yaml"
    message:
        "{rule}: hard filter vcf {wildcards.file} based on annotations"
    script:
        "../scripts/filter_vcf.py"
