# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule background_filter:
    input:
        vcf="filtering/artifact_filter/{sample}_{type}.artifact_filter.vcf",
        background=config["reference"]["background"],
    output:
        vcf=temp("filtering/background_filter/{sample}_{type}.background_filter.vcf"),
    params:
        nr_min_sd=config.get("background_filter", {}).get("nr_min_sd", 5),
    log:
        "filtering/background_filter/{sample}_{type}.background_filter.vcf.log",
    benchmark:
        repeat(
            "filtering/background_filter/{sample}_{type}.background_filter.vcf.benchmark.tsv",
            config.get("background_filter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("background_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("background_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("background_filter", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("background_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("background_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("background_filter", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("background_filter", {}).get("container", config["default_container"])
    conda:
        "../envs/background_filter.yaml"
    message:
        "{rule}: Background filter vcf in filtering/background_filter/{wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/background_filter.py"
