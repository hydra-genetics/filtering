# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule artifact_filter:
    input:
        vcf="filtering/soft_filter_vcf/{sample}_{type}.soft_filt.vcf",
        artifacts=config["reference"]["artifacts"],
    output:
        vcf=temp("filtering/artifact_filter/{sample}_{type}.artifact_filter.vcf"),
    log:
        "filtering/artifact_filter/{sample}_{type}.artifact_filter.vcf.log",
    benchmark:
        repeat(
            "filtering/artifact_filter/{sample}_{type}.artifact_filter.vcf.benchmark.tsv",
            config.get("artifact_filter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("artifact_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("artifact_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("artifact_filter", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("artifact_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("artifact_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("artifact_filter", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("artifact_filter", {}).get("container", config["default_container"])
    conda:
        "../envs/artifact_filter.yaml"
    message:
        "{rule}: Artifact filter vcf in filtering/artifact_filter/{wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/artifact_filter.py"
