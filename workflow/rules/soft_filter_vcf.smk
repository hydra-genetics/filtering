# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule soft_filter_vcf:
    input:
        vcf="snv_indels/ensemble_vcf/HD832.HES45_T.ensembled.vep_annotated.vcf.gz",
        vcf_index="snv_indels/ensemble_vcf/HD832.HES45_T.ensembled.vep_annotated.vcf.gz.tbi",
    output:
        vcf=temp("filtering/soft_filter_vcf/{sample}_{type}.soft_filt.vcf"),
    params:
        filter_config=config.get("soft_filter_vcf", {})["filter_config"],
    log:
        "filtering/soft_filter_vcf/{sample}_{type}.log",
    benchmark:
        repeat(
            "filtering/soft_filter_vcf/{sample}_{type}.benchmark.tsv",
            config.get("soft_filter_vcf", {}).get("benchmark_repeats", 1),
        )
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
        "../envs/soft_filter_vcf.yaml"
    message:
        "{rule}: Filter vcf filtering/soft_filter_vcf/{wildcards.sample}_{wildcards.type} based on annotations"
    script:
        "../scripts/soft_filter_vcf.py"
