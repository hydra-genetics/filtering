__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2022, Jonas Almlöf, Patrik Smeds"
__email__ = "jonas.almlof@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule filter_vcf:
    input:
        vcf="{file}.vcf.gz",
        vcf_index="{file}.vcf.gz.tbi",
    output:
        vcf=temp("{file}.filter.{tag}.vcf"),
    params:
        filter_config=lambda wildcards: config["filter_vcf"][wildcards.tag],
    log:
        "{file}.filter.{tag}.log",
    benchmark:
        repeat("{file}.filter.{tag}.benchmark.tsv", config.get("filter_vcf", {}).get("benchmark_repeats", 1))
    threads: config.get("filter_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("filter_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("filter_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("filter_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("filter_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("filter_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("filter_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/filter_vcf.yaml"
    message:
        "{rule}: filter vcf {input.vcf} based on annotations"
    script:
        "../scripts/filter_vcf.py"
