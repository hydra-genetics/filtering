__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcftools_filter_include_region:
    input:
        vcf="{file}.vcf.gz",
        tabix="{file}.vcf.gz.tbi",
    output:
        vcf=temp("{file}.include.{tag}.vcf.gz"),
    params:
        filter=lambda wildcards: "-R %s" % config["bcftools_filter_include_region"][wildcards.tag],
        extra=config.get("bcftools_filter_include_region", {}).get("extra", ""),
    log:
        "{file}.include.{tag}.log",
    benchmark:
        repeat(
            "{file}.include.{tag}.benchmark.tsv",
            config.get("bcftools_filter_include_region", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_filter_include_region", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_filter_include_region", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_filter_include_region", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("bcftools_filter_include_region", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_filter_include_region", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_filter_include_region", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_filter_include_region", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools.yaml"
    message:
        "{rule}: Use bedtools to include variants in vcf overlapping bed: {output.vcf}"
    shell:
        "(bcftools filter "
        "{params.filter} "
        "{params.extra} "
        "{input.vcf} "
        "-o {output.vcf}) &> {log}"


rule bcftools_filter_exclude_region:
    input:
        vcf="{file}.vcf.gz",
        tabix="{file}.vcf.gz.tbi",
    output:
        vcf=temp("{file}.exclude.{tag}.vcf.gz"),
    params:
        filter=lambda wildcards: "-T ^%s" % config["bcftools_filter_exclude_region"][wildcards.tag],
        extra=config.get("bcftools_filter_exclude_region", {}).get("extra", ""),
    log:
        "{file}.exclude.{tag}.log",
    benchmark:
        repeat(
            "{file}.exclude.{tag}.benchmark.tsv",
            config.get("bcftools_filter_exclude_region", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_filter_exclude_region", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_filter_exclude_region", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_filter_exclude_region", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("bcftools_filter_exclude_region", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_filter_exclude_region", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_filter_exclude_region", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_filter_exclude_region", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools.yaml"
    message:
        "{rule}: use bedtools to exclude variants in vcf overlapping bed: {output.vcf}"
    shell:
        "(bcftools filter "
        "{params.filter} "
        "{params.extra} "
        "{input.vcf} "
        "-o {output.vcf}) &> {log}"


rule bcftools_view:
    input:
        vcf="{file}.vcf.gz",
    output:
        vcf=temp("{file}.bcftools_view.vcf.gz"),
    params:
        extra=config.get("bcftools_view", {}).get("extra", ""),
    log:
        "{file}.bcftools_view.vcf.log",
    benchmark:
        repeat(
            "{file}.bcftools_view.vcf.benchmark.tsv",
            config.get("bcftools_view", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_view", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_view", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_view", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_view", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_view", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_view", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_view", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools.yaml"
    message:
        "{rule}: Use bcftools view to get subset or filter {input.vcf}"
    shell:
        "v1.25.0/bio/bcftools/view"
