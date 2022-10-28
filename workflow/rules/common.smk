__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("7.8.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def compile_output_list(wildcards):
    files = {
        "annotation/background_annotation": [
            ".background_annotation.include.noexon1.vcf.gz",
            ".background_annotation.exclude.noexon1.vcf.gz",
            ".background_annotation.filter.snv_soft_filter.vcf",
            ".background_annotation.filter.snv_hard_filter.vcf",
        ]
    }
    output_files = [
        "%s/%s_%s%s" % (prefix, sample, unit_type, suffix)
        for prefix in files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix]
    ]
    return output_files
