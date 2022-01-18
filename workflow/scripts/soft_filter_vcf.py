import logging
import yaml

from collections import OrderedDict
from hydra_genetics.utils.io import utils
from pysam import VariantFile

in_vcf = snakemake.input.vcf
out_vcf = snakemake.output.vcf
filter_yaml_file = snakemake.params.filter_config

variants = VariantFile(in_vcf)
log = logging.getLogger()


def is_float(element) -> bool:
    try:
        float(element)
        return True
    except ValueError:
        return False


filters = {"filters": []}
if filter_yaml_file is not None:
    log.info("Process yaml for: {}".format(filter_yaml_file))
    with open(filter_yaml_file) as file:
        filters = yaml.load(file, Loader=yaml.FullLoader)

log.info("Checking yaml file parameters")
for filter in filters["filters"]:
    fields_len = 0
    tests_len = 0
    test_operators_len = 0
    if "annotation" in filters["filters"][filter]:
        if filters["filters"][filter]["annotation"] not in ["VEP", "FORMAT"]:
            raise Exception("Unknown annotation: %s" % filters["filters"][filter]["annotation"])
    else:
        raise Exception("No annotation entry for %s" % filter)
    if "fields" in filters["filters"][filter]:
        fields_len = len(filters["filters"][filter]["fields"])
    else:
        raise Exception("No fields entry for %s" % filter)
    if "tests" in filters["filters"][filter]:
        tests_len = len(filters["filters"][filter]["tests"])
    else:
        raise Exception("No tests entry for %s" % filter)
    if "test_operators" in filters["filters"][filter]:
        test_operators_len = len(filters["filters"][filter]["test_operators"])
        for operator in filters["filters"][filter]["test_operators"]:
            if operator not in ["<", ">", "=", "!="]:
                raise Exception("Unknown test operator: %s" % operator)
    else:
        raise Exception("No test_operators entry for %s" % filter)
    if not(fields_len == tests_len and tests_len == test_operators_len):
        raise Exception("Different length for fields, tests, and test_operators %s" % filter)
    if "logical_operator" in filters["filters"][filter]:
        if filters["filters"][filter]["logical_operator"] not in ["or", "and"]:
            raise Exception("Unknown logical operator: %s" % logical_operator)
    if "filter_string" in filters["filters"][filter]:
        filter_text = "Failed %s filter" % filter
        if "description" in filters["filters"][filter]:
            filter_text = filters["filters"][filter]["description"]
        variants.header.filters.add(filters["filters"][filter]["filter_string"], None, None, filter_text)
    else:
        raise Exception("No filter_string entry for %s" % filter)

log.info("Process vcf header: {}".format(in_vcf))
for record in variants.header.records:
    if record.type == "INFO":
        if record['ID'] == "CSQ":
            log.info(" -- found vep information: {}".format(in_vcf))
            log.debug(" -- -- {}".format(record['Description'].split("Format: ")[1].split("|")))
            vep_fields = {v: c for c, v in enumerate(record['Description'].split("Format: ")[1].split("|"))}
            annotation_extractor = utils.get_annoation_data_vep(vep_fields)

vcf_out = VariantFile(out_vcf, 'w', header=variants.header)

log.info("Processing variants")
for variant in variants:
    for filter in filters["filters"]:
        i = 0
        fail_filter_list = []
        for field in filters["filters"][filter]["fields"]:
            data = ""
            if filters["filters"][filter]["annotation"] == "VEP":
                data = annotation_extractor(variant, field)
            elif filters["filters"][filter]["annotation"] == "FORMAT":
                if (
                    type(variant.samples[0][field]) is tuple and
                    "indexes" in filters["filters"][filter] and
                    filters["filters"][filter]["indexes"][i] != ""
                ):
                    data = variant.samples[0][field][filters["filters"][filter]["indexes"][i]]
                else:
                    data = variant.samples[0][field]
            if filters["filters"][filter]["test_operators"][i] == "<":
                if data == "":
                    fail_filter_list.append(True)
                else:
                    fail_filter_list.append(float(data) < float(filters["filters"][filter]["tests"][i]))
            elif filters["filters"][filter]["test_operators"][i] == ">":
                if data == "":
                    fail_filter_list.append(True)
                else:
                    fail_filter_list.append(float(data) > float(filters["filters"][filter]["tests"][i]))
            elif filters["filters"][filter]["test_operators"][i] == "=":
                if is_float(data):
                    fail_filter_list.append(float(data) == float(filters["filters"][filter]["tests"][i]))
                else:
                    fail_filter_list.append(data == filters["filters"][filter]["tests"][i])
            elif filters["filters"][filter]["test_operators"][i] == "!=":
                if is_float(data):
                    fail_filter_list.append(float(data) != float(filters["filters"][filter]["tests"][i]))
                else:
                    fail_filter_list.append(data != filters["filters"][filter]["tests"][i])
            i += 1
        fail_filter = fail_filter_list[0]
        if len(fail_filter_list) > 1:
            if filters["filters"][filter]["logical_operator"] == "or":
                fail_filter = any(fail_filter_list)
            elif filters["filters"][filter]["logical_operator"] == "and":
                fail_filter = all(fail_filter_list)
        if fail_filter:
            variant.filter.add(filters["filters"][filter]["filter_string"])
    vcf_out.write(variant)
vcf_out.close()
