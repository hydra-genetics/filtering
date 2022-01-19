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

def parse_filter_expression(filter):
    def parse_helper(f_iterator):
        """
        function used to parse a string with parantes and create a nested fail_filter_list

        :param f_iterator: filter string that will be parsed
        :type f_iterator: string iterator

        :return: return a nested list
        """
        data = []
        filter_string = ''
        for item in f_iterator:
            if item == '(':
                result, closing_parantes = parse_helper(f_iterator)
                if len(result) == 1:
                    result = result[0]
                data.append(result)
                if not closing_parantes:
                    raise SyntaxError("Missing closing parantes!")
            elif item == ')':
                if len(filter_string) != 0:
                    data.append(filter_string)
                return data, True
            else:
                filter_string += item
                if filter_string.startswith(' or '):
                    data.append(lambda v1, v2: v1 or v2)
                    result, p_found = parse_helper(f_iterator)
                    if len(result) == 1:
                        result = result[0]
                    data.append(result)
                    return data, p_found
                elif filter_string.startswith(' and '):
                    data.append(lambda v1, v2: v1 and v2)
                    result, p_found = parse_helper(f_iterator)
                    if len(result) == 1:
                        result = result[0]
                    data.append(result)
                    return data, p_found
                elif filter_string.endswith(' or '):
                    data.append(filter_string[:-4])
                    data.append(lambda v1, v2: v1 or v2 )
                    result, p_found = parse_helper(f_iterator)
                    if len(result) == 1:
                        result = result[0]
                    data.append(result)
                    return data, p_found
                elif filter_string.endswith(' and '):
                    data.append(filter_string[:-5])
                    data.append(lambda v1, v2: v1 and v2)
                    result, p_found = parse_helper(f_iterator)
                    if len(result) == 1:
                        result = result[0]
                    data.append(result)
                    return data, p_found
        if len(filter_string) != 0:
            data.append(filter_string)
        return data, False

    return parse_helper(iter(filter))


    def process_string(data):
        """
        dummy function used to parse filter strings

        test1 = "( true or false) and (true and true)"
        Evaluate to False
        test2 = "((true and false) and true)"

        moved to test later
        """
        data = data.rstrip(' ').lstrip(' ')
        if 'true' == data:
            return True
        elif 'false' == data:
            return False
        else:
            raise Exception("Unhandled filter substring: {}!".format(data))

    def evaluate_expression(data, process_function=process_string):
        """
        function used to go through a nested list and evaluate the content using a helper process_function

        :param data: a nested list with strings that should be evaluated
        :type data: nested list
        :param process_function: function used to evaluate string
        :type process_function: a function
        :param

        :return: same as process function returns
        """
        if isinstance(data, str):
            return process_function(data)
        elif isinstance(data, list):
            if len(data) == 3:
                return data[1](evaluate_expression(data[0]),evaluate_expression(data[2]))
            elif len(data) == 1:
                return evaluate_expression(data[0])
            else:
                raise Exception("Unhandled data structure case {}, unexpected number of items {}!".format(data, len(data)))
        else:
            print(str(isinstance(data, str)) + " " + str(len(data)) + " " + str(data))
            raise Exception("Unhandled data structure case {}!".format(data))




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
