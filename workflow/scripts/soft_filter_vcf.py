import logging
import re
import yaml

from collections import OrderedDict
from hydra_genetics.utils.io import utils
from pysam import VariantFile




def is_float(element) -> bool:
    try:
        float(element)
        return True
    except ValueError:
        return False

_and_function = lambda v1, v2: v1 and v2
_or_function = lambda v1, v2: v1 or v2

def _parse_helper(f_iterator):
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
            result, closing_parantes = _parse_helper(f_iterator)
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
                data.append(_or_function)
                result, p_found = _parse_helper(f_iterator)
                if len(result) == 1:
                    result = result[0]
                data.append(result)
                return data, p_found
            elif filter_string.startswith(' and '):
                data.append(_and_function)
                result, p_found = _parse_helper(f_iterator)
                if len(result) == 1:
                    result = result[0]
                data.append(result)
                return data, p_found
            elif filter_string.endswith(' or '):
                data.append(filter_string[:-4])
                data.append(_or_function)
                result, p_found = _parse_helper(f_iterator)
                if len(result) == 1:
                    result = result[0]
                data.append(result)
                return data, p_found
            elif filter_string.endswith(' and '):
                data.append(filter_string[:-5])
                data.append(_and_function)
                result, p_found = _parse_helper(f_iterator)
                if len(result) == 1:
                    result = result[0]
                data.append(result)
                return data, p_found
    if len(filter_string) != 0:
        data.append(filter_string)
    return data, False

def _convert_string(data, process_function):
    """
    converts a nested list of string into functions
    :params data: data structure with string
    :type list: nested list with strings
    :params process_function: function used to convert string to lambda function
    :type list: function

    :return: nested list with lambda functions
    """
    if isinstance(data, str):
        return process_function(data)
    elif isinstance(data, list):
        if len(data) == 3:
            return [_convert_string(data[0], process_function), data[1], _convert_string(data[2], process_function)]
        elif len(data) == 1:
            return _convert_string(data[0], process_function)
        else:
            raise Exception("Unhandled data structure case {}, unexpected number of items {}!".format(data, len(data)))
    else:
        raise Exception("Unhandled data structure case {}\ntype: {}\nlength: {}".format(data, type(data), len(data)))

def _evaluate_expression(data, variant):
    """
    function used to evaluate variant using a nested list of functions

    :param data: a nested list with functions that will be evaluated
    :type data: nested list
    :param variant: variant that will be evaluated
    :type process_function: a function
    :param

    :return: boolean
    """
    if callable(data):
        return data(variant)
    elif isinstance(data, list):
        if len(data) == 3:
            return data[1](_evaluate_expression(data[0], variant), _evaluate_expression(data[2], variant))
        elif len(data) == 1:
            return _evaluate_expression(data[0], variant)
        else:
            raise Exception("Unhandled data structure case {}, unexpected number of items {}!".format(data, len(data)))
    else:
        raise Exception("Unhandled data structure case {}!".format(data))

def create_variant_filter(filter, string_parser):
    data, _ = _parse_helper(iter(filter))
    data = _convert_string(data, string_parser)
    return lambda variant: _evaluate_expression(data, variant)


def create_convert_expresion_function(annotation_extractors):
    """

    :param annotation_extractors: dict with functions that can extract annotatiors
    :type annotation_extractors: dict
    """
    def compare_data(comparison, value1, value2, index1=None, index2=None):
        if isinstance(value1, float):
            if value2 == "" or value2 is None:
                return False
            if index2 is not None:
                value2 = value2[index2]
            return comparison(value1, float(value2))
        elif isinstance(value2, float):
            if value1 == "" or value1 is None:
                return False
            if index1 is not None:
                value1 = value1[index1]
            return comparison(float(value1), value2)
        else:
            if index2 is not None:
                value2 = value2[index2]
            if index1 is not None:
                value1 = value2[index1]
            return comparison(value1, value2)

    def convert_to_expresion(expression):
        """
        :params expression: VEP:AF
        :type expression: string
        """
        comparison = {
            ">": lambda value1, value2: value1 > value2,
            "<": lambda value1, value2: value1 < value2,
            "=": lambda value1, value2: value1 == value2,
            "!=": lambda value1, value2: value1 != value2
        }
        data = re.split("[ ]([<>=!]+)[ ]", expression)
        if len(data) != 3:
            raise Exception("Invalid expression: " + expression)
        regex_string = "[ ]*(VEP|FORMAT|INFO):([A-Za-z0-9_.]+):*([0-9]*)"

        if "VEP:" in data[0] or "FORMAT:" in data[0] or  "INFO:" in data[0]:
            source, field, index = re.search(regex_string, data[0]).groups()
            if len(index) == 0:
                index = None
            else:
                index = int(index)
            try:
                data[2] = data[2].rstrip(" ").lstrip(" ")
                value2 = float(data[2])
                return lambda variant: compare_data(comparison[data[1]], annotation_extractors[source](variant, field), value2, index1=index)
            except ValueError:
                return lambda variant: compare_data(comparison[data[1]], annotation_extractors[source](variant, field), data[2], index1=index)
        elif "VEP:" in data[2] or "FORMAT:" in data[2] or  "INFO:" in data[2]:
            source, field, index = re.search(regex_string, data[2]).groups()
            if len(index) == 0:
                index = None
            else:
                index = int(index)
            try:
                data[0] = data[0].rstrip(" ").lstrip(" ")
                value1 = float(data[0])
                return lambda variant: compare_data(comparison[data[1]], value1, annotation_extractors[source](variant, field), index2=index)
            except ValueError:
                return lambda variant: compare_data(comparison[data[1]], data[0], annotation_extractors[source](variant, field), index2=index)
        else:
            raise Exception("Could not find comparison field in: " + expression)

    return convert_to_expresion

# ToDo Move to tools ? and maybe check correct sample?
def extract_format_data(variant, field):
    return variant.samples[0][field]

if __name__ == "__main__":
    in_vcf = snakemake.input.vcf
    out_vcf = snakemake.output.vcf
    filter_yaml_file = snakemake.params.filter_config

    variants = VariantFile(in_vcf)
    log = logging.getLogger()

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
