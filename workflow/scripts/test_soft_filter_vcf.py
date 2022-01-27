import unittest
import pytest
import yaml

from hydra_genetics.utils.io import utils
from pysam import VariantFile

from soft_filter_vcf import create_variant_filter
from soft_filter_vcf import _parse_helper
from soft_filter_vcf import _convert_string
from soft_filter_vcf import _and_function
from soft_filter_vcf import _or_function
from soft_filter_vcf import create_convert_expresion_function
from soft_filter_vcf import extract_format_data

class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_data_structur_generation(self):
        """
        function used to convert filter string into a function that can evaluate objects
        """
        true_fun = lambda variant: True
        false_fun = lambda variant: False
        def process_string(data):
            """
            dummy function used to parse filter strings
            """
            data = data.rstrip(' ').lstrip(' ')
            if 'true' == data:
                return true_fun
            elif 'false' == data:
                return false_fun
            else:
                raise Exception("Unhandled filter substring: {}!".format(data))


        test1 = "(true or false) and (true and true)"
        data, _ = _parse_helper(iter(test1))
        self.assertEqual(data, ([['true', _or_function, 'false'], _and_function, ['true', _and_function, 'true']]))

        data = _convert_string(data, process_string)
        self.assertEqual(data, [[true_fun, _or_function, false_fun], _and_function, [true_fun, _and_function, true_fun]])

        variant_filter = create_variant_filter(test1, process_string)
        self.assertEqual(variant_filter(""), True)

        test2 = "((true and false) and true)"
        data, _ = _parse_helper(iter(test2))
        self.assertEqual(data, ([[['true', _and_function, 'false'], _and_function, 'true']]))

        data = _convert_string(data, process_string)
        self.assertEqual(data, [[true_fun, _and_function, false_fun], _and_function, true_fun])

        variant_filter = create_variant_filter(test2, process_string)
        self.assertEqual(variant_filter(""), False)

    def test_create_convert_expresion(self):
        test_table = {
                        "chr1:934486-934487": True, # SAMD11
                        "chr1:935221-935222": True, # SAMD11
                        "chr1:935338-935339": True, # SAMD11
                        "chr1:2460943-2460947": False, # PLCH2
                        "chr1:2460960-2460961": False, # PLCH2
                        "chr1:2461206-2461207": False # PLCH2
                      }

        variants = VariantFile(".tests/integration/snv_indels/ensemble_vcf/HD832.HES45_T.ensembled.vep_annotated.vcf.gz")
        annotation_extractor = {}
        for record in variants.header.records:
            if record.type == "INFO":
                if record['ID'] == "CSQ":
                    vep_fields = {v: c for c, v in enumerate(record['Description'].split("Format: ")[1].split("|"))}
                    annotation_extractor["VEP"] = utils.get_annoation_data_vep(vep_fields)
        annotation_extractor['FORMAT'] = extract_format_data
        expresion_converter = create_convert_expresion_function(annotation_extractor)
        filters = {"filters": []}
        with open(".tests/integration/config_soft_filter.yaml") as file:
            filters = yaml.load(file, Loader=yaml.FullLoader)
        vcf_filters = []
        for filter, value in filters["filters"].items():
            #data, _ = _parse_helper(iter(value['expression']))
            #vcf_filters.append(_convert_string(data, expresion_converter))
            vcf_filters.append(create_variant_filter(value['expression'], expresion_converter))
        for variant in variants:
            for vcf_filter in vcf_filters:
                try:
                    self.assertEqual(test_table["{}:{}-{}".format(variant.chrom, variant.start, variant.stop)], vcf_filter(variant))
                except AssertionError as e:
                    print("Failed validation gene: " + str(variant))
                    raise e

        test_table = {
                        "chr1:934486-934487": False, # SAMD11
                        "chr1:935221-935222": False, # SAMD11
                        "chr1:935338-935339": False, # SAMD11
                        "chr1:2460943-2460947": False, # PLCH2 -
                        "chr1:2460960-2460961": True, # PLCH2 T
                        "chr1:2461206-2461207": False # PLCH2 A
                      }
        filters = {"filters": []}
        with open(".tests/integration/config_soft_filter_2.yaml") as file:
            filters = yaml.load(file, Loader=yaml.FullLoader)
        vcf_filters = []
        for filter, value in filters["filters"].items():
            # Create filter from expression
            vcf_filters.append(create_variant_filter(value['expression'], expresion_converter))

        variants = VariantFile(".tests/integration/snv_indels/ensemble_vcf/HD832.HES45_T.ensembled.vep_annotated.vcf.gz")
        for variant in variants:
            for vcf_filter in vcf_filters:
                try:
                    self.assertEqual(test_table["{}:{}-{}".format(variant.chrom, variant.start, variant.stop)], vcf_filter(variant))
                except AssertionError as e:
                    print("Failed validation gene: " + str(variant))
                    raise e

        test_table = {
                        "chr1:934486-934487": False, # SAMD11
                        "chr1:935221-935222": False, # SAMD11
                        "chr1:935338-935339": False, # SAMD11
                        "chr1:2460943-2460947": False, # PLCH2 -
                        "chr1:2460960-2460961": True, # PLCH2 T
                        "chr1:2461206-2461207": True # PLCH2 A
                      }
        filters = {"filters": []}
        with open(".tests/integration/config_soft_filter_3.yaml") as file:
            filters = yaml.load(file, Loader=yaml.FullLoader)
        vcf_filters = []
        for filter, value in filters["filters"].items():
            # Create filter from expression
            vcf_filters.append(create_variant_filter(value['expression'], expresion_converter))

        variants = VariantFile(".tests/integration/snv_indels/ensemble_vcf/HD832.HES45_T.ensembled.vep_annotated.vcf.gz")
        for variant in variants:
            for vcf_filter in vcf_filters:
                try:
                    self.assertEqual(test_table["{}:{}-{}".format(variant.chrom, variant.start, variant.stop)], vcf_filter(variant))
                except AssertionError as e:
                    print("Failed validation gene: " + str(variant))
                    raise e

        test_table = {
                        "chr1:934486-934487": True, # 1108,595
                        "chr1:935221-935222": True, # 936,910
                        "chr1:935338-935339": False, # 1301,127
                        "chr1:2460943-2460947": False, # 278,9
                        "chr1:2460960-2460961": False, # 314,5
                        "chr1:2461206-2461207": False # 657,59
                      }
        filters = {"filters": []}
        with open(".tests/integration/config_soft_filter_4.yaml") as file:
            filters = yaml.load(file, Loader=yaml.FullLoader)
        vcf_filters = []
        for filter, value in filters["filters"].items():
            # Create filter from expression
            vcf_filters.append(create_variant_filter(value['expression'], expresion_converter))

        variants = VariantFile(".tests/integration/snv_indels/ensemble_vcf/HD832.HES45_T.ensembled.vep_annotated.vcf.gz")
        for variant in variants:
            for vcf_filter in vcf_filters:
                try:
                    self.assertEqual(test_table["{}:{}-{}".format(variant.chrom, variant.start, variant.stop)], vcf_filter(variant))
                except AssertionError as e:
                    print("Failed validation gene: " + str(variant))
                    raise e

        test_table = {
                        "chr1:934486-934487": False, # SAMD11 1108,595
                        "chr1:935221-935222": False, # SAMD11 936,910
                        "chr1:935338-935339": False, # SAMD11 1301,127
                        "chr1:2460943-2460947": False, # PLCH2 278,9
                        "chr1:2460960-2460961": False, # PLCH2 314,5
                        "chr1:2461206-2461207": True # PLCH2 657,59
                      }
        filters = {"filters": []}
        with open(".tests/integration/config_soft_filter_5.yaml") as file:
            filters = yaml.load(file, Loader=yaml.FullLoader)
        vcf_filters = []
        for filter, value in filters["filters"].items():
            # Create filter from expression
            vcf_filters.append(create_variant_filter(value['expression'], expresion_converter))

        variants = VariantFile(".tests/integration/snv_indels/ensemble_vcf/HD832.HES45_T.ensembled.vep_annotated.vcf.gz")
        for variant in variants:
            for vcf_filter in vcf_filters:
                try:
                    self.assertEqual(test_table["{}:{}-{}".format(variant.chrom, variant.start, variant.stop)], vcf_filter(variant))
                except AssertionError as e:
                    print("Failed validation gene: " + str(variant))
                    raise e
