import unittest
import os
import pytest
import tempfile
import yaml

from hydra_genetics.utils.io import utils
from pysam import VariantFile

from soft_filter_vcf import create_variant_filter
from soft_filter_vcf import _parse_helper
from soft_filter_vcf import _convert_string
from soft_filter_vcf import _and_function
from soft_filter_vcf import _or_function
from soft_filter_vcf import create_convert_expression_function
from soft_filter_vcf import extract_format_data
from soft_filter_vcf import extract_info_data
from soft_filter_vcf import soft_filter_variants


class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        self.in_vcf = ".tests/integration/snv_indels/ensemble_vcf/HD832.HES45_T.ensembled.vep_annotated.vcf.gz"

    def tearDown(self):
        pass

    def test_data_structur_generation(self):
        """
        function used to convert filter string into a function that can evaluate objects
        """
        def true_fun(variant): return True
        def false_fun(variant): return False

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

    def test_create_convert_expression(self):
        variants = VariantFile(self.in_vcf)
        annotation_extractor = {}
        for record in variants.header.records:
            if record.type == "INFO":
                if record['ID'] == "CSQ":
                    vep_fields = {v: c for c, v in enumerate(record['Description'].split("Format: ")[1].split("|"))}
                    annotation_extractor["VEP"] = utils.get_annoation_data_vep(vep_fields)
        annotation_extractor['FORMAT'] = extract_format_data
        annotation_extractor['INFO'] = extract_info_data
        expression_converter = create_convert_expression_function(annotation_extractor)

        def creater_filter(yaml_file):
            filters = {"filters": []}
            with open(yaml_file) as file:
                filters = yaml.load(file, Loader=yaml.FullLoader)
            vcf_filters = []
            for filter, value in filters["filters"].items():
                vcf_filters.append(create_variant_filter(value['expression'], expression_converter))
            return vcf_filters

        def test_filters(test_table, variants, vcf_filters):
            for variant in variants:
                for vcf_filter in vcf_filters:
                    try:
                        self.assertEqual(test_table["{}:{}-{}".format(variant.chrom, variant.start, variant.stop)],
                                         vcf_filter(variant))
                    except AssertionError as e:
                        print("Failed validation gene: " + str(variant))
                        raise e

        test_table = {
                        "chr1:934486-934487": True,  # SAMD11
                        "chr1:935221-935222": True,  # SAMD11
                        "chr1:935338-935339": True,  # SAMD11
                        "chr1:2460943-2460947": False,  # PLCH2
                        "chr1:2460960-2460961": False,  # PLCH2
                        "chr1:2461206-2461207": False  # PLCH2
                      }

        test_filters(test_table,
                     VariantFile(self.in_vcf),
                     creater_filter(".tests/integration/config_soft_filter_unittest_1.yaml"))
        test_table = {
                        "chr1:934486-934487": False,  # SAMD11
                        "chr1:935221-935222": False,  # SAMD11
                        "chr1:935338-935339": False,  # SAMD11
                        "chr1:2460943-2460947": False,  # PLCH2 -
                        "chr1:2460960-2460961": True,  # PLCH2 T
                        "chr1:2461206-2461207": False  # PLCH2 A
                      }

        test_filters(test_table,
                     VariantFile(self.in_vcf),
                     creater_filter(".tests/integration/config_soft_filter_unittest_2.yaml"))

        test_table = {
                        "chr1:934486-934487": False,  # SAMD11
                        "chr1:935221-935222": False,  # SAMD11
                        "chr1:935338-935339": False,  # SAMD11
                        "chr1:2460943-2460947": False,  # PLCH2 -
                        "chr1:2460960-2460961": True,  # PLCH2 T
                        "chr1:2461206-2461207": True  # PLCH2 A
                      }

        test_filters(test_table,
                     VariantFile(self.in_vcf),
                     creater_filter(".tests/integration/config_soft_filter_unittest_3.yaml"))

        test_table = {
                        "chr1:934486-934487": True,  # 1108,595
                        "chr1:935221-935222": True,  # 936,910
                        "chr1:935338-935339": False,  # 1301,127
                        "chr1:2460943-2460947": False,  # 278,9
                        "chr1:2460960-2460961": False,  # 314,5
                        "chr1:2461206-2461207": False  # 657,59
                      }

        test_filters(test_table,
                     VariantFile(self.in_vcf),
                     creater_filter(".tests/integration/config_soft_filter_unittest_4.yaml"))

        test_table = {
                        "chr1:934486-934487": False,  # SAMD11 1108,595
                        "chr1:935221-935222": False,  # SAMD11 936,910
                        "chr1:935338-935339": False,  # SAMD11 1301,127
                        "chr1:2460943-2460947": False,  # PLCH2 278,9
                        "chr1:2460960-2460961": False,  # PLCH2 314,5
                        "chr1:2461206-2461207": True  # PLCH2 657,59
                      }

        test_filters(test_table,
                     VariantFile(self.in_vcf),
                     creater_filter(".tests/integration/config_soft_filter_unittest_5.yaml"))

        test_table = {
                        "chr1:934486-934487": False,  # -
                        "chr1:935221-935222": False,  # -
                        "chr1:935338-935339": False,  # -
                        "chr1:2460943-2460947": True,  # XM_017002872
                        "chr1:2460960-2460961": True,  # XM_017002872
                        "chr1:2461206-2461207": True  # XM_017002872
                      }

        test_filters(test_table,
                     VariantFile(self.in_vcf),
                     creater_filter(".tests/integration/config_soft_filter_unittest_8.yaml"))

        test_table = {
                        "chr1:934486-934487": True,  # -
                        "chr1:935221-935222": True,  # -
                        "chr1:935338-935339": True,  # -
                        "chr1:2460943-2460947": False,  # XM_017002872
                        "chr1:2460960-2460961": False,  # XM_017002872
                        "chr1:2461206-2461207": False  # XM_017002872
                      }

        test_filters(test_table,
                     VariantFile(self.in_vcf),
                     creater_filter(".tests/integration/config_soft_filter_unittest_9.yaml"))

        test_table = {
                        "chr1:934486-934487": False,  # -
                        "chr1:935221-935222": False,  # -
                        "chr1:935338-935339": False,  # -
                        "chr1:2460943-2460947": False,  # XM_017002872
                        "chr1:2460960-2460961": False,  # XM_017002872 mutect2
                        "chr1:2461206-2461207": True  # XM_017002872
                      }

        test_filters(test_table,
                     VariantFile(self.in_vcf),
                     creater_filter(".tests/integration/config_soft_filter_unittest_10.yaml"))

    def test_soft_and_hard_filtering(self):
        tempdir = tempfile.mkdtemp()
        vcf = os.path.join(tempdir, "test.vcf")
        with open(vcf, 'w', encoding="ascii") as out_vcf:
            soft_filter_variants(self.in_vcf, out_vcf, ".tests/integration/config_soft_filter_unittest_6.yaml")
        with VariantFile(vcf) as variants:
            test_table = {
                            "chr1:934486-934487": "SAMD11",  # SAMD11 1108,595
                            "chr1:935221-935222": "SAMD11",  # SAMD11 936,910
                            "chr1:935338-935339": "SAMD11",  # SAMD11 1301,127
                            "chr1:2460943-2460947": "AD50",  # PLCH2 278,9
                            "chr1:2460960-2460961": "AD50",  # PLCH2 314,5
                            "chr1:2461206-2461207": "PASS"  # PLCH2 657,59
                          }
            counter = 0
            for variant in variants:
                try:
                    counter += 1
                    self.assertEqual(test_table["{}:{}-{}".format(variant.chrom, variant.start, variant.stop)],
                                     ",".join([f for f in variant.filter]))
                except AssertionError as e:
                    print("Failed validation gene: " + str(variant))
                    raise e
            self.assertEqual(counter, 6)

        with open(vcf, 'w', encoding="ascii") as out_vcf:
            soft_filter_variants(self.in_vcf, out_vcf, ".tests/integration/config_soft_filter_unittest_7.yaml")
        with VariantFile(vcf) as variants:
            test_table = {
                            "chr1:934486-934487": "SAMD11",  # SAMD11 1108,595
                            "chr1:935221-935222": "SAMD11",  # SAMD11 936,910
                            "chr1:935338-935339": "SAMD11",  # SAMD11 1301,127
                            "chr1:2461206-2461207": "PASS"  # PLCH2 657,59
                          }
            counter = 0
            for variant in variants:
                try:
                    counter += 1
                    self.assertEqual(test_table["{}:{}-{}".format(variant.chrom, variant.start, variant.stop)],
                                     ",".join([f for f in variant.filter]))
                except AssertionError as e:
                    print("Failed validation gene: " + str(variant))
                    raise
            self.assertEqual(counter, 4)
