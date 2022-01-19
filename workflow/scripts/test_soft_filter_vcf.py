import unittest
import pytest


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
        from soft_filter_vcf import create_variant_filter
        from soft_filter_vcf import _parse_helper
        from soft_filter_vcf import _convert_string
        from soft_filter_vcf import _and_function
        from soft_filter_vcf import _or_function


        test1 = "(true or false) and (true and true)"
        data, _ = _parse_helper(iter(test1))
        self.assertEqual(data, ([['true', _or_function, 'false'], _and_function, ['true', _and_function, 'true']]))

        data = _convert_string(data, process_string)
        print(data)
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
