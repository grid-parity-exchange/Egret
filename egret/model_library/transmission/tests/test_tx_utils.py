import unittest
from egret.model_library.transmission import tx_utils
import logging
import copy


def _example_quadratic(p):
    return 0.05 * p ** 2 + p + 3


def example_pw_curve():
    curve = dict()
    curve['data_type'] = 'cost_curve'
    curve['cost_curve_type'] = 'piecewise'
    curve['values'] = [(10, _example_quadratic(10)),
                       (30, _example_quadratic(30)),
                       (50, _example_quadratic(50)),
                       (70, _example_quadratic(70)),
                       (90, _example_quadratic(90))]
    return curve


def example_poly_curve():
    curve = dict()
    curve['data_type'] = 'cost_curve'
    curve['cost_curve_type'] = 'polynomial'
    curve['values'] = {0: 3, 1: 1, 2: 0.05}
    return curve


class TestValidateCostCurves(unittest.TestCase):
    def test_pw_simple(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=10,
                                                                p_max=90,
                                                                gen_name='foo',
                                                                t=None)
        self.assertEqual(cleaned_values, curve['values'])
        self.assertIsNot(cleaned_values, curve['values'])

    def test_wrong_curve_type(self):
        curve = example_pw_curve()
        curve['data_type'] = 'blah'
        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=10,
                                                                    p_max=90,
                                                                    gen_name='foo',
                                                                    t=None)

    def test_pw_no_values(self):
        curve = example_pw_curve()
        curve['values'] = list()
        with self.assertLogs('egret.model_library.transmission.tx_utils', level=logging.WARNING) as cm:
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=10,
                                                                    p_max=90,
                                                                    gen_name='foo',
                                                                    t=None)
        self.assertEqual(cleaned_values, curve['values'])
        self.assertIsNot(cleaned_values, curve['values'])
        self.assertEqual(cm.output, ['WARNING:egret.model_library.transmission.tx_utils:WARNING: Generator foo has no cost information associated with it'])

    def test_pw_repeat_value_and_cost(self):
        curve = example_pw_curve()
        orig_values = copy.deepcopy(curve['values'])
        curve['values'].insert(2, (30, _example_quadratic(30)))
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=10,
                                                                p_max=90,
                                                                gen_name='foo',
                                                                t=None)
        self.assertEqual(cleaned_values, orig_values)

    def test_pw_repeat_value(self):
        curve = example_pw_curve()
        curve['values'].insert(2, (30, _example_quadratic(40)))
        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=10,
                                                                    p_max=90,
                                                                    gen_name='foo',
                                                                    t=None)

    def test_pw_nonconvex(self):
        curve = example_pw_curve()
        curve['values'].insert(2, (35, _example_quadratic(20)))
        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=10,
                                                                    p_max=90,
                                                                    gen_name='foo',
                                                                    t=None)

    def test_pw_low_p_min(self):
        curve = example_pw_curve()
        expected_values = copy.deepcopy(curve['values'])
        expected_values.pop(0)
        expected_values.insert(0, (5, 3))
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=5,
                                                                p_max=90,
                                                                gen_name='foo',
                                                                t=None)
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(cleaned_values, curve['values'])

    def test_pw_high_p_max(self):
        curve = example_pw_curve()
        expected_values = copy.deepcopy(curve['values'])
        expected_values.pop(-1)
        expected_values.append((95, 543))
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=10,
                                                                p_max=95,
                                                                gen_name='foo',
                                                                t=None)
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(cleaned_values, curve['values'])

    def test_extra_pw_pieces_below_pmin(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=30,
                                                                p_max=90,
                                                                gen_name='foo',
                                                                t=None)
        self.assertEqual(cleaned_values, curve['values'][1:])
        self.assertIsNot(cleaned_values, curve['values'])

    def test_extra_pw_pieces_above_pmax(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=10,
                                                                p_max=70,
                                                                gen_name='foo',
                                                                t=None)
        self.assertEqual(cleaned_values, curve['values'][:-1])
        self.assertIsNot(cleaned_values, curve['values'])

    def test_pw_repeated_slope(self):
        curve = dict()
        curve['data_type'] = 'cost_curve'
        curve['cost_curve_type'] = 'piecewise'
        curve['values'] = [(10, 20),
                           (30, 40),
                           (50, 60),
                           (70, 80),
                           (90, 100)]
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=10,
                                                                p_max=90,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(10, 20), (90, 100)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])

    def test_poly_simple(self):
        curve = example_poly_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=15,
                                                                p_max=85,
                                                                gen_name='foo',
                                                                t=None)
        self.assertEqual(cleaned_values, curve['values'])
        self.assertIsNot(cleaned_values, curve['values'])

    def test_poly_nonconvex(self):
        curve = example_poly_curve()
        curve['values'][2] = -1
        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=15,
                                                                    p_max=85,
                                                                    gen_name='foo',
                                                                    t=None)

    def test_poly_cubic(self):
        curve = example_poly_curve()
        curve['values'][3] = 1
        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=15,
                                                                    p_max=85,
                                                                    gen_name='foo',
                                                                    t=None)
