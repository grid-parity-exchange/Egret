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
        with self.assertLogs('egret.model_library.transmission.tx_utils', level=logging.WARNING) as cm:
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=5,
                                                                    p_max=90,
                                                                    gen_name='foo',
                                                                    t=None)
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(cleaned_values, curve['values'])
        self.assertEqual(cm.output, ['WARNING:egret.model_library.transmission.tx_utils:WARNING: Extending piecewise linear cost curve beyond p_min and/or p_max for generator foo (and perhaps others)'])
        # reset for next test
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pw_high_p_max(self):
        curve = example_pw_curve()
        expected_values = copy.deepcopy(curve['values'])
        expected_values.pop(-1)
        expected_values.append((95, 543))
        with self.assertLogs('egret.model_library.transmission.tx_utils', level=logging.WARNING) as cm:
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=10,
                                                                    p_max=95,
                                                                    gen_name='foo',
                                                                    t=None)
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(cleaned_values, curve['values'])
        self.assertEqual(cm.output, ['WARNING:egret.model_library.transmission.tx_utils:WARNING: Extending piecewise linear cost curve beyond p_min and/or p_max for generator foo (and perhaps others)'])
        # reset for next test
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pw_high_p_max_low_p_min_debug(self):
        curve = example_pw_curve()
        # evoke warning once
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=10,
                                                                p_max=95,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = example_pw_curve()['values']
        expected_values.pop(-1)
        expected_values.append((95, 543))
        expected_values.pop(0)
        expected_values.insert(0, (5, 3))
        with self.assertLogs('egret.model_library.transmission.tx_utils', level=logging.DEBUG) as cm:
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=5,
                                                                    p_max=95,
                                                                    gen_name='foo',
                                                                    t=None)
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(cleaned_values, curve['values'])
        # debug output this time
        self.assertEqual(cm.output, ['DEBUG:egret.model_library.transmission.tx_utils:WARNING: Extending piecewise linear cost curve beyond p_min and/or p_max for generator foo'])
        # reset for next test
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

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

    def test_extra_pw_pieces_below_pmin2(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=85,
                                                                p_max=90,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(85, 453), (90, 498)]
        self.assertEqual(cleaned_values, expected_values)
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

    def test_extra_pw_pieces_above_pmax2(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=10,
                                                                p_max=15,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(10,18), (15,33)]
        self.assertEqual(cleaned_values, expected_values)
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

    def test_pw_repeated_slope2(self):
        curve = dict()
        curve['data_type'] = 'cost_curve'
        curve['cost_curve_type'] = 'piecewise'
        curve['values'] = [(10, 20),
                           (30, 40),
                           (50, 60),
                           (70, 90),
                           (90, 120)]
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=5,
                                                                p_max=100,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(5, 15), (50,60), (100, 135)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pw_pmin_is_pmax_on_curve(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=40,
                                                                p_max=40,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(40, 128)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])

    def test_pw_pmin_is_pmax_on_curve_2(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=50,
                                                                p_max=50,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(50, 178.0)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])

    def test_pw_pmin_is_pmax_single_point(self):
        curve = dict()
        curve['data_type'] = 'cost_curve'
        curve['cost_curve_type'] = 'piecewise'
        curve['values'] = [(40,128)]
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=40,
                                                                p_max=40,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(40, 128)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])

    def test_pmax_less_than_first_point(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=-5,
                                                                p_max=5,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(-5, -27), (5, 3)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmax_less_than_first_point2(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=5,
                                                                p_max=5,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(5, 3)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmax_less_than_first_point3(self):
        curve = example_pw_curve()
        curve['values'] = 2*[curve['values'][0]] + curve['values']
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=-5,
                                                                p_max=5,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(-5, -27), (5, 3)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmax_is_first_point(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=0,
                                                                p_max=10,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(0, -12), (10, 18)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmax_is_first_point2(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=10,
                                                                p_max=10,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(10, 18)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmin_greater_than_last_point(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=100,
                                                                p_max=110,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(100, 588), (110, 678)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmin_greater_than_last_point2(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=100,
                                                                p_max=100,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(100, 588)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmin_greater_than_last_point3(self):
        curve = example_pw_curve()
        curve['values'].append(curve['values'][-1])
        curve['values'].append(curve['values'][-1])
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=100,
                                                                p_max=110,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(100, 588), (110, 678)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmin_is_last_point(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=90,
                                                                p_max=110,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(90, 498), (110, 678)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pmin_is_last_point2(self):
        curve = example_pw_curve()
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=90,
                                                                p_max=90,
                                                                gen_name='foo',
                                                                t=None)
        expected_values = [(90, 498)]
        self.assertEqual(cleaned_values, expected_values)
        self.assertIsNot(expected_values, curve['values'])
        tx_utils.validate_and_clean_cost_curve._printed_warning = False

    def test_pw_single_point_raises_value_error(self):
        curve = example_pw_curve()
        curve = dict()
        curve['data_type'] = 'cost_curve'
        curve['cost_curve_type'] = 'piecewise'
        curve['values'] = [(40,128)]
        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=30,
                                                                    p_max=30,
                                                                    gen_name='foo',
                                                                    t=None)

        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=30,
                                                                    p_max=40,
                                                                    gen_name='foo',
                                                                    t=None)
        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=30,
                                                                    p_max=50,
                                                                    gen_name='foo',
                                                                    t=None)

        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=40,
                                                                    p_max=50,
                                                                    gen_name='foo',
                                                                    t=None)

        with self.assertRaises(ValueError):
            cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                    curve_type='cost_curve',
                                                                    p_min=45,
                                                                    p_max=50,
                                                                    gen_name='foo',
                                                                    t=None)

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

class TestTxUtils(unittest.TestCase):

    def test_element_types(self):
        etypes = list(tx_utils.element_types())
        self.assertNotEqual(len(etypes), 0)
        self.assertEqual(len(etypes), len(set(etypes)))
