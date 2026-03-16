from __future__ import annotations

import unittest

import numpy as np

from mf12_python.compare import compare_fields, tolerances_pass


class CompareTests(unittest.TestCase):
    def test_compare_fields_zero_error(self) -> None:
        ref = np.array([[1.0, 2.0], [3.0, 4.0]])
        metrics = compare_fields(ref, ref, "eta")
        self.assertEqual(metrics["eta_max_abs_err"], 0.0)
        self.assertEqual(metrics["eta_rms_err"], 0.0)
        self.assertEqual(metrics["eta_relative_l2_err"], 0.0)

    def test_tolerances(self) -> None:
        self.assertTrue(tolerances_pass({"eta_max_abs_err": 1e-10}, {"eta_max_abs_err": 1e-9}))
        self.assertFalse(tolerances_pass({"eta_max_abs_err": 1e-7}, {"eta_max_abs_err": 1e-9}))


if __name__ == "__main__":
    unittest.main()
