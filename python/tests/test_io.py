from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from mf12_python.io import load_case


class IoTests(unittest.TestCase):
    def test_load_case_reads_arrays(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "inputs").mkdir()
            np.savetxt(root / "inputs" / "a.csv", np.array([1.0, 2.0]), delimiter=",")
            np.savetxt(root / "inputs" / "b.csv", np.array([3.0, 4.0]), delimiter=",")
            np.savetxt(root / "inputs" / "kx.csv", np.array([0.1, 0.2]), delimiter=",")
            np.savetxt(root / "inputs" / "ky.csv", np.array([0.0, 0.1]), delimiter=",")
            manifest = {
                "case_id": "tmp_case",
                "inputs": {"order": 2, "g": 9.81, "h": 40.0, "Ux": 0.0, "Uy": 0.0, "Lx": 10.0, "Ly": 10.0, "Nx": 4, "Ny": 4, "t": 0.0},
                "arrays": {"a": "inputs/a.csv", "b": "inputs/b.csv", "kx": "inputs/kx.csv", "ky": "inputs/ky.csv"},
            }
            (root / "case.json").write_text(json.dumps(manifest), encoding="utf-8")
            loaded = load_case(root)
            self.assertEqual(loaded["case_id"], "tmp_case")
            self.assertEqual(loaded["arrays"]["a"].shape, (2,))
            self.assertAlmostEqual(float(loaded["arrays"]["ky"][1]), 0.1)


if __name__ == "__main__":
    unittest.main()
