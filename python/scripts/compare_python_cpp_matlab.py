from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from mf12_python.compare import compare_fields
from mf12_python.io import load_case
from mf12_python.spectral import run_case


def load_cpp_result(case_id: str, root: Path) -> dict[str, np.ndarray | dict]:
    result_dir = root / case_id
    return {
        "eta": np.loadtxt(result_dir / "eta.csv", delimiter=","),
        "phi": np.loadtxt(result_dir / "phi.csv", delimiter=","),
        "x": np.loadtxt(result_dir / "x.csv", delimiter=",").reshape(-1),
        "y": np.loadtxt(result_dir / "y.csv", delimiter=",").reshape(-1),
        "meta": json.loads((result_dir / "result.json").read_text(encoding="utf-8")),
    }


def max_abs_between(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare MATLAB, Python, and C++ MF12 shared cases")
    parser.add_argument("case_dirs", nargs="+", help="Case directories under cross_language_comparison/cases/")
    parser.add_argument("--cpp-root", default="outputs/cross_language_comparison/verify_cpp")
    parser.add_argument("--output-dir", default="outputs/cross_language_comparison/comparisons")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cpp_root = Path(args.cpp_root)

    rows: list[dict[str, object]] = []
    for case_dir in args.case_dirs:
        case = load_case(case_dir)
        ref = case["reference"]["loaded_arrays"]
        py = run_case(case, repeats=1, warmup=False)
        cpp = load_cpp_result(case["case_id"], cpp_root)

        py_eta = compare_fields(py["eta"], ref["eta"], "eta")
        py_phi = compare_fields(py["phi"], ref["phi"], "phi")
        cpp_eta = compare_fields(cpp["eta"], ref["eta"], "eta")
        cpp_phi = compare_fields(cpp["phi"], ref["phi"], "phi")

        row = {
            "case_id": case["case_id"],
            "python_eta_max_abs_err": py_eta["eta_max_abs_err"],
            "python_phi_max_abs_err": py_phi["phi_max_abs_err"],
            "cpp_eta_max_abs_err": cpp_eta["eta_max_abs_err"],
            "cpp_phi_max_abs_err": cpp_phi["phi_max_abs_err"],
            "python_cpp_eta_max_abs_diff": max_abs_between(py["eta"], cpp["eta"]),
            "python_cpp_phi_max_abs_diff": max_abs_between(py["phi"], cpp["phi"]),
            "python_total_s": py["runtime"]["mean_total_s"],
            "cpp_total_s": cpp["meta"]["runtime"]["mean_total_s"],
        }
        rows.append(row)

    csv_path = out_dir / "python_cpp_matlab_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    json_path = out_dir / "python_cpp_matlab_summary.json"
    json_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(json.dumps({"csv": str(csv_path.resolve()), "json": str(json_path.resolve()), "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
