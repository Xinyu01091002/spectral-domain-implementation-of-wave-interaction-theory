from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from .compare import compare_result_to_reference, tolerances_pass
from .io import load_case, save_result_bundle
from .spectral import run_case


def main() -> None:
    parser = argparse.ArgumentParser(description="MF12 Python spectral-only tools")
    sub = parser.add_subparsers(dest="command", required=True)

    verify = sub.add_parser("verify", help="Run Python reconstruction against MATLAB reference outputs")
    verify.add_argument("case_dir")
    verify.add_argument("--repeats", type=int, default=1)
    verify.add_argument("--warmup", action="store_true")
    verify.add_argument("--output-dir", default=None)

    bench = sub.add_parser("benchmark", help="Benchmark one or more case directories")
    bench.add_argument("case_dirs", nargs="+")
    bench.add_argument("--repeats", type=int, default=3)
    bench.add_argument("--warmup", action="store_true")
    bench.add_argument("--output-dir", default="outputs/cross_language_comparison/benchmarks")

    args = parser.parse_args()
    if args.command == "verify":
        verify_case(args.case_dir, args.repeats, args.warmup, args.output_dir)
    else:
        benchmark_cases(args.case_dirs, args.repeats, args.warmup, args.output_dir)


def verify_case(case_dir: str, repeats: int, warmup: bool, output_dir: str | None) -> None:
    case = load_case(case_dir)
    result = run_case(case, repeats=repeats, warmup=warmup)
    result["case_id"] = case["case_id"]
    result["comparison"] = compare_result_to_reference(result, case)
    if output_dir is None:
        output_dir = str(Path("outputs/cross_language_comparison/verify") / f"{case['case_id']}_{timestamp()}")
    save_result_bundle(output_dir, result)
    summary = {
        "case_id": case["case_id"],
        "output_dir": str(Path(output_dir).resolve()),
        "pass": tolerances_pass(result["comparison"], case.get("tolerances", {})),
        "metrics": result["comparison"],
        "tolerances": case.get("tolerances", {}),
    }
    print(json.dumps(summary, indent=2))
    if not summary["pass"]:
        raise SystemExit(1)


def benchmark_cases(case_dirs: list[str], repeats: int, warmup: bool, output_dir: str) -> None:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for case_dir in case_dirs:
        case = load_case(case_dir)
        result = run_case(case, repeats=repeats, warmup=warmup)
        metrics = compare_result_to_reference(result, case) if case.get("reference") else {}
        rows.append(
            {
                "case_id": case["case_id"],
                "mean_coefficient_s": result["runtime"]["mean_coefficient_s"],
                "mean_reconstruction_s": result["runtime"]["mean_reconstruction_s"],
                "mean_kinematics_s": result["runtime"].get("mean_kinematics_s", ""),
                "mean_total_s": result["runtime"]["mean_total_s"],
                "best_total_s": result["runtime"]["best_total_s"],
                "eta_max_abs_err": metrics.get("eta_max_abs_err", ""),
                "phi_max_abs_err": metrics.get("phi_max_abs_err", ""),
                "u_max_abs_err": metrics.get("u_max_abs_err", ""),
                "v_max_abs_err": metrics.get("v_max_abs_err", ""),
                "w_max_abs_err": metrics.get("w_max_abs_err", ""),
                "p_max_abs_err": metrics.get("p_max_abs_err", ""),
                "phi_vol_max_abs_err": metrics.get("phi_vol_max_abs_err", ""),
                "uV_max_abs_err": metrics.get("uV_max_abs_err", ""),
                "vV_max_abs_err": metrics.get("vV_max_abs_err", ""),
                "a_x_max_abs_err": metrics.get("a_x_max_abs_err", ""),
                "a_y_max_abs_err": metrics.get("a_y_max_abs_err", ""),
                "speedup_vs_matlab_total": metrics.get("speedup_vs_matlab_total", ""),
                "speedup_vs_matlab_reconstruction": metrics.get("speedup_vs_matlab_reconstruction", ""),
            }
        )

    stamp = timestamp()
    csv_path = out_dir / f"python_benchmark_{stamp}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    json_path = out_dir / f"python_benchmark_{stamp}.json"
    json_path.write_text(json.dumps({"generated_at": timestamp(), "rows": rows}, indent=2), encoding="utf-8")
    print(json.dumps({"csv": str(csv_path.resolve()), "json": str(json_path.resolve()), "rows": rows}, indent=2))


def timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


if __name__ == "__main__":
    main()
