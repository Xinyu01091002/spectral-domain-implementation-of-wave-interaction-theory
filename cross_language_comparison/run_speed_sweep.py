from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from mf12_python.io import load_case  # noqa: E402
from mf12_python.spectral import run_case  # noqa: E402


PALETTE = {"matlab": "#12355b", "python": "#d95d39", "cpp": "#2a9d8f"}
LINESTYLES = {"matlab": "-", "python": "-.", "cpp": "--"}
MARKERS = {"matlab": "o", "python": "s", "cpp": "^"}


@dataclass(frozen=True)
class SweepPoint:
    component_count: int
    grid_size: int


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark MATLAB, Python, and C++ across retained components and square grid sizes.")
    parser.add_argument("--base-case", default="cross_language_comparison/cases/benchmark_medium")
    parser.add_argument("--component-counts", default="20,30,40")
    parser.add_argument("--grid-sizes", default="64,128,256,512,1024")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmup", action="store_true")
    parser.add_argument("--languages", default="matlab,python,cpp")
    parser.add_argument("--output-dir", default="outputs/cross_language_comparison/speed_sweep")
    parser.add_argument("--cpp-exe", default="cpp/build/mf12_cpp.exe")
    parser.add_argument("--matlab-command", default="matlab")
    args = parser.parse_args()

    component_counts = parse_int_list(args.component_counts)
    grid_sizes = parse_int_list(args.grid_sizes)
    requested_languages = [item.strip().lower() for item in args.languages.split(",") if item.strip()]

    base_case = load_case(REPO_ROOT / args.base_case)
    available_components = int(base_case["arrays"]["a"].size)
    max_requested = max(component_counts)
    if max_requested > available_components:
        raise ValueError(
            f"Requested up to {max_requested} retained components, but {args.base_case} only contains {available_components}. "
            "Create or export a denser benchmark case before sweeping that range."
        )

    output_dir = REPO_ROOT / args.output_dir / timestamp()
    cases_dir = output_dir / "generated_cases"
    cases_dir.mkdir(parents=True, exist_ok=True)

    available_languages, skipped = resolve_languages(requested_languages, args)
    points = [SweepPoint(component_count=n, grid_size=grid) for grid in grid_sizes for n in component_counts]

    rows: list[dict[str, Any]] = []
    for point in points:
        case_dir = write_sweep_case(base_case, point, cases_dir)
        for language in available_languages:
            summary = run_language(language, case_dir, args)
            rows.append(
                {
                    "language": language,
                    "case_id": summary["case_id"],
                    "component_count": int(summary["component_count"]),
                    "Lx": float(summary["domain"]["Lx"]),
                    "Ly": float(summary["domain"]["Ly"]),
                    "Nx": int(summary["grid"]["Nx"]),
                    "Ny": int(summary["grid"]["Ny"]),
                    "mean_coefficient_s": float(summary["runtime"]["mean_coefficient_s"]),
                    "mean_reconstruction_s": float(summary["runtime"]["mean_reconstruction_s"]),
                    "mean_total_s": float(summary["runtime"]["mean_total_s"]),
                    "best_total_s": float(summary["runtime"]["best_total_s"]),
                    "repeats": int(summary["runtime"]["repeats"]),
                }
            )

    csv_path = output_dir / "speed_sweep_summary.csv"
    json_path = output_dir / "speed_sweep_summary.json"
    write_rows(csv_path, rows)
    json_path.write_text(
        json.dumps(
            {
                "generated_at": timestamp(),
                "base_case": str((REPO_ROOT / args.base_case).resolve()),
                "component_counts": component_counts,
                "grid_sizes": grid_sizes,
                "fixed_domain": {"Lx": float(base_case["inputs"]["Lx"]), "Ly": float(base_case["inputs"]["Ly"])},
                "languages": available_languages,
                "skipped_languages": skipped,
                "assumption": "Retained components are selected by descending sqrt(a^2 + b^2) from the base case. The physical domain is kept fixed while Nx = Ny is swept.",
                "rows": rows,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    figure_paths = make_figures(rows, component_counts, grid_sizes, available_languages, output_dir)
    print(
        json.dumps(
            {
                "output_dir": str(output_dir.resolve()),
                "csv": str(csv_path.resolve()),
                "json": str(json_path.resolve()),
                "figures": {key: str(path.resolve()) for key, path in figure_paths.items()},
                "languages": available_languages,
                "skipped_languages": skipped,
                "available_components_in_base_case": available_components,
            },
            indent=2,
        )
    )


def parse_int_list(text: str) -> list[int]:
    values = sorted({int(item.strip()) for item in text.split(",") if item.strip()})
    if not values:
        raise ValueError("Expected at least one integer value.")
    return values


def resolve_languages(requested_languages: list[str], args: argparse.Namespace) -> tuple[list[str], list[dict[str, str]]]:
    available: list[str] = []
    skipped: list[dict[str, str]] = []
    for language in requested_languages:
        if language == "python":
            available.append(language)
            continue
        if language == "cpp":
            cpp_exe = REPO_ROOT / args.cpp_exe
            if cpp_exe.exists():
                available.append(language)
            else:
                skipped.append({"language": language, "reason": f"Missing executable: {cpp_exe}"})
            continue
        if language == "matlab":
            if shutil.which(args.matlab_command):
                available.append(language)
            else:
                skipped.append({"language": language, "reason": f"'{args.matlab_command}' was not found on PATH"})
            continue
        skipped.append({"language": language, "reason": "Unsupported language tag"})
    if not available:
        raise RuntimeError("No runnable languages were available for the requested sweep.")
    return available, skipped


def write_sweep_case(base_case: dict[str, Any], point: SweepPoint, cases_dir: Path) -> Path:
    sorted_idx = np.argsort(-np.hypot(base_case["arrays"]["a"], base_case["arrays"]["b"]))
    keep = sorted_idx[: point.component_count]
    case_id = f"{base_case['case_id']}_c{point.component_count}_N{point.grid_size}"
    case_dir = cases_dir / case_id
    inputs_dir = case_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)

    arrays = {}
    for name in ("a", "b", "kx", "ky"):
        values = np.asarray(base_case["arrays"][name], dtype=np.float64).reshape(-1)[keep]
        np.savetxt(inputs_dir / f"{name}.csv", values.reshape(-1, 1), delimiter=",")
        arrays[name] = f"inputs/{name}.csv"

    manifest = {
        "case_id": case_id,
        "description": (
            f"Swept benchmark derived from {base_case['case_id']} with {point.component_count} retained components "
            f"and Nx = Ny = {point.grid_size}."
        ),
        "purpose": "benchmark_sweep",
        "inputs": {
            **base_case["inputs"],
            "Nx": point.grid_size,
            "Ny": point.grid_size,
        },
        "arrays": arrays,
    }
    (case_dir / "case.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return case_dir


def run_language(language: str, case_dir: Path, args: argparse.Namespace) -> dict[str, Any]:
    if language == "python":
        case = load_case(case_dir)
        result = run_case(case, repeats=args.repeats, warmup=args.warmup)
        return {
            "case_id": case["case_id"],
            "component_count": int(case["arrays"]["a"].size),
            "domain": {"Lx": float(case["inputs"]["Lx"]), "Ly": float(case["inputs"]["Ly"])},
            "grid": {"Nx": int(case["inputs"]["Nx"]), "Ny": int(case["inputs"]["Ny"])},
            "runtime": result["runtime"],
        }
    if language == "cpp":
        cmd = [str((REPO_ROOT / args.cpp_exe).resolve()), "benchmark", str(case_dir.resolve()), str(args.repeats), "1" if args.warmup else "0"]
        return json.loads(run_command(cmd))
    if language == "matlab":
        matlab_json = case_dir / "matlab_benchmark.json"
        addpath = matlab_quote(str((REPO_ROOT / "matlab" / "repro").resolve()).replace("\\", "/"))
        case_arg = matlab_quote(str(case_dir.resolve()).replace("\\", "/"))
        out_arg = matlab_quote(str(matlab_json.resolve()).replace("\\", "/"))
        batch = f"addpath({addpath}); benchmark_case_from_manifest({case_arg}, {out_arg}, {args.repeats});"
        run_command([args.matlab_command, "-batch", batch], cwd=REPO_ROOT)
        return json.loads(matlab_json.read_text(encoding="utf-8"))
    raise ValueError(f"Unsupported language: {language}")


def matlab_quote(text: str) -> str:
    return "'" + text.replace("'", "''") + "'"


def run_command(cmd: list[str], cwd: Path | None = None) -> str:
    completed = subprocess.run(
        cmd,
        cwd=cwd or REPO_ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    return completed.stdout.strip()


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("No benchmark rows were collected.")
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_figures(rows: list[dict[str, Any]], component_counts: list[int], grid_sizes: list[int], languages: list[str], output_dir: Path) -> dict[str, Path]:
    figures: dict[str, Path] = {}
    path_by_grid = output_dir / "runtime_vs_components_by_grid.png"
    plot_runtime_vs_components(rows, component_counts, grid_sizes, languages, path_by_grid)
    figures["runtime_vs_components_by_grid"] = path_by_grid

    path_by_components = output_dir / "runtime_vs_grid_by_components.png"
    plot_runtime_vs_grid(rows, component_counts, grid_sizes, languages, path_by_components)
    figures["runtime_vs_grid_by_components"] = path_by_components
    return figures


def plot_runtime_vs_components(rows: list[dict[str, Any]], component_counts: list[int], grid_sizes: list[int], languages: list[str], output_path: Path) -> None:
    apply_plot_style()
    ncols = min(2, len(grid_sizes))
    nrows = int(math.ceil(len(grid_sizes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.6 * ncols, 5.2 * nrows), constrained_layout=True, squeeze=False)

    for ax, grid_size in zip(axes.flat, grid_sizes):
        style_axes(ax)
        matlab_y: list[float] | None = None
        cpp_y: list[float] | None = None
        for language in languages:
            y = [lookup_metric(rows, language, component_count, grid_size, "mean_total_s") for component_count in component_counts]
            if language == "matlab":
                matlab_y = y
            if language == "cpp":
                cpp_y = y
            ax.plot(
                component_counts,
                y,
                marker=MARKERS.get(language, "o"),
                linestyle=LINESTYLES.get(language, "-"),
                linewidth=2.6,
                markersize=6.5,
                label=language.title(),
                color=PALETTE.get(language),
                markerfacecolor="white",
                markeredgewidth=1.4,
            )
            annotate_series(ax, component_counts, y, PALETTE.get(language), language, x_axis="components")
        add_cpp_matlab_band(ax, component_counts, matlab_y, cpp_y)
        ax.set_title(f"Nx = Ny = {grid_size}")
        ax.set_xlabel("Retained components")
        ax.set_ylabel("Mean total runtime (s)")
        ax.set_yscale("log")
        add_ratio_note(ax, matlab_y, cpp_y)

    for ax in axes.flat[len(grid_sizes) :]:
        ax.axis("off")

    add_shared_legend(fig, languages)
    fig.suptitle("Runtime vs Retained Components at Fixed Grid Sizes", fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_runtime_vs_grid(rows: list[dict[str, Any]], component_counts: list[int], grid_sizes: list[int], languages: list[str], output_path: Path) -> None:
    apply_plot_style()
    ncols = min(2, len(component_counts))
    nrows = int(math.ceil(len(component_counts) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.6 * ncols, 5.2 * nrows), constrained_layout=True, squeeze=False)

    for ax, component_count in zip(axes.flat, component_counts):
        style_axes(ax)
        matlab_y: list[float] | None = None
        cpp_y: list[float] | None = None
        for language in languages:
            y = [lookup_metric(rows, language, component_count, grid_size, "mean_total_s") for grid_size in grid_sizes]
            if language == "matlab":
                matlab_y = y
            if language == "cpp":
                cpp_y = y
            ax.plot(
                grid_sizes,
                y,
                marker=MARKERS.get(language, "o"),
                linestyle=LINESTYLES.get(language, "-"),
                linewidth=2.6,
                markersize=6.5,
                label=language.title(),
                color=PALETTE.get(language),
                markerfacecolor="white",
                markeredgewidth=1.4,
            )
            annotate_series(ax, grid_sizes, y, PALETTE.get(language), language, x_axis="grid")
        add_cpp_matlab_band(ax, grid_sizes, matlab_y, cpp_y)
        ax.set_title(f"{component_count} retained components")
        ax.set_xlabel("Grid size Nx = Ny")
        ax.set_ylabel("Mean total runtime (s)")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        add_ratio_note(ax, matlab_y, cpp_y)

    for ax in axes.flat[len(component_counts) :]:
        ax.axis("off")

    add_shared_legend(fig, languages)
    fig.suptitle("Runtime vs Grid Size at Fixed Retained-Component Counts", fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def lookup_metric(rows: list[dict[str, Any]], language: str, component_count: int, grid_size: int, metric: str) -> float:
    for row in rows:
        if row["language"] != language:
            continue
        if int(row["component_count"]) != int(component_count):
            continue
        if int(row["Nx"]) != int(grid_size):
            continue
        return float(row[metric])
    return float("nan")


def annotate_series(ax: Any, x_values: list[int], y_values: list[float], color: str | None, language: str, x_axis: str) -> None:
    offsets = {
        ("matlab", "components"): (0, 7),
        ("python", "components"): (0, -13),
        ("cpp", "components"): (0, 7),
        ("matlab", "grid"): (0, 7),
        ("python", "grid"): (0, -13),
        ("cpp", "grid"): (0, 7),
    }
    x_off, y_off = offsets.get((language, x_axis), (0, 6))
    for x, y in zip(x_values, y_values):
        if not np.isfinite(y):
            continue
        ax.annotate(
            f"{y:.2e}",
            xy=(x, y),
            xytext=(x_off, y_off),
            textcoords="offset points",
            ha="center",
            va="bottom" if y_off >= 0 else "top",
            fontsize=7.5,
            color=color or "black",
            bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": "none", "alpha": 0.78},
        )


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10.5,
            "axes.facecolor": "#fcfcf8",
            "figure.facecolor": "white",
            "axes.edgecolor": "#d8d2c2",
            "grid.color": "#d7d2c8",
            "grid.linestyle": ":",
            "grid.linewidth": 0.8,
        }
    )


def style_axes(ax: Any) -> None:
    ax.grid(True, which="major", alpha=0.5)
    ax.grid(True, which="minor", alpha=0.18)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#c9c3b6")
    ax.spines["bottom"].set_color("#c9c3b6")


def add_shared_legend(fig: Any, languages: list[str]) -> None:
    handles = []
    labels = []
    for language in languages:
        handle = plt.Line2D(
            [0],
            [0],
            color=PALETTE.get(language),
            linestyle=LINESTYLES.get(language, "-"),
            marker=MARKERS.get(language, "o"),
            linewidth=2.6,
            markersize=6.5,
            markerfacecolor="white",
            markeredgewidth=1.4,
        )
        handles.append(handle)
        labels.append(language.title())
    legend = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=len(languages), frameon=True)
    legend.get_frame().set_facecolor("#f7f4ec")
    legend.get_frame().set_edgecolor("#d7d2c8")
    legend.get_frame().set_alpha(0.95)


def add_cpp_matlab_band(ax: Any, x_values: list[int], matlab_y: list[float] | None, cpp_y: list[float] | None) -> None:
    if matlab_y is None or cpp_y is None:
        return
    y1 = np.asarray(matlab_y, dtype=float)
    y2 = np.asarray(cpp_y, dtype=float)
    mask = np.isfinite(y1) & np.isfinite(y2)
    if not np.any(mask):
        return
    x = np.asarray(x_values, dtype=float)[mask]
    lower = np.minimum(y1[mask], y2[mask])
    upper = np.maximum(y1[mask], y2[mask])
    ax.fill_between(x, lower, upper, color="#7db8ad", alpha=0.12, zorder=0)


def add_ratio_note(ax: Any, matlab_y: list[float] | None, cpp_y: list[float] | None) -> None:
    if matlab_y is None or cpp_y is None:
        return
    y1 = np.asarray(matlab_y, dtype=float)
    y2 = np.asarray(cpp_y, dtype=float)
    mask = np.isfinite(y1) & np.isfinite(y2) & (y2 > 0.0)
    if not np.any(mask):
        return
    ratio = y1[mask] / y2[mask]
    text = f"MATLAB/C++ ratio: {np.min(ratio):.2f}x to {np.max(ratio):.2f}x"
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        color="#3d4b55",
        bbox={"boxstyle": "round,pad=0.24", "fc": "#f7f4ec", "ec": "#d7d2c8", "alpha": 0.92},
    )


def timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


if __name__ == "__main__":
    main()
