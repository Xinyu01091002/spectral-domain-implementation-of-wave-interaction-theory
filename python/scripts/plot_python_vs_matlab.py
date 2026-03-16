from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mf12_python.compare import compare_fields
from mf12_python.io import load_case
from mf12_python.spectral import run_case


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Python vs MATLAB MF12 case comparisons")
    parser.add_argument("case_dirs", nargs="+", help="Case directories under cross_language_comparison/cases/")
    parser.add_argument("--output-dir", default="outputs/cross_language_comparison/figures")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    for case_dir in args.case_dirs:
        case = load_case(case_dir)
        result = run_case(case, repeats=1, warmup=False)
        summary = make_case_figure(case, result, out_dir)
        summaries.append(summary)

    summary_path = out_dir / "python_vs_matlab_summary.json"
    summary_path.write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    print(json.dumps({"output_dir": str(out_dir.resolve()), "summary": str(summary_path.resolve()), "cases": summaries}, indent=2))


def make_case_figure(case: dict, result: dict, out_dir: Path) -> dict:
    ref = case["reference"]["loaded_arrays"]
    x = ref["x"].reshape(-1)
    y = ref["y"].reshape(-1)
    eta_py = result["eta"]
    phi_py = result["phi"]
    eta_m = ref["eta"]
    phi_m = ref["phi"]

    eta_metrics = compare_fields(eta_py, eta_m, "eta")
    phi_metrics = compare_fields(phi_py, phi_m, "phi")

    fig, axes = plt.subplots(3, 4, figsize=(16, 10), constrained_layout=True)
    cmap_main = "viridis"
    cmap_diff = "coolwarm"
    ix_mid = len(x) // 2
    iy_mid = len(y) // 2

    draw_field(axes[0, 0], x, y, eta_m, "MATLAB eta", cmap_main)
    draw_field(axes[0, 1], x, y, eta_py, "Python eta", cmap_main)
    draw_field(axes[0, 2], x, y, eta_py - eta_m, "eta diff", cmap_diff, symmetric=True)
    draw_line(axes[0, 3], x, eta_m[iy_mid, :], eta_py[iy_mid, :], "eta centerline", eta_metrics["eta_max_abs_err"])

    draw_field(axes[1, 0], x, y, phi_m, "MATLAB phi", cmap_main)
    draw_field(axes[1, 1], x, y, phi_py, "Python phi", cmap_main)
    draw_field(axes[1, 2], x, y, phi_py - phi_m, "phi diff", cmap_diff, symmetric=True)
    draw_line(axes[1, 3], x, phi_m[iy_mid, :], phi_py[iy_mid, :], "phi centerline", phi_metrics["phi_max_abs_err"])

    draw_line(axes[2, 0], y, eta_m[:, ix_mid], eta_py[:, ix_mid], "eta cross-section", eta_metrics["eta_rms_err"], xlabel="y (m)")
    draw_line(axes[2, 1], y, phi_m[:, ix_mid], phi_py[:, ix_mid], "phi cross-section", phi_metrics["phi_rms_err"], xlabel="y (m)")
    draw_hist(axes[2, 2], eta_py - eta_m, phi_py - phi_m)
    draw_text_panel(axes[2, 3], case, result, eta_metrics, phi_metrics)

    fig.suptitle(f"Python vs MATLAB: {case['case_id']}", fontsize=16)
    png_path = out_dir / f"{case['case_id']}_python_vs_matlab.png"
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    return {
        "case_id": case["case_id"],
        "figure": str(png_path.resolve()),
        "eta_max_abs_err": eta_metrics["eta_max_abs_err"],
        "phi_max_abs_err": phi_metrics["phi_max_abs_err"],
    }


def draw_field(ax, x, y, z, title, cmap, symmetric: bool = False) -> None:
    if symmetric:
        vmax = float(np.max(np.abs(z)))
        vmin = -vmax
    else:
        vmin = float(np.min(z))
        vmax = float(np.max(z))
    im = ax.imshow(z, origin="lower", extent=[x[0], x[-1], y[0], y[-1]], aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    plt.colorbar(im, ax=ax, shrink=0.85)


def draw_line(ax, x, matlab_line, python_line, title, metric, xlabel: str = "x (m)") -> None:
    ax.plot(x, matlab_line, "-", color="#0a4a8a", linewidth=2.0, label="MATLAB")
    ax.plot(x, python_line, "--", color="#c24a0a", linewidth=2.0, label="Python")
    ax.set_title(f"{title}\nmax/rms metric = {metric:.3e}")
    ax.set_xlabel(xlabel)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)


def draw_hist(ax, eta_diff, phi_diff) -> None:
    ax.hist(eta_diff.ravel(), bins=40, alpha=0.65, label="eta diff", color="#2a7f62")
    ax.hist(phi_diff.ravel(), bins=40, alpha=0.50, label="phi diff", color="#7b4fa3")
    ax.set_title("Difference distribution")
    ax.set_xlabel("Python - MATLAB")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)


def draw_text_panel(ax, case, result, eta_metrics, phi_metrics) -> None:
    ax.axis("off")
    text = "\n".join(
        [
            f"case_id: {case['case_id']}",
            f"purpose: {case.get('purpose', 'n/a')}",
            f"order: {case['inputs']['order']}",
            f"grid: {case['inputs']['Nx']} x {case['inputs']['Ny']}",
            f"N components: {case['arrays']['a'].size}",
            "",
            f"eta max abs err: {eta_metrics['eta_max_abs_err']:.3e}",
            f"eta rel L2 err : {eta_metrics['eta_relative_l2_err']:.3e}",
            f"phi max abs err: {phi_metrics['phi_max_abs_err']:.3e}",
            f"phi rel L2 err : {phi_metrics['phi_relative_l2_err']:.3e}",
            "",
            f"python coeff mean: {result['runtime']['mean_coefficient_s']:.3e} s",
            f"python recon mean: {result['runtime']['mean_reconstruction_s']:.3e} s",
            f"python total mean: {result['runtime']['mean_total_s']:.3e} s",
        ]
    )
    ax.text(0.0, 1.0, text, va="top", ha="left", family="monospace", fontsize=10)


if __name__ == "__main__":
    main()
