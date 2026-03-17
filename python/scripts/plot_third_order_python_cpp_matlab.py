from __future__ import annotations

import argparse
import json
from pathlib import Path
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import numpy as np

from mf12_python.io import load_case
from mf12_python.spectral import spectral_coefficients, spectral_surface


def generated_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot third-superharmonic lines for MATLAB vs Python vs C++ with error curves")
    parser.add_argument("--case-dir", default="cross_language_comparison/cases/wavegroup_regression")
    parser.add_argument("--cpp-root", default="outputs/cross_language_comparison/verify_cpp_orders")
    parser.add_argument("--output", default="outputs/cross_language_comparison/figures/wavegroup_regression_third_order_python_cpp_matlab.png")
    args = parser.parse_args()

    case = load_case(args.case_dir)
    matlab = load_matlab_components(Path(args.case_dir))
    python = compute_python_third(case, matlab["s"])
    cpp = load_cpp_third(case, matlab["s"], Path(args.cpp_root))
    make_figure(case, matlab, python, cpp, Path(args.output))


def load_matlab_components(case_dir: Path) -> dict:
    comp_dir = case_dir / "reference" / "matlab_components"
    meta = json.loads((comp_dir / "components.json").read_text(encoding="utf-8"))
    center = np.loadtxt(comp_dir / "centerline_phi_components.csv", delimiter=",")
    diag = np.loadtxt(comp_dir / "diagonal_phi_components.csv", delimiter=",")
    return {
        "center": center[3],
        "diag": diag[3],
        "x": np.loadtxt(comp_dir / "x.csv", delimiter=",").reshape(-1),
        "s": np.loadtxt(comp_dir / "s.csv", delimiter=",").reshape(-1),
        "k2_cut": float(meta["k2_cut"]),
    }


def compute_python_third(case: dict, s_ref: np.ndarray) -> dict:
    inp = case["inputs"]
    arr = case["arrays"]
    coeff2 = spectral_coefficients(2, inp["g"], inp["h"], arr["a"], arr["b"], arr["kx"], arr["ky"], inp["Ux"], inp["Uy"], {"enable_subharmonic": False})
    coeff3 = spectral_coefficients(3, inp["g"], inp["h"], arr["a"], arr["b"], arr["kx"], arr["ky"], inp["Ux"], inp["Uy"], {"enable_subharmonic": False})
    _, phi2, xg, yg = spectral_surface(coeff2, inp["Lx"], inp["Ly"], inp["Nx"], inp["Ny"], inp["t"])
    _, phi3, _, _ = spectral_surface(coeff3, inp["Lx"], inp["Ly"], inp["Nx"], inp["Ny"], inp["t"])
    phi_third = phi3 - phi2
    x_axis = xg[0, :]
    y_axis = yg[:, 0]
    return extract_lines(inp, x_axis, y_axis, phi_third, s_ref)


def load_cpp_third(case: dict, s_ref: np.ndarray, cpp_root: Path) -> dict:
    inp = case["inputs"]
    phi2 = np.loadtxt(cpp_root / "order2" / "phi.csv", delimiter=",")
    phi3 = np.loadtxt(cpp_root / "order3" / "phi.csv", delimiter=",")
    x_axis = np.loadtxt(cpp_root / "order3" / "x.csv", delimiter=",").reshape(-1)
    y_axis = np.loadtxt(cpp_root / "order3" / "y.csv", delimiter=",").reshape(-1)
    phi_third = phi3 - phi2
    return extract_lines(inp, x_axis, y_axis, phi_third, s_ref)


def extract_lines(inp: dict, x_axis: np.ndarray, y_axis: np.ndarray, phi_third: np.ndarray, s_ref: np.ndarray) -> dict:
    iyc = int(np.argmin(np.abs(y_axis - inp["Ly"] / 2.0)))
    s = np.asarray(s_ref, dtype=float).reshape(-1)
    center = phi_third[iyc, :]
    diag = interp_line(x_axis, y_axis, phi_third, s)
    return {"center": center, "diag": diag, "x": x_axis, "s": s}


def interp_line(x_axis: np.ndarray, y_axis: np.ndarray, field: np.ndarray, s: np.ndarray) -> np.ndarray:
    out = np.zeros_like(s)
    for i, val in enumerate(s):
        out[i] = bilinear_interp(x_axis, y_axis, field, val, val)
    return out


def bilinear_interp(x_axis: np.ndarray, y_axis: np.ndarray, field: np.ndarray, x: float, y: float) -> float:
    if x < x_axis[0] or x > x_axis[-1] or y < y_axis[0] or y > y_axis[-1]:
        return float("nan")
    ix = np.searchsorted(x_axis, x, side="right") - 1
    iy = np.searchsorted(y_axis, y, side="right") - 1
    ix = np.clip(ix, 0, len(x_axis) - 2)
    iy = np.clip(iy, 0, len(y_axis) - 2)
    x1, x2 = x_axis[ix], x_axis[ix + 1]
    y1, y2 = y_axis[iy], y_axis[iy + 1]
    q11 = field[iy, ix]
    q21 = field[iy, ix + 1]
    q12 = field[iy + 1, ix]
    q22 = field[iy + 1, ix + 1]
    tx = 0.0 if x2 == x1 else (x - x1) / (x2 - x1)
    ty = 0.0 if y2 == y1 else (y - y1) / (y2 - y1)
    return (1 - tx) * (1 - ty) * q11 + tx * (1 - ty) * q21 + (1 - tx) * ty * q12 + tx * ty * q22


def finite_max(values: np.ndarray) -> float:
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    return float(np.max(np.abs(finite)))


def fmt_max(values: np.ndarray) -> str:
    value = finite_max(values)
    if not np.isfinite(value):
        return "n/a"
    return f"{value:.2e}"


def make_figure(case: dict, matlab: dict, python: dict, cpp: dict, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    blue = "#0a4a8a"
    orange = "#c24a0a"
    green = "#1b7f3a"
    red = "#b22222"

    py_center_err = python["center"] - matlab["center"]
    cpp_center_err = cpp["center"] - matlab["center"]
    py_diag_err = python["diag"] - matlab["diag"]
    cpp_diag_err = cpp["diag"] - matlab["diag"]

    axes[0, 0].plot(matlab["x"], matlab["center"], "-", color=blue, linewidth=2.2, label="MATLAB spectral")
    axes[0, 0].plot(python["x"], python["center"], "--", color=orange, linewidth=2.0, label="Python spectral")
    axes[0, 0].plot(cpp["x"], cpp["center"], ":", color=green, linewidth=2.4, label="C++ spectral")
    axes[0, 0].set_title("Third Superharmonic: Centerline")
    axes[0, 0].set_xlabel("Centerline x (m)")
    axes[0, 0].set_ylabel(r"$\phi_s^{(3)}$")
    axes[0, 0].grid(True, alpha=0.25)

    diag_mask = np.isfinite(matlab["diag"]) & np.isfinite(python["diag"]) & np.isfinite(cpp["diag"])
    axes[0, 1].plot(matlab["s"][diag_mask], matlab["diag"][diag_mask], "-", color=blue, linewidth=2.2, label="MATLAB spectral")
    axes[0, 1].plot(python["s"][diag_mask], python["diag"][diag_mask], "--", color=orange, linewidth=2.0, label="Python spectral")
    axes[0, 1].plot(cpp["s"][diag_mask], cpp["diag"][diag_mask], ":", color=green, linewidth=2.4, label="C++ spectral")
    axes[0, 1].set_title("Third Superharmonic: Diagonal")
    axes[0, 1].set_xlabel("Diagonal distance s (m)")
    axes[0, 1].set_ylabel(r"$\phi_s^{(3)}$")
    axes[0, 1].grid(True, alpha=0.25)

    axes[1, 0].plot(python["x"], py_center_err, "--", color=orange, linewidth=2.0, label="Python - MATLAB")
    axes[1, 0].plot(cpp["x"], cpp_center_err, "-", color=red, linewidth=2.0, label="C++ - MATLAB")
    axes[1, 0].set_title(
        f"Centerline Error\n"
        f"Py max|diff|={fmt_max(py_center_err)}   "
        f"C++ max|diff|={fmt_max(cpp_center_err)}"
    )
    axes[1, 0].set_xlabel("Centerline x (m)")
    axes[1, 0].set_ylabel("Difference")
    axes[1, 0].grid(True, alpha=0.25)

    axes[1, 1].plot(python["s"], py_diag_err, "--", color=orange, linewidth=2.0, label="Python - MATLAB")
    axes[1, 1].plot(cpp["s"], cpp_diag_err, "-", color=red, linewidth=2.0, label="C++ - MATLAB")
    axes[1, 1].set_title(
        f"Diagonal Error\n"
        f"Py max|diff|={fmt_max(py_diag_err)}   "
        f"C++ max|diff|={fmt_max(cpp_diag_err)}"
    )
    axes[1, 1].set_xlabel("Diagonal distance s (m)")
    axes[1, 1].set_ylabel("Difference")
    axes[1, 1].grid(True, alpha=0.25)

    inp = case["inputs"]
    n_comp = int(case["arrays"]["a"].size)
    summary = (
        f"case_id: {case['case_id']}\n"
        f"order={inp['order']}  grid={inp['Nx']}x{inp['Ny']}  N_components={n_comp}\n"
        f"domain=({inp['Lx']:.0f} m, {inp['Ly']:.0f} m)  h={inp['h']:.6g} m  t={inp['t']:.6g} s\n"
        f"directional setup: Tp=12 s, theta1=+25 deg, theta2=-35 deg\n"
        f"spreading: sigma_theta1=20 deg, sigma_theta2=24 deg, weights=(0.55, 0.45)\n"
        f"global Py max|diff|={max(finite_max(py_center_err), finite_max(py_diag_err)):.2e}  "
        f"global C++ max|diff|={max(finite_max(cpp_center_err), finite_max(cpp_diag_err)):.2e}"
    )
    fig.text(0.5, 0.01, summary, ha="center", va="bottom", family="monospace", fontsize=10)
    fig.text(
        0.995,
        0.995,
        f"Generated: {generated_timestamp()}",
        ha="right",
        va="top",
        fontsize=9,
        color="#555555",
    )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    err_handles, err_labels = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles + err_handles, labels + err_labels, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=5, frameon=True, fontsize=10)
    fig.suptitle("Wavegroup Third-Superharmonic Comparison: MATLAB vs Python vs C++", fontsize=15)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    print(json.dumps({"figure": str(output_path.resolve())}, indent=2))


if __name__ == "__main__":
    main()
