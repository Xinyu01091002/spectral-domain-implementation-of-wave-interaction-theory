from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mf12_python.io import load_case
from mf12_python.spectral import spectral_coefficients, spectral_surface


LABELS = [
    "First Harmonic",
    "Second Superharmonic",
    "Second Subharmonic",
    "Third Superharmonic",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot centerline and diagonal component lines for MATLAB vs Python vs C++")
    parser.add_argument("--case-dir", default="cross_language_comparison/cases/wavegroup_regression")
    parser.add_argument("--cpp-root", default="outputs/cross_language_comparison/verify_cpp_orders")
    parser.add_argument("--output", default="outputs/cross_language_comparison/figures/wavegroup_regression_component_lines_python_cpp_matlab.png")
    args = parser.parse_args()

    case = load_case(args.case_dir)
    matlab = load_matlab_components(Path(args.case_dir))
    python = compute_python_components(case, matlab["k2_cut"], matlab["s"])
    cpp = load_cpp_components(case, matlab["k2_cut"], matlab["s"], Path(args.cpp_root))
    make_figure(case, matlab, python, cpp, Path(args.output))


def load_matlab_components(case_dir: Path) -> dict:
    comp_dir = case_dir / "reference" / "matlab_components"
    meta = json.loads((comp_dir / "components.json").read_text(encoding="utf-8"))
    return {
        "center": np.loadtxt(comp_dir / "centerline_phi_components.csv", delimiter=","),
        "diag": np.loadtxt(comp_dir / "diagonal_phi_components.csv", delimiter=","),
        "x": np.loadtxt(comp_dir / "x.csv", delimiter=",").reshape(-1),
        "s": np.loadtxt(comp_dir / "s.csv", delimiter=",").reshape(-1),
        "labels": meta["labels"],
        "k2_cut": float(meta["k2_cut"]),
    }


def compute_python_components(case: dict, k2_cut: float, s_ref: np.ndarray) -> dict:
    inp = case["inputs"]
    arr = case["arrays"]
    coeff1 = spectral_coefficients(1, inp["g"], inp["h"], arr["a"], arr["b"], arr["kx"], arr["ky"], inp["Ux"], inp["Uy"], {"enable_subharmonic": False})
    coeff2 = spectral_coefficients(2, inp["g"], inp["h"], arr["a"], arr["b"], arr["kx"], arr["ky"], inp["Ux"], inp["Uy"], {"enable_subharmonic": False})
    coeff3 = spectral_coefficients(3, inp["g"], inp["h"], arr["a"], arr["b"], arr["kx"], arr["ky"], inp["Ux"], inp["Uy"], {"enable_subharmonic": False})

    phi1, x_axis, y_axis = phi_only(coeff1, inp)
    phi2_total, _, _ = phi_only(coeff2, inp)
    phi3_total, _, _ = phi_only(coeff3, inp)

    return split_components(inp, x_axis, y_axis, phi1, phi2_total, phi3_total, k2_cut, s_ref)


def load_cpp_components(case: dict, k2_cut: float, s_ref: np.ndarray, cpp_root: Path) -> dict:
    order1 = cpp_root / "order1"
    order2 = cpp_root / "order2"
    order3 = cpp_root / "order3"
    phi1 = np.loadtxt(order1 / "phi.csv", delimiter=",")
    phi2_total = np.loadtxt(order2 / "phi.csv", delimiter=",")
    phi3_total = np.loadtxt(order3 / "phi.csv", delimiter=",")
    x_axis = np.loadtxt(order3 / "x.csv", delimiter=",").reshape(-1)
    y_axis = np.loadtxt(order3 / "y.csv", delimiter=",").reshape(-1)
    inp = case["inputs"]
    return split_components(inp, x_axis, y_axis, phi1, phi2_total, phi3_total, k2_cut, s_ref)


def split_components(inp: dict, x_axis: np.ndarray, y_axis: np.ndarray, phi1: np.ndarray, phi2_total: np.ndarray, phi3_total: np.ndarray, k2_cut: float, s_ref: np.ndarray) -> dict:
    phi2_inc = phi2_total - phi1
    phi2sup = split_by_wavenumber(phi2_inc, inp["Lx"], inp["Ly"], k2_cut, "high")
    phi2sub = split_by_wavenumber(phi2_inc, inp["Lx"], inp["Ly"], k2_cut, "low")
    phi3sup = phi3_total - phi2_total

    iyc = int(np.argmin(np.abs(y_axis - inp["Ly"] / 2.0)))
    s = np.asarray(s_ref, dtype=float).reshape(-1)

    center = np.vstack([phi1[iyc, :], phi2sup[iyc, :], phi2sub[iyc, :], phi3sup[iyc, :]])
    diag = np.vstack(
        [
            interp_line(x_axis, y_axis, phi1, s),
            interp_line(x_axis, y_axis, phi2sup, s),
            interp_line(x_axis, y_axis, phi2sub, s),
            interp_line(x_axis, y_axis, phi3sup, s),
        ]
    )
    return {"center": center, "diag": diag, "x": x_axis, "s": s}


def phi_only(coeffs: dict, inp: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    _, phi, xg, yg = spectral_surface(coeffs, inp["Lx"], inp["Ly"], inp["Nx"], inp["Ny"], inp["t"])
    return phi, xg[0, :], yg[:, 0]


def split_by_wavenumber(phi: np.ndarray, lx: float, ly: float, k_cut: float, mode: str) -> np.ndarray:
    ny, nx = phi.shape
    dkx = 2.0 * np.pi / lx
    dky = 2.0 * np.pi / ly
    kx_idx = np.concatenate((np.arange(0, nx // 2), np.arange(-nx // 2, 0)))
    ky_idx = np.concatenate((np.arange(0, ny // 2), np.arange(-ny // 2, 0)))
    kx, ky = np.meshgrid(kx_idx * dkx, ky_idx * dky)
    kval = np.hypot(kx, ky)
    spec = np.fft.fft2(phi)
    mask = kval >= k_cut if mode == "high" else kval < k_cut
    return np.fft.ifft2(spec * mask).real


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
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    return float(np.max(finite))


def fmt_max(values: np.ndarray) -> str:
    value = finite_max(values)
    if not np.isfinite(value):
        return "n/a"
    return f"{value:.2e}"


def build_summary_text(case: dict, matlab: dict, python: dict, cpp: dict) -> str:
    inp = case["inputs"]
    n_comp = int(case["arrays"]["a"].size)
    py_center_err = np.abs(python["center"] - matlab["center"])
    py_diag_err = np.abs(python["diag"] - matlab["diag"])
    cpp_center_err = np.abs(cpp["center"] - matlab["center"])
    cpp_diag_err = np.abs(cpp["diag"] - matlab["diag"])

    lines = [
        f"case_id: {case['case_id']}",
        f"description: {case.get('description', 'n/a')}",
        f"order={inp['order']}  grid={inp['Nx']}x{inp['Ny']}  N_components={n_comp}",
        f"domain=({inp['Lx']:.0f} m, {inp['Ly']:.0f} m)  h={inp['h']:.6g} m  t={inp['t']:.6g} s",
        f"subharmonic_mode={inp.get('subharmonic_mode', 'skip')}",
    ]
    if case["case_id"] == "wavegroup_regression":
        lines.extend(
            [
                "directional setup: Tp=12 s, theta1=+25 deg, theta2=-35 deg",
                "spreading: sigma_theta1=20 deg, sigma_theta2=24 deg, weights=(0.55, 0.45)",
                "amplitude setup: A0=0.25, focused at domain center",
            ]
        )

    lines.append("")
    lines.append("max|implementation - MATLAB| by component:")
    for j, label in enumerate(LABELS):
        lines.append(
            f"  {label:<21} "
            f"Py center={finite_max(py_center_err[j]):.2e} diag={finite_max(py_diag_err[j]):.2e}   "
            f"C++ center={finite_max(cpp_center_err[j]):.2e} diag={finite_max(cpp_diag_err[j]):.2e}"
        )

    py_all = np.concatenate([py_center_err[np.isfinite(py_center_err)], py_diag_err[np.isfinite(py_diag_err)]])
    cpp_all = np.concatenate([cpp_center_err[np.isfinite(cpp_center_err)], cpp_diag_err[np.isfinite(cpp_diag_err)]])
    lines.extend(
        [
            "",
            f"global Python line max|diff| = {np.max(py_all):.2e}",
            f"global Python line rms diff  = {np.sqrt(np.mean(py_all * py_all)):.2e}",
            f"global C++ line max|diff|    = {np.max(cpp_all):.2e}",
            f"global C++ line rms diff     = {np.sqrt(np.mean(cpp_all * cpp_all)):.2e}",
        ]
    )
    return "\n".join(lines)


def make_figure(case: dict, matlab: dict, python: dict, cpp: dict, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(17, 10), constrained_layout=True)
    gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.0, 0.7])
    axes = np.empty((2, 4), dtype=object)
    for row in range(2):
        for col in range(4):
            axes[row, col] = fig.add_subplot(gs[row, col])
    info_ax = fig.add_subplot(gs[2, :])

    blue = "#0a4a8a"
    orange = "#c24a0a"
    green = "#1b7f3a"

    for j, label in enumerate(LABELS):
        ax = axes[0, j]
        ax.plot(matlab["x"], matlab["center"][j], "-", color=blue, linewidth=2.2, label="MATLAB spectral")
        ax.plot(python["x"], python["center"][j], "--", color=orange, linewidth=2.0, label="Python spectral")
        ax.plot(cpp["x"], cpp["center"][j], ":", color=green, linewidth=2.4, label="C++ spectral")
        ax.set_title(
            f"{label}\n"
            f"Py max|diff|={fmt_max(python['center'][j]-matlab['center'][j])}   "
            f"C++ max|diff|={fmt_max(cpp['center'][j]-matlab['center'][j])}"
        )
        ax.set_xlabel("Centerline x (m)")
        if j == 0:
            ax.set_ylabel(r"Centerline $\phi_s$")
        ax.grid(True, alpha=0.25)

    for j, label in enumerate(LABELS):
        ax = axes[1, j]
        diag_mask = np.isfinite(matlab["diag"][j]) & np.isfinite(python["diag"][j]) & np.isfinite(cpp["diag"][j])
        ax.plot(matlab["s"][diag_mask], matlab["diag"][j][diag_mask], "-", color=blue, linewidth=2.2)
        ax.plot(python["s"][diag_mask], python["diag"][j][diag_mask], "--", color=orange, linewidth=2.0)
        ax.plot(cpp["s"][diag_mask], cpp["diag"][j][diag_mask], ":", color=green, linewidth=2.4)
        ax.set_title(
            f"diag: Py={fmt_max(python['diag'][j]-matlab['diag'][j])}   "
            f"C++={fmt_max(cpp['diag'][j]-matlab['diag'][j])}"
        )
        ax.set_xlabel("Diagonal distance s (m)")
        if j == 0:
            ax.set_ylabel(r"Diagonal $\phi_s$")
        ax.grid(True, alpha=0.25)

    info_ax.axis("off")
    info_ax.text(0.01, 0.98, build_summary_text(case, matlab, python, cpp), va="top", ha="left", family="monospace", fontsize=10)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    legend = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.985), ncol=3, frameon=True, fontsize=11)
    legend.get_frame().set_edgecolor("#cccccc")
    legend.get_frame().set_linewidth(0.8)
    legend.get_frame().set_alpha(0.95)
    fig.suptitle(
        "Wavegroup Component Summary: MATLAB vs Python vs C++\n"
        "Top row: centerline, middle row: diagonal, bottom row: directional setup and error summary",
        fontsize=15,
    )
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    print(json.dumps({"figure": str(output_path.resolve())}, indent=2))


if __name__ == "__main__":
    main()
