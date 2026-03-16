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
    parser = argparse.ArgumentParser(description="Plot centerline and diagonal component lines for MATLAB vs Python")
    parser.add_argument("--case-dir", default="cross_language_comparison/cases/wavegroup_regression")
    parser.add_argument("--output", default="outputs/cross_language_comparison/figures/wavegroup_regression_component_lines_python_vs_matlab.png")
    args = parser.parse_args()

    case = load_case(args.case_dir)
    matlab = load_matlab_components(Path(args.case_dir))
    python = compute_python_components(case, matlab["k2_cut"])
    make_figure(case, matlab, python, Path(args.output))


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


def compute_python_components(case: dict, k2_cut: float) -> dict:
    inp = case["inputs"]
    arr = case["arrays"]
    coeff1 = spectral_coefficients(1, inp["g"], inp["h"], arr["a"], arr["b"], arr["kx"], arr["ky"], inp["Ux"], inp["Uy"], {"enable_subharmonic": False})
    coeff2 = spectral_coefficients(2, inp["g"], inp["h"], arr["a"], arr["b"], arr["kx"], arr["ky"], inp["Ux"], inp["Uy"], {"enable_subharmonic": False})

    phi1, _, x, y = phi_only(coeff1, inp)
    phi2_total, _, _, _ = phi_only(coeff2, inp)
    phi2_inc = phi2_total - phi1
    phi2sup = split_by_wavenumber(phi2_inc, inp["Lx"], inp["Ly"], k2_cut, "high")
    phi2sub = split_by_wavenumber(phi2_inc, inp["Lx"], inp["Ly"], k2_cut, "low")

    try:
        coeff3 = spectral_coefficients(3, inp["g"], inp["h"], arr["a"], arr["b"], arr["kx"], arr["ky"], inp["Ux"], inp["Uy"], {"enable_subharmonic": False})
        phi3_total, _, _, _ = phi_only(coeff3, inp)
        phi3 = phi3_total - phi2_total
    except NotImplementedError:
        phi3 = np.zeros_like(phi1)

    x_axis = case["reference"]["loaded_arrays"]["x"].reshape(-1)
    y_axis = case["reference"]["loaded_arrays"]["y"].reshape(-1)
    x_grid, y_grid = np.meshgrid(x_axis, y_axis)
    iyc = int(np.argmin(np.abs(y_axis - inp["Ly"] / 2.0)))
    s = np.linspace(0.0, min(inp["Lx"], inp["Ly"]), inp["Nx"])

    center = np.vstack([phi1[iyc, :], phi2sup[iyc, :], phi2sub[iyc, :], phi3[iyc, :]])
    diag = np.vstack(
        [
            interp_line(x_grid, y_grid, phi1, s),
            interp_line(x_grid, y_grid, phi2sup, s),
            interp_line(x_grid, y_grid, phi2sub, s),
            interp_line(x_grid, y_grid, phi3, s),
        ]
    )
    return {"center": center, "diag": diag, "x": x_axis, "s": s}


def phi_only(coeffs: dict, inp: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    eta, phi, xg, yg = spectral_surface(coeffs, inp["Lx"], inp["Ly"], inp["Nx"], inp["Ny"], inp["t"])
    return phi, eta, xg, yg


def split_by_wavenumber(phi: np.ndarray, lx: float, ly: float, k_cut: float, mode: str) -> np.ndarray:
    ny, nx = phi.shape
    dkx = 2.0 * np.pi / lx
    dky = 2.0 * np.pi / ly
    kx_idx = np.concatenate((np.arange(0, nx // 2), np.arange(-nx // 2, 0)))
    ky_idx = np.concatenate((np.arange(0, ny // 2), np.arange(-ny // 2, 0)))
    kx, ky = np.meshgrid(kx_idx * dkx, ky_idx * dky)
    kval = np.hypot(kx, ky)
    spec = np.fft.fft2(phi)
    if mode == "high":
        mask = kval >= k_cut
    else:
        mask = kval < k_cut
    return np.fft.ifft2(spec * mask).real


def interp_line(x_grid: np.ndarray, y_grid: np.ndarray, field: np.ndarray, s: np.ndarray) -> np.ndarray:
    x_axis = x_grid[0, :]
    y_axis = y_grid[:, 0]
    out = np.zeros_like(s)
    for i, val in enumerate(s):
        out[i] = bilinear_interp(x_axis, y_axis, field, val, val)
    return out


def bilinear_interp(x_axis: np.ndarray, y_axis: np.ndarray, field: np.ndarray, x: float, y: float) -> float:
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


def make_figure(case: dict, matlab: dict, python: dict, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.0, 0.55])
    axes = np.empty((2, 4), dtype=object)
    for row in range(2):
        for col in range(4):
            axes[row, col] = fig.add_subplot(gs[row, col])
    info_ax = fig.add_subplot(gs[2, :])
    blue = "#0a4a8a"
    orange = "#c24a0a"
    eps = np.finfo(float).eps

    center_err = np.abs(python["center"] - matlab["center"])
    diag_err = np.abs(python["diag"] - matlab["diag"])
    center_max = np.array([finite_max(center_err[j]) for j in range(4)])
    diag_max = np.array([finite_max(diag_err[j]) for j in range(4)])

    for j, label in enumerate(LABELS):
        ax = axes[0, j]
        ax.plot(matlab["x"], matlab["center"][j], "-", color=blue, linewidth=2.2, label="MATLAB spectral")
        ax.plot(python["x"], python["center"][j], "--", color=orange, linewidth=2.2, label="Python spectral")
        ax.set_title(f"{label}\ncenter max|diff|={center_max[j]:.2e}")
        ax.set_xlabel("Centerline x (m)")
        if j == 0:
            ax.set_ylabel(r"Centerline $\phi_s$")
        ax.grid(True, alpha=0.25)

    for j, label in enumerate(LABELS):
        ax = axes[1, j]
        ax.plot(matlab["s"], matlab["diag"][j], "-", color=blue, linewidth=2.2, label="MATLAB spectral")
        ax.plot(python["s"], python["diag"][j], "--", color=orange, linewidth=2.2, label="Python spectral")
        ax.set_title(f"diag max|diff|={diag_max[j]:.2e}")
        ax.set_xlabel("Diagonal distance s (m)")
        if j == 0:
            ax.set_ylabel(r"Diagonal $\phi_s$")
        ax.grid(True, alpha=0.25)

    info_ax.axis("off")
    summary_text = build_summary_text(case, matlab, python, center_err, diag_err, center_max, diag_max, eps)
    info_ax.text(0.01, 0.98, summary_text, va="top", ha="left", family="monospace", fontsize=10)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    legend = fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=2,
        frameon=True,
        fontsize=11,
    )
    legend.get_frame().set_edgecolor("#cccccc")
    legend.get_frame().set_linewidth(0.8)
    legend.get_frame().set_alpha(0.95)
    fig.suptitle(
        "Wavegroup Component Summary: MATLAB spectral vs Python spectral\n"
        "Top row: centerline, middle row: diagonal; bottom row: setup and machine-precision summary",
        fontsize=15,
    )
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    print(json.dumps({"figure": str(output_path.resolve())}, indent=2))


def build_summary_text(case: dict, matlab: dict, python: dict, center_err: np.ndarray, diag_err: np.ndarray, center_max: np.ndarray, diag_max: np.ndarray, eps: float) -> str:
    inp = case["inputs"]
    n_comp = int(case["arrays"]["a"].size)
    labels = LABELS
    wavegroup_lines = [
        f"case_id: {case['case_id']}",
        f"description: {case.get('description', 'n/a')}",
        f"order={inp['order']}  grid={inp['Nx']}x{inp['Ny']}  N_components={n_comp}",
        f"domain=({inp['Lx']:.0f} m, {inp['Ly']:.0f} m)  h={inp['h']:.6g} m  t={inp['t']:.6g} s",
        f"subharmonic_mode={inp.get('subharmonic_mode', 'skip')}",
    ]

    if case["case_id"] == "wavegroup_regression":
        wavegroup_lines.extend(
            [
                "directional setup: Tp=12 s, theta1=+25 deg, theta2=-35 deg",
                "spreading: sigma_theta1=20 deg, sigma_theta2=24 deg, weights=(0.55, 0.45)",
                "amplitude setup: A0=0.25, focused at domain center",
            ]
        )

    max_center = [float(center_max[j]) for j in range(4)]
    max_diag = [float(diag_max[j]) for j in range(4)]
    phi_all = np.concatenate([center_err[np.isfinite(center_err)], diag_err[np.isfinite(diag_err)]])
    max_all = float(np.max(phi_all))
    rms_all = float(np.sqrt(np.mean(phi_all * phi_all)))
    diff_lines = [
        "",
        "component max|Python - MATLAB|:",
    ]
    for idx, label in enumerate(labels):
        diff_lines.append(
            f"  {label:<21} center={max_center[idx]:.2e}   diag={max_diag[idx]:.2e}"
        )

    diff_lines.extend(
        [
            "",
            f"global line max|diff| = {max_all:.2e}",
            f"global line rms diff  = {rms_all:.2e}",
            f"machine epsilon       = {eps:.2e}",
        ]
    )

    return "\n".join(wavegroup_lines + diff_lines)


def finite_max(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    return float(np.max(finite))


if __name__ == "__main__":
    main()
