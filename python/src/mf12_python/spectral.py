from __future__ import annotations

import math
import platform
import socket
import sys
import time
from typing import Any

import numpy as np


def spectral_coefficients(
    order: int,
    g: float,
    h: float,
    a: np.ndarray,
    b: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    ux: float,
    uy: float,
    opts: dict[str, Any] | None = None,
) -> dict[str, Any]:
    opts = dict(opts or {})
    enable_subharmonic = bool(opts.get("enable_subharmonic", False))
    if enable_subharmonic:
        raise NotImplementedError("Python v0.1 supports the superharmonic-only spectral path. Use subharmonic_mode='skip'.")

    a = np.asarray(a, dtype=np.float64).reshape(-1)
    b = np.asarray(b, dtype=np.float64).reshape(-1)
    kx = np.asarray(kx, dtype=np.float64).reshape(-1)
    ky = np.asarray(ky, dtype=np.float64).reshape(-1)

    n_comp = a.size
    kappa = np.hypot(kx, ky)
    omega1 = np.sqrt(g * kappa * np.tanh(h * kappa))
    omega = kx * ux + ky * uy + omega1
    f = -omega1 / (kappa * np.sinh(h * kappa))
    mu = f * np.cosh(h * kappa)
    kappa_2 = 2.0 * kappa
    c = np.sqrt(a * a + b * b)
    mu_star = np.zeros_like(a)

    coeffs: dict[str, Any] = {
        "g": float(g),
        "h": float(h),
        "N": int(n_comp),
        "a": a,
        "b": b,
        "kx": kx,
        "ky": ky,
        "Ux": float(ux),
        "Uy": float(uy),
        "kappa": kappa,
        "omega1": omega1,
        "omega": omega.copy(),
        "mu": mu,
        "muStar": mu_star,
        "F": f,
        "c": c,
        "kappa_2": kappa_2,
        "superharmonic_only": True,
    }

    if order >= 2:
        g_2 = 0.5 * h * kappa * (2.0 + np.cosh(2.0 * h * kappa)) * coth(h * kappa) / (np.sinh(h * kappa) ** 2)
        f_2 = -0.75 * h * omega1 / (np.sinh(h * kappa) ** 4)
        a_2 = (a * a - b * b) / (2.0 * h)
        b_2 = (a * b) / h
        mu_2 = f_2 * np.cosh(h * kappa_2) - h * omega1
        m_flux = c * c * omega1 / (2.0 * kappa) * coth(h * kappa)

        num_pairs = n_comp * (n_comp - 1) // 2
        len2 = 2 * num_pairs
        omega_npm = np.zeros(len2)
        kx_npm = np.zeros(len2)
        ky_npm = np.zeros(len2)
        kappa_npm = np.zeros(len2)
        gamma_npm = np.zeros(len2)
        f_npm = np.zeros(len2)
        g_npm = np.zeros(len2)
        mu_npm = np.zeros(len2)
        a_npm = np.zeros(len2)
        b_npm = np.zeros(len2)

        pair_count = 0
        for n in range(n_comp):
            for m in range(n + 1, n_comp):
                pair_count += 1
                idx_plus = 2 * pair_count - 2
                idx_minus = idx_plus + 1
                vals_plus = pair_terms(1.0, n, m, omega1, kx, ky, kappa, g, h)
                vals_minus = pair_terms(-1.0, n, m, omega1, kx, ky, kappa, g, h)
                omega_npm[idx_plus], kx_npm[idx_plus], ky_npm[idx_plus], kappa_npm[idx_plus], gamma_npm[idx_plus], f_npm[idx_plus], g_npm[idx_plus], mu_npm[idx_plus] = vals_plus
                omega_npm[idx_minus], kx_npm[idx_minus], ky_npm[idx_minus], kappa_npm[idx_minus], gamma_npm[idx_minus], f_npm[idx_minus], g_npm[idx_minus], mu_npm[idx_minus] = vals_minus
                a_npm[idx_plus] = (a[n] * a[m] - b[n] * b[m]) / h
                b_npm[idx_plus] = (a[m] * b[n] + a[n] * b[m]) / h
                a_npm[idx_minus] = (a[n] * a[m] + b[n] * b[m]) / h
                b_npm[idx_minus] = (a[m] * b[n] - a[n] * b[m]) / h

        coeffs.update(
            {
                "A_2": a_2,
                "B_2": b_2,
                "F_2": f_2,
                "G_2": g_2,
                "mu_2": mu_2,
                "F_npm": f_npm,
                "G_npm": g_npm,
                "A_npm": a_npm,
                "B_npm": b_npm,
                "mu_npm": mu_npm,
                "kappa_npm": kappa_npm,
                "omega_npm": omega_npm,
                "kx_2": 2.0 * kx,
                "ky_2": 2.0 * ky,
                "kx_npm": kx_npm,
                "ky_npm": ky_npm,
                "M": m_flux,
                "gamma_npm": gamma_npm,
            }
        )

    if order == 3:
        coeffs.update(third_order_terms(coeffs, g, h))

    return coeffs


def spectral_surface(coeffs: dict[str, Any], lx: float, ly: float, nx: int, ny: int, t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dx = lx / nx
    dy = ly / ny
    x_axis = np.arange(nx, dtype=np.float64) * dx
    y_axis = np.arange(ny, dtype=np.float64) * dy
    x_grid, y_grid = np.meshgrid(x_axis, y_axis)
    dkx = 2.0 * np.pi / lx
    dky = 2.0 * np.pi / ly
    spec_eta = np.zeros((ny, nx), dtype=np.complex128)
    spec_phi = np.zeros((ny, nx), dtype=np.complex128)

    def accumulate_spectrum(k_in_x: np.ndarray, k_in_y: np.ndarray, values: np.ndarray, target: str) -> None:
        spec = spec_eta if target == "eta" else spec_phi
        ux_idx = np.asarray(k_in_x, dtype=np.float64).reshape(-1) / dkx
        uy_idx = np.asarray(k_in_y, dtype=np.float64).reshape(-1) / dky
        vals = np.asarray(values, dtype=np.complex128).reshape(-1)
        valid = np.isfinite(ux_idx) & np.isfinite(uy_idx) & np.isfinite(vals.real) & np.isfinite(vals.imag)
        ux_loc = ux_idx[valid]
        uy_loc = uy_idx[valid]
        vals_loc = vals[valid]
        if vals_loc.size == 0:
            return
        ix0 = np.floor(ux_loc).astype(np.int64)
        iy0 = np.floor(uy_loc).astype(np.int64)
        fx = ux_loc - ix0
        fy = uy_loc - iy0
        np.add.at(spec, (np.mod(iy0, ny), np.mod(ix0, nx)), vals_loc * (1.0 - fx) * (1.0 - fy))
        np.add.at(spec, (np.mod(iy0, ny), np.mod(ix0 + 1, nx)), vals_loc * fx * (1.0 - fy))
        np.add.at(spec, (np.mod(iy0 + 1, ny), np.mod(ix0, nx)), vals_loc * (1.0 - fx) * fy)
        np.add.at(spec, (np.mod(iy0 + 1, ny), np.mod(ix0 + 1, nx)), vals_loc * fx * fy)

    z_lin = (coeffs["a"] + 1j * coeffs["b"]) * np.exp(-1j * coeffs["omega"] * t)
    accumulate_spectrum(coeffs["kx"], coeffs["ky"], z_lin, "eta")
    accumulate_spectrum(coeffs["kx"], coeffs["ky"], z_lin * (coeffs["mu"] + coeffs["muStar"]) * 1j, "phi")

    if "G_2" in coeffs:
        z_2 = (coeffs["A_2"] + 1j * coeffs["B_2"]) * np.exp(-1j * (2.0 * coeffs["omega"]) * t)
        accumulate_spectrum(2.0 * coeffs["kx"], 2.0 * coeffs["ky"], z_2 * coeffs["G_2"], "eta")
        accumulate_spectrum(2.0 * coeffs["kx"], 2.0 * coeffs["ky"], z_2 * coeffs["mu_2"] * 1j, "phi")
        z_npm = (coeffs["A_npm"] + 1j * coeffs["B_npm"]) * np.exp(-1j * coeffs["omega_npm"] * t)
        accumulate_spectrum(coeffs["kx_npm"], coeffs["ky_npm"], z_npm * coeffs["G_npm"], "eta")
        accumulate_spectrum(coeffs["kx_npm"], coeffs["ky_npm"], z_npm * coeffs["mu_npm"] * 1j, "phi")

    if "G_3" in coeffs:
        z_3 = (coeffs["A_3"] + 1j * coeffs["B_3"]) * np.exp(-1j * (3.0 * coeffs["omega"]) * t)
        accumulate_spectrum(3.0 * coeffs["kx"], 3.0 * coeffs["ky"], z_3 * coeffs["G_3"], "eta")
        accumulate_spectrum(3.0 * coeffs["kx"], 3.0 * coeffs["ky"], z_3 * coeffs["mu_3"] * 1j, "phi")

        z_np2m = (coeffs["A_np2m"] + 1j * coeffs["B_np2m"]) * np.exp(-1j * coeffs["omega_np2m"] * t)
        accumulate_spectrum(coeffs["kx_np2m"], coeffs["ky_np2m"], z_np2m * coeffs["G_np2m"], "eta")
        accumulate_spectrum(coeffs["kx_np2m"], coeffs["ky_np2m"], z_np2m * coeffs["mu_np2m"] * 1j, "phi")

        z_2npm = (coeffs["A_2npm"] + 1j * coeffs["B_2npm"]) * np.exp(-1j * coeffs["omega_2npm"] * t)
        accumulate_spectrum(coeffs["kx_2npm"], coeffs["ky_2npm"], z_2npm * coeffs["G_2npm"], "eta")
        accumulate_spectrum(coeffs["kx_2npm"], coeffs["ky_2npm"], z_2npm * coeffs["mu_2npm"] * 1j, "phi")

        z_npmpp = 2.0 * (coeffs["A_npmpp"] + 1j * coeffs["B_npmpp"]) * np.exp(-1j * coeffs["omega_npmpp"] * t)
        accumulate_spectrum(coeffs["kx_npmpp"], coeffs["ky_npmpp"], z_npmpp * coeffs["G_npmpp"], "eta")
        accumulate_spectrum(coeffs["kx_npmpp"], coeffs["ky_npmpp"], z_npmpp * coeffs["mu_npmpp"] * 1j, "phi")

    eta = np.fft.ifft2(spec_eta).real * (nx * ny)
    phi_wave = np.fft.ifft2(spec_phi).real * (nx * ny)
    phi_s = phi_wave + coeffs["Ux"] * x_grid + coeffs["Uy"] * y_grid
    return eta, phi_s, x_grid, y_grid


def run_case(case: dict[str, Any], repeats: int = 1, warmup: bool = False) -> dict[str, Any]:
    inputs = case["inputs"]
    arrays = case["arrays"]
    coeff_times: list[float] = []
    recon_times: list[float] = []
    if warmup:
        coeffs = spectral_coefficients(inputs["order"], inputs["g"], inputs["h"], arrays["a"], arrays["b"], arrays["kx"], arrays["ky"], inputs["Ux"], inputs["Uy"], {"enable_subharmonic": False})
        spectral_surface(coeffs, inputs["Lx"], inputs["Ly"], inputs["Nx"], inputs["Ny"], inputs["t"])
    eta = phi = x_grid = y_grid = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        coeffs = spectral_coefficients(inputs["order"], inputs["g"], inputs["h"], arrays["a"], arrays["b"], arrays["kx"], arrays["ky"], inputs["Ux"], inputs["Uy"], {"enable_subharmonic": False})
        coeff_times.append(time.perf_counter() - t0)
        t1 = time.perf_counter()
        eta, phi, x_grid, y_grid = spectral_surface(coeffs, inputs["Lx"], inputs["Ly"], inputs["Nx"], inputs["Ny"], inputs["t"])
        recon_times.append(time.perf_counter() - t1)
    return {
        "eta": eta,
        "phi": phi,
        "x": x_grid[0, :],
        "y": y_grid[:, 0],
        "runtime": {
            "repeats": repeats,
            "warmup": warmup,
            "coefficient_s": coeff_times,
            "reconstruction_s": recon_times,
            "total_s": [c + r for c, r in zip(coeff_times, recon_times)],
            "mean_coefficient_s": float(np.mean(coeff_times)),
            "mean_reconstruction_s": float(np.mean(recon_times)),
            "mean_total_s": float(np.mean(np.asarray(coeff_times) + np.asarray(recon_times))),
            "best_total_s": float(np.min(np.asarray(coeff_times) + np.asarray(recon_times))),
        },
        "metadata": runtime_metadata(),
    }


def runtime_metadata() -> dict[str, Any]:
    return {
        "language": "python",
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "fft_backend": "numpy.fft",
        "implementation": "mf12_python spectral-only v0.1",
    }


def coth(x: np.ndarray) -> np.ndarray:
    return 1.0 / np.tanh(x)


def pair_terms(pm: float, n: int, m: int, omega1: np.ndarray, kx: np.ndarray, ky: np.ndarray, kappa: np.ndarray, g: float, h: float) -> tuple[float, ...]:
    omega_out = omega1[n] + pm * omega1[m]
    kx_out = kx[n] + pm * kx[m]
    ky_out = ky[n] + pm * ky[m]
    kappa_out = math.hypot(kx_out, ky_out)
    gamma_out = kappa_out * math.sinh(h * kappa_out)
    beta_out = omega_out * omega_out * math.cosh(h * kappa_out) - g * kappa_out * math.sinh(h * kappa_out)
    f_out = gamma2(omega1[n], kx[n], ky[n], kappa[n], pm * omega1[m], pm * kx[m], pm * ky[m], kappa[m], omega_out, beta_out, g, h)
    g_out = lambda2(omega1[n], kx[n], ky[n], kappa[n], pm * omega1[m], pm * kx[m], pm * ky[m], kappa[m], omega_out, omega_out * math.cosh(h * kappa_out), gamma_out, beta_out, g, h)
    mu_out = f_out * math.cosh(h * kappa_out) - 0.5 * h * (omega1[n] + pm * omega1[m])
    return omega_out, kx_out, ky_out, kappa_out, gamma_out, f_out, g_out, mu_out


def lambda2(omega1n: float, knx: float, kny: float, kappan: float, omega1m: float, kmx: float, kmy: float, kappam: float, omega_npm: float, alpha_npm: float, gamma_npm: float, beta_npm: float, g: float, h: float) -> float:
    knkm = knx * kmx + kny * kmy
    return h / (2.0 * omega1n * omega1m * beta_npm) * (g * alpha_npm * (omega1n * (kappam * kappam + knkm) + omega1m * (kappan * kappan + knkm)) + gamma_npm * (g * g * knkm + omega1n * omega1n * omega1m * omega1m - omega1n * omega1m * omega_npm * omega_npm))


def gamma2(omega1n: float, knx: float, kny: float, kappan: float, omega1m: float, kmx: float, kmy: float, kappam: float, omega_npm: float, beta_npm: float, g: float, h: float) -> float:
    knkm = knx * kmx + kny * kmy
    return h / (2.0 * omega1n * omega1m * beta_npm) * (omega1n * omega1m * omega_npm * (omega_npm * omega_npm - omega1n * omega1m) - g * g * omega1n * (kappam * kappam + 2.0 * knkm) - g * g * omega1m * (kappan * kappan + 2.0 * knkm))


def third_order_terms(coeffs: dict[str, Any], g: float, h: float) -> dict[str, Any]:
    a = coeffs["a"]
    b = coeffs["b"]
    kx = coeffs["kx"]
    ky = coeffs["ky"]
    kappa = coeffs["kappa"]
    omega1 = coeffs["omega1"]
    omega = coeffs["omega"]
    c = coeffs["c"]
    g_2 = coeffs["G_2"]
    f_2 = coeffs["F_2"]
    kappa_2 = coeffs["kappa_2"]
    gamma_2 = kappa_2 * np.sinh(h * kappa_2)
    n_comp = coeffs["N"]

    upsilon = omega1 * kappa * (-13.0 + 24.0 * np.cosh(2.0 * h * kappa) + np.cosh(4.0 * h * kappa)) / (64.0 * np.sinh(h * kappa) ** 5)
    xi = (omega1 * g_2 + f_2 * kappa_2 * np.sinh(h * kappa_2) - g * h * kappa * kappa / (2.0 * omega1)) / (4.0 * h)
    f13 = c * c * upsilon
    mu_star = c * c * xi
    omega_cap = (8.0 + np.cosh(4.0 * h * kappa)) / (16.0 * np.sinh(h * kappa) ** 4)
    omega3 = c * c * kappa * kappa * omega_cap

    for n in range(n_comp):
        for m in list(range(n)) + list(range(n + 1, n_comp)):
            omega_np, _, _, kappa_np, gamma_np, f_np, g_np, _ = pair_terms(1.0, n, m, omega1, kx, ky, kappa, g, h)
            omega_nm, _, _, kappa_nm, gamma_nm, f_nm, g_nm, _ = pair_terms(-1.0, n, m, omega1, kx, ky, kappa, g, h)
            f13[n] += c[m] * c[m] * upsilon_nm(omega1[n], kx[n], ky[n], kappa[n], omega1[m], kx[m], ky[m], kappa[m], f_np, f_nm, g_np, g_nm, kappa_np, kappa_nm, g, h)
            mu_star[n] += c[m] * c[m] * xi_nm(omega1[n], kappa[n], omega1[m], g_np, f_np, gamma_np, g_nm, f_nm, gamma_nm, h, g)
            omega3[n] += c[m] * c[m] * kappa[m] * kappa[m] * omega_nm_fn(omega1[n], kx[n], ky[n], omega1[m], kx[m], ky[m], kappa[m], f_np, f_nm, g_np, g_nm, kappa_np, kappa_nm, g, h)

    omega = omega + omega3 * omega1
    mu_star = mu_star + f13 * np.cosh(h * kappa)

    a_3 = np.zeros(n_comp)
    b_3 = np.zeros(n_comp)
    f_3 = np.zeros(n_comp)
    g_3 = np.zeros(n_comp)
    mu_3 = np.zeros(n_comp)
    for n in range(n_comp):
        a_3[n] = 0.5 * theta_a(a[n], b[n], a[n], b[n], a[n], b[n], h)
        b_3[n] = 0.5 * theta_b(a[n], b[n], a[n], b[n], a[n], b[n], h)
        f_3[n] = (h * h * kappa[n] * omega1[n] / (32.0 * np.sinh(h * kappa[n]) ** 7)) * (-11.0 + 2.0 * np.cosh(2.0 * h * kappa[n]))
        g_3[n] = (3.0 * h * h * kappa[n] * kappa[n] / (128.0 * np.sinh(h * kappa[n]) ** 6)) * (
            14.0 + 15.0 * np.cosh(2.0 * h * kappa[n]) + 6.0 * np.cosh(4.0 * h * kappa[n]) + np.cosh(6.0 * h * kappa[n])
        )
        mu_3[n] = f_3[n] * np.cosh(h * 3.0 * kappa[n]) - g * h * h * kappa[n] * kappa[n] / (4.0 * omega1[n]) + 0.5 * h * (f_2[n] * gamma_2[n] - omega1[n] * g_2[n])

    m_nm_row_odd = np.zeros((n_comp, n_comp), dtype=np.int64)
    m_nm_col_odd = np.zeros((n_comp, n_comp), dtype=np.int64)
    pair_count = 0
    for n in range(n_comp):
        for m in range(n + 1, n_comp):
            pair_count += 1
            m_nm_row_odd[n, m] = 2 * pair_count - 1
    pair_count = 0
    for col in range(1, n_comp):
        for row in range(col):
            pair_count += 1
            m_nm_col_odd[row, col] = 2 * pair_count - 1

    num_pairs = n_comp * (n_comp - 1) // 2
    omega_np2m = np.zeros(num_pairs)
    kx_np2m = np.zeros(num_pairs)
    ky_np2m = np.zeros(num_pairs)
    kappa_np2m = np.zeros(num_pairs)
    a_np2m = np.zeros(num_pairs)
    b_np2m = np.zeros(num_pairs)
    f_np2m = np.zeros(num_pairs)
    g_np2m = np.zeros(num_pairs)
    mu_np2m = np.zeros(num_pairs)

    omega_2npm = np.zeros(num_pairs)
    kx_2npm = np.zeros(num_pairs)
    ky_2npm = np.zeros(num_pairs)
    kappa_2npm = np.zeros(num_pairs)
    a_2npm = np.zeros(num_pairs)
    b_2npm = np.zeros(num_pairs)
    f_2npm = np.zeros(num_pairs)
    g_2npm = np.zeros(num_pairs)
    mu_2npm = np.zeros(num_pairs)

    pair_idx = 0
    for n in range(n_comp):
        for m in range(n + 1, n_comp):
            idx_sum_nm = m_nm_row_odd[n, m] - 1

            omega_np2m[pair_idx] = omega1[n] + 2.0 * omega1[m]
            kx_np2m[pair_idx] = kx[n] + 2.0 * kx[m]
            ky_np2m[pair_idx] = ky[n] + 2.0 * ky[m]
            kappa_np2m[pair_idx] = math.hypot(kx_np2m[pair_idx], ky_np2m[pair_idx])
            alpha_np2m = omega_np2m[pair_idx] * math.cosh(h * kappa_np2m[pair_idx])
            gamma_np2m = kappa_np2m[pair_idx] * math.sinh(h * kappa_np2m[pair_idx])
            beta_np2m = omega_np2m[pair_idx] ** 2 * math.cosh(h * kappa_np2m[pair_idx]) - g * kappa_np2m[pair_idx] * math.sinh(h * kappa_np2m[pair_idx])
            a_np2m[pair_idx] = 0.5 * theta_a(a[n], b[n], a[m], b[m], a[m], b[m], h)
            b_np2m[pair_idx] = 0.5 * theta_b(a[n], b[n], a[m], b[m], a[m], b[m], h)
            g_np2m[pair_idx] = lambda3(
                omega1[n], kx[n], ky[n], kappa[n],
                omega1[m], kx[m], ky[m], kappa[m],
                omega1[m], kx[m], ky[m], kappa[m],
                coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                kappa_2[m], gamma_2[m], g_2[m], f_2[m],
                omega_np2m[pair_idx], alpha_np2m, gamma_np2m, beta_np2m, g, h,
            )
            f_np2m[pair_idx] = gamma3(
                omega1[n], kx[n], ky[n], kappa[n],
                omega1[m], kx[m], ky[m], kappa[m],
                omega1[m], kx[m], ky[m], kappa[m],
                coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                kappa_2[m], gamma_2[m], g_2[m], f_2[m],
                omega_np2m[pair_idx], beta_np2m, g, h,
            )
            mu_np2m[pair_idx] = pi_fn(
                omega1[n], kappa[n], omega1[m], kappa[m], omega1[m], kappa[m],
                coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                gamma_2[m], g_2[m], f_2[m],
                f_np2m[pair_idx], kappa_np2m[pair_idx], g, h,
            )

            omega_2npm[pair_idx] = 2.0 * omega1[n] + omega1[m]
            kx_2npm[pair_idx] = 2.0 * kx[n] + kx[m]
            ky_2npm[pair_idx] = 2.0 * ky[n] + ky[m]
            kappa_2npm[pair_idx] = math.hypot(kx_2npm[pair_idx], ky_2npm[pair_idx])
            alpha_2npm = omega_2npm[pair_idx] * math.cosh(h * kappa_2npm[pair_idx])
            gamma_2npm = kappa_2npm[pair_idx] * math.sinh(h * kappa_2npm[pair_idx])
            beta_2npm = omega_2npm[pair_idx] ** 2 * math.cosh(h * kappa_2npm[pair_idx]) - g * kappa_2npm[pair_idx] * math.sinh(h * kappa_2npm[pair_idx])
            a_2npm[pair_idx] = 0.5 * theta_a(a[n], b[n], a[n], b[n], a[m], b[m], h)
            b_2npm[pair_idx] = 0.5 * theta_b(a[n], b[n], a[n], b[n], a[m], b[m], h)
            g_2npm[pair_idx] = lambda3(
                omega1[n], kx[n], ky[n], kappa[n],
                omega1[n], kx[n], ky[n], kappa[n],
                omega1[m], kx[m], ky[m], kappa[m],
                kappa_2[n], gamma_2[n], g_2[n], f_2[n],
                coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                omega_2npm[pair_idx], alpha_2npm, gamma_2npm, beta_2npm, g, h,
            )
            f_2npm[pair_idx] = gamma3(
                omega1[n], kx[n], ky[n], kappa[n],
                omega1[n], kx[n], ky[n], kappa[n],
                omega1[m], kx[m], ky[m], kappa[m],
                kappa_2[n], gamma_2[n], g_2[n], f_2[n],
                coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                omega_2npm[pair_idx], beta_2npm, g, h,
            )
            mu_2npm[pair_idx] = pi_fn(
                omega1[n], kappa[n], omega1[n], kappa[n], omega1[m], kappa[m],
                gamma_2[n], g_2[n], f_2[n],
                coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                f_2npm[pair_idx], kappa_2npm[pair_idx], g, h,
            )
            pair_idx += 1

    num_triplets = n_comp * (n_comp - 1) * (n_comp - 2) // 6
    omega_npmpp = np.zeros(num_triplets)
    kx_npmpp = np.zeros(num_triplets)
    ky_npmpp = np.zeros(num_triplets)
    kappa_npmpp = np.zeros(num_triplets)
    a_npmpp = np.zeros(num_triplets)
    b_npmpp = np.zeros(num_triplets)
    f_npmpp = np.zeros(num_triplets)
    g_npmpp = np.zeros(num_triplets)
    mu_npmpp = np.zeros(num_triplets)

    c3 = 0
    for n in range(n_comp):
        for m in range(n + 1, n_comp):
            idx_sum_nm = m_nm_row_odd[n, m] - 1
            for p in range(m + 1, n_comp):
                idx_sum_np = m_nm_col_odd[n, p] - 1
                idx_sum_mp = m_nm_col_odd[m, p] - 1
                omega_npmpp[c3] = omega1[n] + omega1[m] + omega1[p]
                kx_npmpp[c3] = kx[n] + kx[m] + kx[p]
                ky_npmpp[c3] = ky[n] + ky[m] + ky[p]
                kappa_npmpp[c3] = math.hypot(kx_npmpp[c3], ky_npmpp[c3])
                alpha_npmpp = omega_npmpp[c3] * math.cosh(h * kappa_npmpp[c3])
                gamma_npmpp = kappa_npmpp[c3] * math.sinh(h * kappa_npmpp[c3])
                beta_npmpp = omega_npmpp[c3] ** 2 * math.cosh(h * kappa_npmpp[c3]) - g * kappa_npmpp[c3] * math.sinh(h * kappa_npmpp[c3])
                a_npmpp[c3] = 0.5 * theta_a(a[n], b[n], a[m], b[m], a[p], b[p], h)
                b_npmpp[c3] = 0.5 * theta_b(a[n], b[n], a[m], b[m], a[p], b[p], h)
                g_npmpp[c3] = lambda3(
                    omega1[n], kx[n], ky[n], kappa[n],
                    omega1[m], kx[m], ky[m], kappa[m],
                    omega1[p], kx[p], ky[p], kappa[p],
                    coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                    coeffs["kappa_npm"][idx_sum_np], coeffs["gamma_npm"][idx_sum_np], coeffs["G_npm"][idx_sum_np], coeffs["F_npm"][idx_sum_np],
                    coeffs["kappa_npm"][idx_sum_mp], coeffs["gamma_npm"][idx_sum_mp], coeffs["G_npm"][idx_sum_mp], coeffs["F_npm"][idx_sum_mp],
                    omega_npmpp[c3], alpha_npmpp, gamma_npmpp, beta_npmpp, g, h,
                )
                f_npmpp[c3] = gamma3(
                    omega1[n], kx[n], ky[n], kappa[n],
                    omega1[m], kx[m], ky[m], kappa[m],
                    omega1[p], kx[p], ky[p], kappa[p],
                    coeffs["kappa_npm"][idx_sum_nm], coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                    coeffs["kappa_npm"][idx_sum_np], coeffs["gamma_npm"][idx_sum_np], coeffs["G_npm"][idx_sum_np], coeffs["F_npm"][idx_sum_np],
                    coeffs["kappa_npm"][idx_sum_mp], coeffs["gamma_npm"][idx_sum_mp], coeffs["G_npm"][idx_sum_mp], coeffs["F_npm"][idx_sum_mp],
                    omega_npmpp[c3], beta_npmpp, g, h,
                )
                mu_npmpp[c3] = pi_fn(
                    omega1[n], kappa[n], omega1[m], kappa[m], omega1[p], kappa[p],
                    coeffs["gamma_npm"][idx_sum_nm], coeffs["G_npm"][idx_sum_nm], coeffs["F_npm"][idx_sum_nm],
                    coeffs["gamma_npm"][idx_sum_np], coeffs["G_npm"][idx_sum_np], coeffs["F_npm"][idx_sum_np],
                    coeffs["gamma_npm"][idx_sum_mp], coeffs["G_npm"][idx_sum_mp], coeffs["F_npm"][idx_sum_mp],
                    f_npmpp[c3], kappa_npmpp[c3], g, h,
                )
                c3 += 1

    omega_npm_corr = coeffs["omega_npm"].copy()
    pair_idx = 0
    for n in range(n_comp):
        for m in range(n + 1, n_comp):
            idx_sum = 2 * pair_idx
            idx_diff = idx_sum + 1
            omega_npm_corr[idx_sum] = omega[n] + omega[m]
            omega_npm_corr[idx_diff] = omega[n] - omega[m]
            omega_np2m[pair_idx] = omega[n] + 2.0 * omega[m]
            omega_2npm[pair_idx] = 2.0 * omega[n] + omega[m]
            pair_idx += 1
    c3 = 0
    for n in range(n_comp):
        for m in range(n + 1, n_comp):
            for p in range(m + 1, n_comp):
                omega_npmpp[c3] = omega[n] + omega[m] + omega[p]
                c3 += 1

    return {
        "omega": omega,
        "muStar": mu_star,
        "F13": f13,
        "A_3": a_3,
        "B_3": b_3,
        "F_3": f_3,
        "G_3": g_3,
        "mu_3": mu_3,
        "kappa_3": 3.0 * kappa,
        "A_np2m": a_np2m,
        "B_np2m": b_np2m,
        "F_np2m": f_np2m,
        "G_np2m": g_np2m,
        "mu_np2m": mu_np2m,
        "kappa_np2m": kappa_np2m,
        "omega_np2m": omega_np2m,
        "kx_np2m": kx_np2m,
        "ky_np2m": ky_np2m,
        "A_2npm": a_2npm,
        "B_2npm": b_2npm,
        "F_2npm": f_2npm,
        "G_2npm": g_2npm,
        "mu_2npm": mu_2npm,
        "kappa_2npm": kappa_2npm,
        "omega_2npm": omega_2npm,
        "kx_2npm": kx_2npm,
        "ky_2npm": ky_2npm,
        "A_npmpp": a_npmpp,
        "B_npmpp": b_npmpp,
        "F_npmpp": f_npmpp,
        "G_npmpp": g_npmpp,
        "mu_npmpp": mu_npmpp,
        "kappa_npmpp": kappa_npmpp,
        "omega_npmpp": omega_npmpp,
        "kx_npmpp": kx_npmpp,
        "ky_npmpp": ky_npmpp,
        "gamma_2": gamma_2,
        "omega_npm": omega_npm_corr,
    }


def upsilon_nm(omega1n: float, knx: float, kny: float, kappan: float, omega1m: float, kmx: float, kmy: float, kappam: float, fnpm: float, fnmm: float, gnpm: float, gnmm: float, kappanpm: float, kappanmm: float, g: float, h: float) -> float:
    knkm = knx * kmx + kny * kmy
    return g / (4.0 * omega1n * omega1m * math.cosh(h * kappan)) * (omega1m * (kappan * kappan - kappam * kappam) - omega1n * knkm) + (gnpm + gnmm) / (4.0 * h * omega1n * omega1n * omega1m * math.cosh(h * kappan)) * (g * g * knkm + omega1m**3 * omega1n) - 1.0 / (4.0 * h * math.cosh(h * kappan)) * (fnpm * kappanpm * math.sinh(h * kappanpm) + fnmm * kappanmm * math.sinh(h * kappanmm)) + g * fnpm * math.cosh(h * kappanpm) / (4.0 * h * omega1n * omega1n * omega1m * math.cosh(h * kappan)) * ((omega1n + omega1m) * (knkm + kappam * kappam) - omega1m * kappanpm * kappanpm) + g * fnmm * math.cosh(h * kappanmm) / (4.0 * h * omega1n * omega1n * omega1m * math.cosh(h * kappan)) * ((omega1n - omega1m) * (knkm - kappam * kappam) - omega1m * kappanmm * kappanmm)


def xi_nm(omega1n: float, kappan: float, omega1m: float, gnpm: float, fnpm: float, gamma_npm: float, gnmm: float, fnmm: float, gamma_nmm: float, h: float, g: float) -> float:
    return 1.0 / (2.0 * h) * (omega1m * (gnpm - gnmm) + fnpm * gamma_npm + fnmm * gamma_nmm - g * h * kappan * kappan / (2.0 * omega1n))


def theta_a(an: float, bn: float, am: float, bm: float, ap: float, bp: float, h: float) -> float:
    return (an * am * ap - bn * bm * ap - bn * am * bp - an * bm * bp) / (h * h)


def theta_b(an: float, bn: float, am: float, bm: float, ap: float, bp: float, h: float) -> float:
    return (bn * am * ap + an * bm * ap + an * am * bp - bn * bm * bp) / (h * h)


def omega_nm_fn(omega1n: float, knx: float, kny: float, omega1m: float, kmx: float, kmy: float, kappam: float, fnpm: float, fnmm: float, gnpm: float, gnmm: float, kappanpm: float, kappanmm: float, g: float, h: float) -> float:
    knkm = knx * kmx + kny * kmy
    return 1.0 / (kappam * kappam) * (((2.0 * omega1m * omega1m + omega1n * omega1n) / (4.0 * omega1n * omega1m)) * knkm + 0.25 * kappam * kappam) + (gnpm + gnmm) / (kappam * kappam) * (g * knkm / (4.0 * h * omega1n * omega1m) - omega1m * omega1m / (4.0 * g * h)) + omega1n / (4.0 * g * h * kappam * kappam) * (fnpm * kappanpm * math.sinh(h * kappanpm) + fnmm * kappanmm * math.sinh(h * kappanmm)) - fnpm * math.cosh(h * kappanpm) / (4.0 * h * omega1n * omega1m * kappam * kappam) * ((omega1n - omega1m) * (kappam * kappam + knkm) + omega1m * kappanpm * kappanpm) + fnmm * math.cosh(h * kappanmm) / (4.0 * h * omega1n * omega1m * kappam * kappam) * ((omega1n + omega1m) * (kappam * kappam - knkm) - omega1m * kappanmm * kappanmm)


def lambda3(
    omega1n: float, knx: float, kny: float, kappan: float,
    omega1m: float, kmx: float, kmy: float, kappam: float,
    omega1p: float, kpx: float, kpy: float, kappap: float,
    kappanpm: float, gammanpm: float, gnpm: float, fnpm: float,
    kappanpp: float, gammanpp: float, gnpp: float, fnpp: float,
    kappampp: float, gammampp: float, gmpp: float, fmpp: float,
    omega_npmpp: float, alpha_npmpp: float, gamma_npmpp: float, beta_npmpp: float,
    g: float, h: float,
) -> float:
    knkm = knx * kmx + kny * kmy
    knkp = knx * kpx + kny * kpy
    kmkp = kmx * kpx + kmy * kpy
    return (
        h * h / (4.0 * beta_npmpp) * (
            alpha_npmpp * (
                omega1n * (knkm + knkp + kappan * kappan)
                + omega1m * (knkm + kmkp + kappam * kappam)
                + omega1p * (knkp + kmkp + kappap * kappap)
            )
            + gamma_npmpp * (
                g / omega1n * (omega1m * knkm + omega1p * knkp - omega_npmpp * kappan * kappan)
                + g / omega1m * (omega1n * knkm + omega1p * kmkp - omega_npmpp * kappam * kappam)
                + g / omega1p * (omega1n * knkp + omega1m * kmkp - omega_npmpp * kappap * kappap)
            )
        )
        - h * fnpm / (2.0 * beta_npmpp) * (
            alpha_npmpp * math.cosh(h * kappanpm) * (knkp + kmkp + kappanpm * kappanpm)
            + gamma_npmpp * (g / omega1p * (knkp + kmkp) * math.cosh(h * kappanpm) - gammanpm * omega_npmpp)
        )
        - h * fnpp / (2.0 * beta_npmpp) * (
            alpha_npmpp * math.cosh(h * kappanpp) * (knkm + kmkp + kappanpp * kappanpp)
            + gamma_npmpp * (g / omega1m * (knkm + kmkp) * math.cosh(h * kappanpp) - gammanpp * omega_npmpp)
        )
        - h * fmpp / (2.0 * beta_npmpp) * (
            alpha_npmpp * math.cosh(h * kappampp) * (knkm + knkp + kappampp * kappampp)
            + gamma_npmpp * (g / omega1n * (knkm + knkp) * math.cosh(h * kappampp) - gammampp * omega_npmpp)
        )
        + h * gnpm / (2.0 * beta_npmpp) * (alpha_npmpp * g / omega1p * (knkp + kmkp + kappap * kappap) - gamma_npmpp * omega1p * omega1p)
        + h * gnpp / (2.0 * beta_npmpp) * (alpha_npmpp * g / omega1m * (knkm + kmkp + kappam * kappam) - gamma_npmpp * omega1m * omega1m)
        + h * gmpp / (2.0 * beta_npmpp) * (alpha_npmpp * g / omega1n * (knkm + knkp + kappan * kappan) - gamma_npmpp * omega1n * omega1n)
    )


def gamma3(
    omega1n: float, knx: float, kny: float, kappan: float,
    omega1m: float, kmx: float, kmy: float, kappam: float,
    omega1p: float, kpx: float, kpy: float, kappap: float,
    kappanpm: float, gammanpm: float, gnpm: float, fnpm: float,
    kappanpp: float, gammanpp: float, gnpp: float, fnpp: float,
    kappampp: float, gammampp: float, gmpp: float, fmpp: float,
    omega_npmpp: float, beta_npmpp: float,
    g: float, h: float,
) -> float:
    knkm = knx * kmx + kny * kmy
    knkp = knx * kpx + kny * kpy
    kmkp = kmx * kpx + kmy * kpy
    return (
        -g * h * h / (4.0 * beta_npmpp) * (
            omega1n * (knkm + knkp + kappan * kappan)
            + omega1m * (knkm + kmkp + kappam * kappam)
            + omega1p * (knkp + kmkp + kappap * kappap)
            + omega_npmpp / omega1n * (omega1m * knkm + omega1p * knkp - omega_npmpp * kappan * kappan)
            + omega_npmpp / omega1m * (omega1n * knkm + omega1p * kmkp - omega_npmpp * kappam * kappam)
            + omega_npmpp / omega1p * (omega1n * knkp + omega1m * kmkp - omega_npmpp * kappap * kappap)
        )
        + h * fnpm / (2.0 * beta_npmpp) * (
            g * math.cosh(h * kappanpm) * ((knkp + kmkp + kappanpm * kappanpm) + omega_npmpp / omega1p * (knkp + kmkp))
            - gammanpm * omega_npmpp * omega_npmpp
        )
        + h * fnpp / (2.0 * beta_npmpp) * (
            g * math.cosh(h * kappanpp) * ((knkm + kmkp + kappanpp * kappanpp) + omega_npmpp / omega1m * (knkm + kmkp))
            - gammanpp * omega_npmpp * omega_npmpp
        )
        + h * fmpp / (2.0 * beta_npmpp) * (
            g * math.cosh(h * kappampp) * ((knkm + knkp + kappampp * kappampp) + omega_npmpp / omega1n * (knkm + knkp))
            - gammampp * omega_npmpp * omega_npmpp
        )
        + h * gnpm / (2.0 * beta_npmpp) * (omega1p * omega1p * omega_npmpp - g * g / omega1p * (knkp + kmkp + kappap * kappap))
        + h * gnpp / (2.0 * beta_npmpp) * (omega1m * omega1m * omega_npmpp - g * g / omega1m * (knkm + kmkp + kappam * kappam))
        + h * gmpp / (2.0 * beta_npmpp) * (omega1n * omega1n * omega_npmpp - g * g / omega1n * (knkm + knkp + kappan * kappan))
    )


def pi_fn(
    omega1n: float, kappan: float,
    omega1m: float, kappam: float,
    omega1p: float, kappap: float,
    gamma_npm: float, gnpm: float, fnpm: float,
    gamma_npp: float, gnpp: float, fnpp: float,
    gamma_mpp: float, gmpp: float, fmpp: float,
    fnpmpp: float, kappa_npmpp: float,
    g: float, h: float,
) -> float:
    return (
        fnpmpp * math.cosh(h * kappa_npmpp)
        - g * h * h / 4.0 * (kappan * kappan / omega1n + kappam * kappam / omega1m + kappap * kappap / omega1p)
        - h / 2.0 * (omega1n * gmpp + omega1m * gnpp + omega1p * gnpm)
        + h / 2.0 * (fnpm * gamma_npm + fnpp * gamma_npp + fmpp * gamma_mpp)
    )
