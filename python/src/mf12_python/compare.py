from __future__ import annotations

from typing import Any

import numpy as np


def compare_fields(candidate: np.ndarray, reference: np.ndarray, prefix: str) -> dict[str, float]:
    cand = np.asarray(candidate, dtype=np.float64)
    ref = np.asarray(reference, dtype=np.float64)
    diff = cand - ref
    denom = np.linalg.norm(ref.ravel())
    return {
        f"{prefix}_max_abs_err": float(np.max(np.abs(diff))),
        f"{prefix}_rms_err": float(np.sqrt(np.mean(diff * diff))),
        f"{prefix}_relative_l2_err": float(np.linalg.norm(diff.ravel()) / max(denom, np.finfo(float).eps)),
    }


def compare_result_to_reference(result: dict[str, Any], case: dict[str, Any]) -> dict[str, float]:
    reference = case["reference"]["loaded_arrays"]
    metrics: dict[str, float] = {}
    metrics.update(compare_fields(result["eta"], reference["eta"], "eta"))
    metrics.update(compare_fields(result["phi"], reference["phi"], "phi"))
    if "kinematics" in result:
        for name, values in result["kinematics"].items():
            if name in reference:
                metrics.update(compare_fields(values, reference[name], name))

    ref_meta = case["reference"].get("loaded_metadata", {})
    matlab_runtime = ref_meta.get("runtime", {})
    py_runtime = result.get("runtime", {})
    matlab_total = float(matlab_runtime.get("mean_total_s", 0.0))
    matlab_recon = float(matlab_runtime.get("mean_reconstruction_s", 0.0))
    py_total = float(py_runtime.get("mean_total_s", 0.0))
    py_recon = float(py_runtime.get("mean_reconstruction_s", 0.0))

    if matlab_total > 0.0 and py_total > 0.0:
        metrics["speedup_vs_matlab_total"] = matlab_total / py_total
    if matlab_recon > 0.0 and py_recon > 0.0:
        metrics["speedup_vs_matlab_reconstruction"] = matlab_recon / py_recon
    return metrics


def tolerances_pass(metrics: dict[str, float], tolerances: dict[str, float]) -> bool:
    for key, limit in tolerances.items():
        if key in metrics and metrics[key] > limit:
            return False
    return True
