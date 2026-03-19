from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


def load_case(case_dir: str | Path) -> dict[str, Any]:
    case_path = Path(case_dir)
    manifest = json.loads(case_path.joinpath("case.json").read_text(encoding="utf-8"))

    arrays: dict[str, np.ndarray] = {}
    for name, rel_path in manifest.get("arrays", {}).items():
        arrays[name] = load_csv(case_path / rel_path)

    reference = manifest.get("reference", {})
    if reference:
        ref_arrays: dict[str, np.ndarray] = {}
        for name, rel_path in reference.get("arrays", {}).items():
            ref_arrays[name] = load_csv(case_path / rel_path)
        manifest["reference"]["loaded_arrays"] = ref_arrays
        meta_rel = reference.get("metadata")
        if meta_rel:
            manifest["reference"]["loaded_metadata"] = json.loads((case_path / meta_rel).read_text(encoding="utf-8"))

    manifest["case_dir"] = str(case_path.resolve())
    manifest["arrays"] = arrays
    return manifest


def load_csv(path: str | Path) -> np.ndarray:
    arr = np.loadtxt(path, delimiter=",")
    if np.ndim(arr) == 0:
        arr = np.asarray([float(arr)], dtype=np.float64)
    return np.asarray(arr, dtype=np.float64)


def save_result_bundle(output_dir: str | Path, result: dict[str, Any]) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    np.savetxt(out / "eta.csv", np.asarray(result["eta"]), delimiter=",")
    np.savetxt(out / "phi.csv", np.asarray(result["phi"]), delimiter=",")
    np.savetxt(out / "x.csv", np.asarray(result["x"]).reshape(-1, 1), delimiter=",")
    np.savetxt(out / "y.csv", np.asarray(result["y"]).reshape(-1, 1), delimiter=",")
    if "kinematics" in result:
        for name, values in result["kinematics"].items():
            np.savetxt(out / f"{name}.csv", np.asarray(values), delimiter=",")
    metadata = {
        "runtime": result.get("runtime", {}),
        "metadata": result.get("metadata", {}),
        "comparison": result.get("comparison", {}),
        "case_id": result.get("case_id"),
    }
    out.joinpath("result.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
