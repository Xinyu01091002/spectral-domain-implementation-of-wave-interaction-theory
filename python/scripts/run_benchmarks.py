from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    config = json.loads((repo_root / "cross_language_comparison" / "cases" / "index.json").read_text(encoding="utf-8"))
    cases = [str((repo_root / case).resolve()) for case in config["cases"]]
    defaults = config.get("python_defaults", {})

    cmd = [
        sys.executable,
        "-m",
        "mf12_python.cli",
        "benchmark",
        *cases,
        "--repeats",
        str(defaults.get("repeats", 3)),
    ]
    if defaults.get("warmup", False):
        cmd.append("--warmup")

    env = dict(os.environ)
    env["PYTHONPATH"] = str((repo_root / "python" / "src").resolve())
    subprocess.run(cmd, cwd=repo_root, env=env, check=True)


if __name__ == "__main__":
    main()
