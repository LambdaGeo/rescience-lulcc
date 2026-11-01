# SPDX-License-Identifier: GPL-3.0-or-later
"""Platform information, as ReScience C asks (https://rescience.github.io/platform/).

Usage: python scripts/platform_info.py [--out results/platform.txt]

Prints the operating system, the Python interpreter, the machine, and the
versions of the packages the replication imports, and writes the same text to
--out. Run it with the interpreter that produced the results (make platform).
The R environment of the reference is in lulcc-goldens/goldens/environment.json.
"""
from __future__ import annotations

import argparse
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

PACKAGES = ["dissmodel", "disslucc", "numpy", "scipy", "pandas", "statsmodels",
            "rasterio", "matplotlib"]


def cpu() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def os_name() -> str:
    try:
        rel = dict(line.split("=", 1) for line in
                   Path("/etc/os-release").read_text().splitlines() if "=" in line)
        return rel.get("PRETTY_NAME", "").strip('"') or platform.platform()
    except OSError:
        return platform.platform()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("results/platform.txt"))
    args = ap.parse_args()

    lines = [
        f"Platform: {sys.platform} ({platform.platform()})",
        f"Operating system: {os_name()}",
        f"Machine and architecture: {platform.machine()} {' '.join(platform.architecture())}",
        f"CPU: {cpu()}",
        f"Python: {sys.version.splitlines()[0]}",
    ]
    for name in PACKAGES:
        try:
            lines.append(f"{name}: {version(name)}")
        except PackageNotFoundError:
            lines.append(f"{name}: not installed")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
