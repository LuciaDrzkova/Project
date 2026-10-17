#!/usr/bin/env python3
"""Run broad catalogue once, then selected difficult families across resolutions.

Run from the repository root after building isosurface_singularity_diagnostics:
  python3 singularity_diagnostics_update/run_resolution_series.py

Defaults:
  * broad 38-surface catalogue + phase regression at resolution 48;
  * eight high-value families at resolutions 48, 64, 96, 128;
  * grid-phase checks on those families at each focused resolution;
  * geometry-interpolation checks for Mixed cubic and the plane example at each resolution;
  * five k values for the high-resolution subset to contain runtime;
  * aggregate CSVs for base meshes, phase variations, and interpolation checks.

Use --exe, --resolutions, --broad-resolution, --out, --focus, and --k-values
for a custom run. This wrapper aggregates rows but never treats topology defects
as mathematical evidence by themselves.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import subprocess
import sys

DEFAULT_FOCUS = [
    "plane_x3_ky", "quartic_cusp_x", "mixed_cusp_x", "mixed_cusp_y",
    "mixed_cubic", "mixed_order_z", "whitney_umbrella", "cubic_family",
]


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_aggregated(path: Path, tagged_runs: list[tuple[str, Path]]):
    rows = []
    base_fields = None
    for tag, source in tagged_runs:
        if not source.exists():
            print(f"WARNING: missing output file, skipped: {source}", file=sys.stderr)
            continue
        with source.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if base_fields is None:
                base_fields = reader.fieldnames or []
            for row in reader:
                row["resolution_run"] = tag
                rows.append(row)
    if base_fields is None:
        print(f"No source CSVs found for {path}", file=sys.stderr)
        return
    fields = ["resolution_run"] + [f for f in base_fields if f != "resolution_run"]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {path} ({len(rows)} rows)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", default="./build/isosurface_singularity_diagnostics",
                        help="diagnostics executable path (relative to repository root by default)")
    parser.add_argument("--resolutions", default="48,64,96,128",
                        help="comma-separated resolutions for focused families")
    parser.add_argument("--broad-resolution", type=int, default=48,
                        help="resolution for the broad catalogue / grid-phase stage")
    parser.add_argument("--out", default="singularity_resolution_series",
                        help="output directory")
    parser.add_argument("--focus", default=",".join(DEFAULT_FOCUS),
                        help="comma-separated surface IDs for the high-resolution subset")
    parser.add_argument("--k-values", default="-0.1,-0.01,0,0.01,0.1",
                        help="comma-separated k values for the focused resolution runs")
    parser.add_argument("--skip-broad", action="store_true", help="skip the broad catalogue stage")
    args = parser.parse_args()

    root = Path.cwd()
    exe = Path(args.exe)
    if not exe.is_absolute():
        exe = (root / exe).resolve()
    if not exe.is_file():
        parser.error(f"executable not found: {exe}. Build the diagnostics target first.")

    out = (root / args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    runs: list[tuple[str, Path]] = []
    phase_runs: list[tuple[str, Path]] = []
    phase_variation_runs: list[tuple[str, Path]] = []
    interpolation_runs: list[tuple[str, Path]] = []
    assessment_runs: list[tuple[str, Path]] = []

    if not args.skip_broad:
        tag = f"broad_N{args.broad_resolution}"
        run_dir = out / tag
        command = [str(exe), "--resolution", str(args.broad_resolution), "--out", str(run_dir)]
        print("+", " ".join(command))
        subprocess.run(command, cwd=root, check=True)
        runs.append((tag, run_dir / "summary.csv"))
        phase_runs.append((tag, run_dir / "phase_summary.csv"))
        phase_variation_runs.append((tag, run_dir / "phase_variation.csv"))
        interpolation_runs.append((tag, run_dir / "interpolation_validation.csv"))
        assessment_runs.append((tag, run_dir / "singularity_assessment.csv"))

    resolutions = [int(x) for x in args.resolutions.split(",") if x.strip()]
    focus = [x.strip() for x in args.focus.split(",") if x.strip()]
    if not focus:
        parser.error("--focus must contain at least one surface ID")

    for res in resolutions:
        tag = f"focus_N{res}"
        run_dir = out / tag
        command = [str(exe), "--resolution", str(res), "--out", str(run_dir),
                   "--k-values", args.k_values,
                   "--interpolation-surface", "mixed_cubic",
                   "--interpolation-surface", "plane_x3_ky",
                   "--interpolation-magnitudes", "0.1,0.01"]
        for surface_id in focus:
            command.extend(["--surface", surface_id])
        print("+", " ".join(command))
        subprocess.run(command, cwd=root, check=True)
        runs.append((tag, run_dir / "summary.csv"))
        phase_runs.append((tag, run_dir / "phase_summary.csv"))
        phase_variation_runs.append((tag, run_dir / "phase_variation.csv"))
        interpolation_runs.append((tag, run_dir / "interpolation_validation.csv"))
        assessment_runs.append((tag, run_dir / "singularity_assessment.csv"))

    write_aggregated(out / "all_resolutions.csv", runs)
    write_aggregated(out / "all_phase_resolutions.csv", phase_runs)
    write_aggregated(out / "all_phase_variation.csv", phase_variation_runs)
    write_aggregated(out / "all_interpolation_validation.csv", interpolation_runs)
    write_aggregated(out / "all_singularity_assessments.csv", assessment_runs)
    print("Done. Keep resolution_run in the merged CSVs to distinguish the broad run from focused runs.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
