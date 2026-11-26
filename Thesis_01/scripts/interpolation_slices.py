#!/usr/bin/env python3
"""Slice figures of the interpolated reconstructions of S_0.

Light blue: the exact S_-eps and S_+eps; grey: the exact S_0; red: the
reconstruction (intersection of its mesh with the plane, from
interpolation_results/slices/*.csv written by
`isosurface_interpolation_experiments interpolation_results --slices`).

    python3 Thesis_01/scripts/interpolation_slices.py        (from the project root)
"""
import csv
import os
import sys

from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact_slices as ex  # noqa: E402

ex.CELLS = 480
SLICES = os.path.join(ex.ROOT, "interpolation_results", "slices")
OUT = os.path.join(ex.ROOT, "Thesis_01", "pictures", "report", "interp")
LIGHT = (150, 185, 235)

CASES = {  # tag: (F, G, eps, plane)
    "A1_1": ("x**2+y**2-z**2+x*y*z", "1", 1e-2, "x"),
    "A1_x": ("x**2+y**2-z**2+x*y*z", "x", 1e-2, "x"),
    "A2a_x": ("x**2+y**2+z**3", "x", 1e-2, "y"),
    "E6a_x": ("x**3+y**3+z**3", "x", 1e-2, "z"),
}
METHODS = ["vertex_mid_from_minus", "vertex_mid_from_plus", "sdf_average", "exact_vertex_mid", "normal_mid",
           "hybrid_cone", "singular_cone", "direct_S0"]
# Close-ups of P (window [-ZOOM, ZOOM]^2) for the extraction of S_0 itself.
ZOOM = 0.04
ZOOM_METHODS = ["direct_S0", "singular_cone"]


def to_grid(u, v):
    return (u + ex.HALF) / (2 * ex.HALF) * ex.CELLS, (v + ex.HALF) / (2 * ex.HALF) * ex.CELLS


def panels(tag, ftext, g, eps, plane, methods, suffix):
    f = ex.compile_poly(ftext)
    gf = ex.compile_poly(g)
    s0 = ex.contour(ex.sample(f, plane))
    sm = ex.contour(ex.sample(lambda x, y, z: f(x, y, z) - eps * gf(x, y, z), plane))
    sp = ex.contour(ex.sample(lambda x, y, z: f(x, y, z) + eps * gf(x, y, z), plane))
    for m in methods:
        segs = []
        with open(os.path.join(SLICES, f"{tag}_{m}.csv")) as fh:
            for r in csv.DictReader(fh):
                u0, v0, u1, v1 = (float(r[k]) for k in ("u0", "v0", "u1", "v1"))
                if max(abs(u0), abs(v0), abs(u1), abs(v1)) <= ex.HALF * 1.05:
                    segs.append((to_grid(u0, v0), to_grid(u1, v1)))
        layers = [(s0, ex.GREY, 2), (segs, ex.RED, 3)]
        if not suffix:  # the perturbed surfaces only in the overview panels
            layers = [(sm, LIGHT, 2), (sp, LIGHT, 2)] + layers
        ex.draw_panel(os.path.join(OUT, f"interp_{tag}_{m}{suffix}.png"), layers, [])


def main():
    os.makedirs(OUT, exist_ok=True)
    for tag, (ftext, g, eps, plane) in CASES.items():
        ex.HALF = 0.15
        panels(tag, ftext, g, eps, plane, METHODS, "")
        ex.HALF = ZOOM
        panels(tag, ftext, g, eps, plane, ZOOM_METHODS, "_zoom")
        print("done", tag)


if __name__ == "__main__":
    main()
