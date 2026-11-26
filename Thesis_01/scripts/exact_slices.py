#!/usr/bin/env python3
"""Exact slices of S_k = {F + kG = 0} through the singular point P = 0.

Draws the panels in pictures/report/exact/ directly from the polynomials
(marching squares on a fine grid, no voxel mesh involved):
    grey  S_0,   blue  S_k,   red  singular points of S_k,   cross  the origin.
Singular points are read from the stationary-point search of the experiments
(experiments_thesis/<surface>/res128/<G>/stationary_points_certified.csv).

Only the Python standard library and Pillow are needed.

    python3 Thesis_01/scripts/exact_slices.py            (from the project root)
"""
import csv
import math
import os
import sys

from PIL import Image, ImageDraw

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "Thesis_01", "pictures", "report", "exact")
DATA = os.path.join(ROOT, "experiments_thesis")

HALF = 0.25      # window [-HALF, HALF]^2 on the slice plane
K = 0.01         # |k| of the panels
CELLS = 560      # marching-squares cells per axis
SIZE = 420       # output pixels
SS = 4           # supersampling factor for anti-aliasing

WHITE, BORDER, AXIS = (255, 255, 255), (160, 160, 160), (228, 228, 228)
GREY, BLACK, BLUE, RED = (200, 200, 200), (0, 0, 0), (30, 80, 200), (220, 30, 30)

# id: (experiment folder, F, slice plane as the fixed coordinate, panel G list)
SURFACES = {
    "A1":  ("A1_x2+y2-z2+xyz", "x**2+y**2-z**2+x*y*z", "x", ["1", "x+y", "x^2+y^2-z^2"]),
    "A2a": ("A2_x2+y2+z3", "x**2+y**2+z**3", "y", ["1", "z", "x"]),
    "A2b": ("A2_x2-y2+z3", "x**2-y**2+z**3", "z", ["1", "x+y", "z"]),
    "A2c": ("A2_x3+y3+z2+x2", "x**3+y**3+z**2+x**2", "x", ["1", "y", "x+z"]),
    "D4a": ("D4_x2+y3+z3", "x**2+y**3+z**3", "x", ["1", "y", "x^2+y^2-z^2"]),
    "D4b": ("D4_x2+y3-3yz2", "x**2+y**3-3*y*z**2", "x", ["1", "y", "z"]),
    "E6a": ("E6_x3+y3+z3", "x**3+y**3+z**3", "z", ["1", "x+y", "x^2+y^2-z^2"]),
    "E6b": ("E6_x3+y3-z3", "x**3+y**3-z**3", "z", ["1", "x+y", "x^2+y^2-z^2"]),
    "E6c": ("E6_x3+y3+z3+xyz", "x**3+y**3+z**3+x*y*z", "z", ["1", "x", "x^2+y^2-z^2"]),
    "E6d": ("E6_x3+2y3+3z3+2xyz", "x**3+2*y**3+3*z**3+2*x*y*z", "z", ["1", "x", "x^2+y^2-z^2"]),
}
FILE_TAG = {"1": "1", "x": "x", "y": "y", "z": "z", "x+y": "xpy", "x+z": "xpz", "x^2+y^2-z^2": "x2py2mz2"}


def compile_poly(text):
    code = compile(text.replace("^", "**"), text, "eval")
    return lambda x, y, z: eval(code, {}, {"x": x, "y": y, "z": z})


def to_xyz(plane, u, v):
    """Slice coordinates (u horizontal, v vertical) -> (x, y, z)."""
    if plane == "x":
        return 0.0, u, v
    if plane == "y":
        return u, 0.0, v
    return u, v, 0.0


def plane_coords(plane, p):
    x, y, z = p
    return {"x": (x, (y, z)), "y": (y, (x, z)), "z": (z, (x, y))}[plane]


def sample(fn, plane):
    n = CELLS + 1
    step = 2 * HALF / CELLS
    return [[fn(*to_xyz(plane, -HALF + i * step, -HALF + j * step)) for i in range(n)] for j in range(n)]


def contour(values):
    """Marching squares: zero-level segments in grid-index coordinates."""
    segs = []
    for j in range(CELLS):
        r0, r1 = values[j], values[j + 1]
        for i in range(CELLS):
            a, b, c, d = r0[i], r0[i + 1], r1[i + 1], r1[i]  # corners counter-clockwise
            pts = []
            for (p, q, fp, fq) in (((i, j), (i + 1, j), a, b), ((i + 1, j), (i + 1, j + 1), b, c),
                                   ((i + 1, j + 1), (i, j + 1), c, d), ((i, j + 1), (i, j), d, a)):
                if (fp < 0) != (fq < 0):
                    t = fp / (fp - fq)
                    pts.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
            if len(pts) == 2:
                segs.append((pts[0], pts[1]))
            elif len(pts) == 4:
                segs.append((pts[0], pts[1]))
                segs.append((pts[2], pts[3]))
    return segs


def pix(gx, gy):
    s = SIZE * SS / CELLS
    return gx * s, (CELLS - gy) * s  # vertical axis upwards


def world_pix(u, v):
    return pix((u + HALF) / (2 * HALF) * CELLS, (v + HALF) / (2 * HALF) * CELLS)


def singular_points(folder, g, k, plane):
    gdir = g.replace("*", "x")
    path = os.path.join(DATA, folder, "res128", gdir, "stationary_points_certified.csv")
    out = []
    if not os.path.exists(path):
        print("  missing", path, file=sys.stderr)
        return out
    with open(path) as f:
        for r in csv.DictReader(f):
            if abs(float(r["k"]) - k) > 1e-12 or r["on_zero_level"] != "1":
                continue
            p = (float(r["x"]), float(r["y"]), float(r["z"]))
            off, (u, v) = plane_coords(plane, p)
            if abs(off) < 1e-6 and abs(u) <= HALF and abs(v) <= HALF:
                out.append((u, v))
    return out


def draw_panel(path, layers, points):
    im = Image.new("RGB", (SIZE * SS, SIZE * SS), WHITE)
    d = ImageDraw.Draw(im)
    c = SIZE * SS / 2
    d.line([(0, c), (SIZE * SS, c)], fill=AXIS, width=SS)
    d.line([(c, 0), (c, SIZE * SS)], fill=AXIS, width=SS)
    for segs, colour, width in layers:
        for (p, q) in segs:
            d.line([pix(*p), pix(*q)], fill=colour, width=width * SS)
    a = 5 * SS
    d.line([(c - a, c), (c + a, c)], fill=BLACK, width=SS)
    d.line([(c, c - a), (c, c + a)], fill=BLACK, width=SS)
    for (u, v) in points:
        x, y = world_pix(u, v)
        r = 5 * SS
        d.ellipse([x - r, y - r, x + r, y + r], fill=RED)
    im = im.resize((SIZE, SIZE), Image.LANCZOS)
    ImageDraw.Draw(im).rectangle([0, 0, SIZE - 1, SIZE - 1], outline=BORDER)
    im.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    only = set(sys.argv[1:])
    for sid, (folder, ftext, plane, glist) in SURFACES.items():
        if only and sid not in only:
            continue
        f = compile_poly(ftext)
        s0 = contour(sample(f, plane))
        draw_panel(os.path.join(OUT, f"exact_{sid}_k0.png"), [(s0, BLACK, 2)], [(0.0, 0.0)])
        for g in glist:
            gf = compile_poly(g)
            for k, tag in ((-K, "neg"), (K, "pos")):
                h = lambda x, y, z, k=k: f(x, y, z) + k * gf(x, y, z)
                sk = contour(sample(h, plane))
                pts = singular_points(folder, g, k, plane)
                draw_panel(os.path.join(OUT, f"exact_{sid}_{FILE_TAG[g]}_{tag}.png"),
                           [(s0, GREY, 2), (sk, BLUE, 3)], pts)
        print("done", sid, flush=True)


if __name__ == "__main__":
    main()
