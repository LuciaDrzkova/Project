#!/usr/bin/env python3
"""Tables and plot data for the interpolation section of the report.

Reads interpolation_results/interpolation_all.csv (written by
isosurface_interpolation_experiments) and prints LaTeX table rows; writes
pgfplots data files to Thesis_01/data/interp_*.dat.

    python3 Thesis_01/scripts/interpolation_tables.py        (from the project root)
"""
import csv
import math
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CSV = os.path.join(ROOT, "interpolation_results", "interpolation_all.csv")
DATA = os.path.join(ROOT, "Thesis_01", "data")

TEX_F = {
    "A1_x2+y2-z2+xyz": "$x^2+y^2-z^2+xyz$",
    "A2_x2+y2+z3": "$x^2+y^2+z^3$",
    "A2_x2-y2+z3": "$x^2-y^2+z^3$",
    "A2_x3+y3+z2+x2": "$x^3+y^3+z^2+x^2$",
    "D4_x2+y3+z3": "$x^2+y^3+z^3$",
    "D4_x2+y3-3yz2": "$x^2+y^3-3yz^2$",
    "E6_x3+y3+z3": "$x^3+y^3+z^3$",
    "E6_x3+y3-z3": "$x^3+y^3-z^3$",
    "E6_x3+y3+z3+xyz": "$x^3+y^3+z^3+xyz$",
    "E6_x3+2y3+3z3+2xyz": "$x^3+2y^3+3z^3+2xyz$",
}
TEX_G = {"1": "$1$", "x": "$x$", "x+y": "$x+y$", "x^2+y^2-z^2": "$x^2+y^2-z^2$"}
TYPE = {"A1": "$A_1$", "A2": "$A_2$", "D4": "$D_4$", "E6~": "$\\tilde E_6$"}

rows = list(csv.DictReader(open(CSV)))
R = {}
for r in rows:
    R[(r["surface"], r["G"], float(r["eps"]), int(r["resolution"]), r["method"])] = r
surfaces = list(dict.fromkeys(r["surface"] for r in rows))
gs = list(dict.fromkeys(r["G"] for r in rows))
eps = sorted({float(r["eps"]) for r in rows}, reverse=True)
types = {r["surface"]: r["type"] for r in rows}


def val(s, g, e, n, m, col):
    v = R[(s, g, e, n, m)][col]
    return float("inf") if v in ("NA", "") else float(v)  # NA: no part of the surface in Q


def f4(x):
    return "empty" if x == float("inf") else f"{x:.4f}"


def endpoint(s, g, e, n, col):
    return max(val(s, g, e, n, "endpoint_minus", col), val(s, g, e, n, "endpoint_plus", col))


def vmid(s, g, e, n, col):
    return max(val(s, g, e, n, "vertex_mid_from_minus", col), val(s, g, e, n, "vertex_mid_from_plus", col))


def slope(xs, ys):
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if 0 < y < float("inf")]
    if len(pts) < 3:
        return float("nan")
    mx = sum(p[0] for p in pts) / len(pts)
    my = sum(p[1] for p in pts) / len(pts)
    return sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)


def f3(x):
    return "--" if not math.isfinite(x) else f"{x:.3f}"


def e2(x):
    if not math.isfinite(x):
        return "--"
    m, ex = f"{x:.1e}".split("e")
    return f"${m}\\cdot10^{{{int(ex)}}}$"


N = 128
print("% ---- Table 1: errors at eps = 1e-2 and 1e-3, N = 128 (H_Q) ----")
for s in surfaces:
    for g in gs:
        cells = []
        for e in (1e-2, 1e-3):
            cells += [endpoint(s, g, e, N, "hausdorff_Q"), vmid(s, g, e, N, "hausdorff_Q"),
                      val(s, g, e, N, "sdf_average", "hausdorff_Q"), val(s, g, e, N, "exact_vertex_mid", "hausdorff_Q")]
        print(f"{TEX_F[s] if g == gs[0] else ''} & {TEX_G[g]} & " + " & ".join(f4(c) for c in cells) + " \\\\")
    print("\\midrule" if s != surfaces[-1] else "")
print("direct S0 errors N=128:", {s: round(val(s, gs[0], eps[0], N, 'direct_S0', 'hausdorff_Q'), 4) for s in surfaces})

print("\n% ---- Table 2: slopes vs eps of exact interpolation (N = 128) ----")
slopes = defaultdict(list)
for s in surfaces:
    for g in gs:
        end = [endpoint(s, g, e, N, "hausdorff_Q") for e in eps]
        exn = [val(s, g, e, N, "exact_vertex_mid", "hausdorff_near_P") for e in eps]
        exf = [val(s, g, e, N, "exact_vertex_mid", "hausdorff_away_from_P") for e in eps]
        endf = [endpoint(s, g, e, N, "hausdorff_away_from_P") for e in eps]
        a, b, c, d = slope(eps, end), slope(eps, exn), slope(eps, endf), slope(eps, exf)
        slopes[g].append((a, b, c, d))
        print(f"{TEX_F[s] if g == gs[0] else ''} & {TEX_G[g]} & {f3(a)} & {f3(b)} & {f3(c)} & {f3(d)} & {e2(exf[0])} & {e2(exf[2])} \\\\")
    print("\\midrule" if s != surfaces[-1] else "")

print("\n% ---- Topology at N = 128 (components / euler / dist P) ----")
topo = defaultdict(int)
for s in surfaces:
    for g in gs:
        for e in eps:
            d0 = R[(s, g, e, N, "direct_S0")]
            for m in ("endpoint_minus", "endpoint_plus", "vertex_mid_from_minus", "vertex_mid_from_plus",
                      "sdf_average", "exact_vertex_mid"):
                r = R[(s, g, e, N, m)]
                same = r["components"] == d0["components"] and r["euler"] == d0["euler"]
                topo[(m, "same")] += same
                topo[(m, "all")] += 1
                topo[(m, "holes")] += int(r["hole_edges"]) > 0
                topo[(m, "nonmanifold")] += int(r["nonmanifold_edges"]) > 0
                topo[(m, "flipped")] += r["flipped_triangles"] not in ("", "0")
for m in ("endpoint_minus", "endpoint_plus", "vertex_mid_from_minus", "vertex_mid_from_plus", "sdf_average", "exact_vertex_mid"):
    print(m, {k: topo[(m, k)] for k in ("same", "all", "holes", "nonmanifold", "flipped")})

print("\n% ---- topology per surface/G at N=128, eps=1e-3: (comp, euler) for S0, -, +, mid-, mid+, sdf ----")
for s in surfaces:
    for g in gs:
        e = 1e-3
        t = [(R[(s, g, e, N, m)]["components"], R[(s, g, e, N, m)]["euler"]) for m in
             ("direct_S0", "endpoint_minus", "endpoint_plus", "vertex_mid_from_minus", "vertex_mid_from_plus", "sdf_average")]
        print(s, g, t)

print("\n% ---- resolution dependence (eps = 1e-3), median over all pairs of H_Q ----")
for n in (32, 64, 128):
    def med(xs):
        xs = sorted(xs)
        return xs[len(xs) // 2]
    line = []
    for m in ("direct_S0", "endpoint", "vertex_mid", "sdf_average", "exact_vertex_mid"):
        xs = []
        for s in surfaces:
            for g in gs:
                if m == "endpoint":
                    xs.append(endpoint(s, g, 1e-3, n, "hausdorff_Q"))
                elif m == "vertex_mid":
                    xs.append(vmid(s, g, 1e-3, n, "hausdorff_Q"))
                else:
                    xs.append(val(s, g, 1e-3, n, m, "hausdorff_Q"))
        line.append(med(xs))
    print(n, " & ".join(f"{x:.4f}" for x in line))

print("\n% ---- dist(P, reconstruction), N=128, eps=1e-3 ----")
for s in surfaces:
    print(s, [f"{g}:{float(R[(s, g, 1e-3, N, 'exact_vertex_mid')]['dist_P_to_mesh']):.1e}/{float(R[(s, g, 1e-3, N, 'sdf_average')]['dist_P_to_mesh']):.1e}" for g in gs])

print("field average max error:", max(float(r["field_average_max_error"]) for r in rows))
print("correspondence failures:", sum(int(r["correspondence_failures"]) for r in rows),
      "distance failures:", sum(int(r["distance_failures"]) for r in rows))

# pgfplots data: error vs eps for a few pairs at N = 128
os.makedirs(DATA, exist_ok=True)
for tag, s, g in (("A1_1", "A1_x2+y2-z2+xyz", "1"), ("A1_x", "A1_x2+y2-z2+xyz", "x"),
                  ("E6a_x", "E6_x3+y3+z3", "x"), ("A2a_x", "A2_x2+y2+z3", "x")):
    with open(os.path.join(DATA, f"interp_{tag}.dat"), "w") as f:
        f.write("eps endpoint vertexmid sdf exact exactnear exactfar endpointfar direct hybridnear conenear\n")
        for e in sorted(eps):
            f.write(" ".join(f"{v:.6e}" for v in (
                e, endpoint(s, g, e, N, "hausdorff_Q"), vmid(s, g, e, N, "hausdorff_Q"),
                val(s, g, e, N, "sdf_average", "hausdorff_Q"), val(s, g, e, N, "exact_vertex_mid", "hausdorff_Q"),
                val(s, g, e, N, "exact_vertex_mid", "hausdorff_near_P"),
                val(s, g, e, N, "exact_vertex_mid", "hausdorff_away_from_P"),
                endpoint(s, g, e, N, "hausdorff_away_from_P"), val(s, g, e, N, "direct_S0", "hausdorff_Q"),
                val(s, g, e, N, "hybrid_cone", "hausdorff_near_P"),
                val(s, g, e, N, "singular_cone", "hausdorff_near_P"))) + "\n")
print("wrote", DATA)


# ============================================================================
# Beyond linear interpolation: gradient-line correspondence, hybrid and
# singular cone patch (src/analysis/Reconstruction.h).
# ============================================================================
def med(xs):
    """Upper median, as in the tables above; an empty reconstruction counts as +inf."""
    xs = sorted(x for x in xs if not math.isnan(x))
    return xs[len(xs) // 2] if xs else float("nan")


def g4(x):
    if not math.isfinite(x):
        return "--"
    return f"{x:.4f}" if x >= 0.001 else f"{x:.5f}"


def contractible(r):
    return (r["components"] == "1" and r["euler"] == "1" and r["hole_edges"] == "0"
            and r["nonmanifold_edges"] == "0")


ALL = ["direct_S0", "endpoint", "vertex_mid", "exact_vertex_mid", "normal_mid", "sdf_average", "hybrid_cone",
       "singular_cone"]
LABEL = {"direct_S0": "direct mesh of $F$", "endpoint": "$S_{\\pm\\varepsilon}$ (no interpolation)",
         "vertex_mid": "vertex interpolation (mesh)", "exact_vertex_mid": "vertex interpolation (exact)",
         "normal_mid": "vertex interpolation along $\\nabla F$", "sdf_average": "SDF interpolation",
         "hybrid_cone": "hybrid: interpolation + cone patch", "singular_cone": "singular cone patch"}


def err(s, g, e, n, m, col):
    if m == "endpoint":
        return endpoint(s, g, e, n, col)
    if m == "vertex_mid":
        return vmid(s, g, e, n, col)
    return val(s, g, e, n, m, col)


print("\n% ---- NEW Table: median near / away for all methods, N = 128 ----")
for m in ALL:
    near = [med([err(s, g, e, N, m, "hausdorff_near_P") for s in surfaces for g in gs]) for e in eps]
    far = [med([err(s, g, e, N, m, "hausdorff_away_from_P") for s in surfaces for g in gs]) for e in eps]
    print(f"{LABEL[m]:45s} & " + " & ".join(g4(x) for x in near + far) + " \\\\")

print("\n% ---- NEW: pairs within twice the direct error (N = 128) ----")
for m in ALL[1:]:
    cnt = []
    for e in (1e-2, 1e-3, 1e-4):
        c = 0
        for s in surfaces:
            for g in gs:
                c += err(s, g, e, N, m, "hausdorff_Q") <= 2 * val(s, g, e, N, "direct_S0", "hausdorff_Q")
        cnt.append(c)
    print(f"{LABEL[m]:45s} & " + " & ".join(str(c) for c in cnt) + " \\\\")

print("\n% ---- NEW: topology over the 200 runs at N = 128 ----")
for m in ["direct_S0", "endpoint_minus", "endpoint_plus", "vertex_mid_from_minus", "vertex_mid_from_plus",
          "exact_vertex_mid", "normal_mid", "sdf_average", "hybrid_cone", "singular_cone"]:
    rs = [R[(s, g, e, N, m)] for s in surfaces for g in gs for e in eps]
    ok = sum(contractible(r) for r in rs)
    flips = sum(1 for r in rs if r["flipped_triangles"] not in ("", "0"))
    nm = sum(1 for r in rs if r["nonmanifold_edges"] != "0")
    holes = sum(1 for r in rs if r["hole_edges"] != "0")
    flipcol = "--" if all(r["flipped_triangles"] == "" for r in rs) else str(flips)
    print(f"{m:24s} & {ok} & {flipcol} & {nm} & {holes} \\\\   (of {len(rs)})")

print("\n% ---- NEW: per surface, direct MC vs singular cone (G-independent; G = 1, eps = 1e-2 rows) ----")
for s in surfaces:
    cells = []
    for n in (32, 64, 128):
        d = R[(s, "1", 1e-2, n, "direct_S0")]
        c = R[(s, "1", 1e-2, n, "singular_cone")]
        cells += [g4(float(d["hausdorff_near_P"])), g4(float(c["hausdorff_near_P"]))]
    d = R[(s, "1", 1e-2, 128, "direct_S0")]
    c = R[(s, "1", 1e-2, 128, "singular_cone")]
    topo = f"{d['components']}/{d['euler']} & {c['components']}/{c['euler']}"
    print(f"{TYPE[types[s]]} & {TEX_F[s]} & " + " & ".join(cells) + f" & {topo} & {c['patch_loops']} \\\\")

print("\n% ---- NEW: orders in h of near-P error (64 -> 128) ----")
for m in ("direct_S0", "singular_cone"):
    xs = []
    for s in surfaces:
        a = float(R[(s, "1", 1e-2, 64, m)]["hausdorff_near_P"])
        b = float(R[(s, "1", 1e-2, 128, m)]["hausdorff_near_P"])
        xs.append(math.log2(a / b))
    print(m, "median order", round(med(xs), 2), [round(x, 2) for x in xs])

print("\n% ---- NEW: resolution dependence (eps = 1e-3), median H_Q and near ----")
for n in (32, 64, 128):
    print(n, " & ".join(g4(med([err(s, g, 1e-3, n, m, "hausdorff_Q") for s in surfaces for g in gs])) for m in ALL))
    print(n, "near", " & ".join(g4(med([err(s, g, 1e-3, n, m, "hausdorff_near_P") for s in surfaces for g in gs])) for m in ALL))

print("\n% ---- NEW: normal_mid vs exact_vertex_mid away from P ----")
b = w = t = 0
for s in surfaces:
    for g in gs:
        for e in eps:
            a1 = val(s, g, e, N, "exact_vertex_mid", "hausdorff_away_from_P")
            a2 = val(s, g, e, N, "normal_mid", "hausdorff_away_from_P")
            if a2 < 0.9 * a1: b += 1
            elif a2 > 1.1 * a1: w += 1
            else: t += 1
print("better", b, "worse", w, "same", t)
print("normal_mid correspondence fallbacks (N=128):",
      sum(int(R[(s, g, e, N, "normal_mid")]["correspondence_failures"]) for s in surfaces for g in gs for e in eps))

print("\n% ---- NEW: hybrid near-P slopes vs eps ----")
for s, g in (("E6_x3+y3+z3", "1"), ("E6_x3+y3+z3", "x"), ("A1_x2+y2-z2+xyz", "1"), ("A2_x2+y2+z3", "x")):
    hy = [val(s, g, e, N, "hybrid_cone", "hausdorff_near_P") for e in eps]
    en = [endpoint(s, g, e, N, "hausdorff_near_P") for e in eps]
    print(s, g, "hybrid", [g4(x) for x in hy], "endpoint", [g4(x) for x in en])
print("hybrid radius examples:", {e: R[("E6_x3+y3+z3", "1", e, N, "hybrid_cone")]["patch_radius"] for e in eps})
