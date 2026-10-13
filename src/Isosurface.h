// Isosurface.h
//
// Extracts the zero level set {F(x,y,z;k) = 0} of a user-supplied scalar
// field on a regular grid.
//
// Surface extraction is done with MARCHING TETRAHEDRA rather than the
// classic 256-case Marching Cubes table: each grid cube is split into 6
// tetrahedra sharing its main diagonal, and each tetrahedron is cut
// in one of only 3 non-trivial ways (1, 2, or 3 of its 4 corners
// "inside"). That is simple enough to implement directly and correctly
// (no large hand-transcribed lookup table, no ambiguous-face cases),
// while being the standard, well-known cube-decomposition variant of
// marching cubes used whenever topological robustness matters more than
// raw triangle count.
//
// Every generated vertex gets its position from linear interpolation
// along the cut tetrahedron edge, but its NORMAL / GRADIENT is obtained
// by re-evaluating the field with Dual3 (exact forward-mode AD) at that
// exact point -- so the gradient you see on the surface is the true
// analytic gradient of F, not a finite-difference or interpolated
// approximation.
#pragma once
#include "Expr.h"
#include <vector>
#include <array>
#include <cstdint>
#include <chrono>
#include <algorithm>

struct SurfaceVertex {
    float px, py, pz;      // position
    float nx, ny, nz;      // unit normal = normalized gradient
    float gradMag;         // |grad F| at this point (un-normalized)
};

struct GridSample {
    double f;
    double gx, gy, gz; // analytic gradient, from Dual3
};

struct IsosurfaceResult {
    std::vector<SurfaceVertex> vertices; // non-indexed triangle soup, 3 per tri
    double buildMillis = 0.0;
    double maxGradMag = 1e-9; // for color normalization
};

// The 6-tetrahedra decomposition of a unit cube, sharing the diagonal
// corner 0 -> corner 6, built from the hexagonal ring of the remaining
// cube corners taken in cube-edge order (1-2-3-7-4-5-1). Cube corner
// numbering/positions:
//   0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
//   4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
static const int kTetCorners[6][4] = {
    {0, 6, 1, 2},
    {0, 6, 2, 3},
    {0, 6, 3, 7},
    {0, 6, 7, 4},
    {0, 6, 4, 5},
    {0, 6, 5, 1},
};
static const int kCubeCornerOffset[8][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
    {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},
};

class Isosurface {
public:
    // domainHalfExtent: grid spans [-h,h]^3. res: number of cells per axis.
    static IsosurfaceResult build(const Node* ast, double k, int res, double domainHalfExtent) {
        auto t0 = std::chrono::steady_clock::now();
        IsosurfaceResult result;
        if (res < 2) res = 2;

        const int N = res + 1; // grid points per axis
        const double h = domainHalfExtent;
        const double step = (2.0 * h) / res;

        // Precompute value + analytic gradient at every grid point.
        std::vector<GridSample> grid(size_t(N) * N * N);
        auto idx = [N](int i, int j, int kk) { return (size_t(i) * N + j) * N + kk; };

        Dual3 kd(k, 0, 0, 0);
        for (int i = 0; i < N; i++) {
            double x0 = -h + i * step;
            for (int j = 0; j < N; j++) {
                double y0 = -h + j * step;
                for (int kk = 0; kk < N; kk++) {
                    double z0 = -h + kk * step;
                    Dual3 xd(x0, 1, 0, 0), yd(y0, 0, 1, 0), zd(z0, 0, 0, 1);
                    Dual3 r = evalExpr<Dual3>(ast, xd, yd, zd, kd);
                    GridSample& s = grid[idx(i, j, kk)];
                    s.f = r.v; s.gx = r.dx; s.gy = r.dy; s.gz = r.dz;
                }
            }
        }

        auto worldPos = [&](int i, int j, int kk) {
            return std::array<double,3>{ -h + i * step, -h + j * step, -h + kk * step };
        };

        // A helper that, given two tet corners (world pos + sample), linearly
        // interpolates the zero-crossing position and then re-evaluates the
        // field there (via AD) for an exact normal/gradient.
        auto cutEdge = [&](const std::array<double,3>& pa, const GridSample& sa,
                            const std::array<double,3>& pb, const GridSample& sb) {
            double t = sa.f / (sa.f - sb.f); // sa.f and sb.f have opposite signs
            double px = pa[0] + t * (pb[0] - pa[0]);
            double py = pa[1] + t * (pb[1] - pa[1]);
            double pz = pa[2] + t * (pb[2] - pa[2]);
            Dual3 xd(px, 1, 0, 0), yd(py, 0, 1, 0), zd(pz, 0, 0, 1);
            Dual3 r = evalExpr<Dual3>(ast, xd, yd, zd, kd);
            SurfaceVertex v{};
            v.px = (float)px; v.py = (float)py; v.pz = (float)pz;
            double gm = std::sqrt(r.dx*r.dx + r.dy*r.dy + r.dz*r.dz);
            v.gradMag = (float)gm;
            if (gm > 1e-12) { v.nx = (float)(r.dx/gm); v.ny=(float)(r.dy/gm); v.nz=(float)(r.dz/gm); }
            else { v.nx = 0; v.ny = 0; v.nz = 1; }
            return v;
        };

        auto emitTri = [&](SurfaceVertex a, SurfaceVertex b, SurfaceVertex c) {
            // Orient winding so the face normal agrees with the (averaged)
            // analytic gradient direction (F increases from inside to outside).
            float ex1 = b.px - a.px, ey1 = b.py - a.py, ez1 = b.pz - a.pz;
            float ex2 = c.px - a.px, ey2 = c.py - a.py, ez2 = c.pz - a.pz;
            float fnx = ey1*ez2 - ez1*ey2;
            float fny = ez1*ex2 - ex1*ez2;
            float fnz = ex1*ey2 - ey1*ex2;
            float gnx = a.nx+b.nx+c.nx, gny = a.ny+b.ny+c.ny, gnz = a.nz+b.nz+c.nz;
            float dot = fnx*gnx + fny*gny + fnz*gnz;
            float gm = std::max({a.gradMag, b.gradMag, c.gradMag});
            result.maxGradMag = std::max(result.maxGradMag, (double)gm);
            if (dot < 0) {
                result.vertices.push_back(a);
                result.vertices.push_back(c);
                result.vertices.push_back(b);
            } else {
                result.vertices.push_back(a);
                result.vertices.push_back(b);
                result.vertices.push_back(c);
            }
        };

        for (int i = 0; i < res; i++) {
            for (int j = 0; j < res; j++) {
                for (int kk = 0; kk < res; kk++) {
                    std::array<double,3> cp[8];
                    const GridSample* cs[8];
                    for (int c = 0; c < 8; c++) {
                        int ci = i + kCubeCornerOffset[c][0];
                        int cj = j + kCubeCornerOffset[c][1];
                        int ck = kk + kCubeCornerOffset[c][2];
                        cp[c] = worldPos(ci, cj, ck);
                        cs[c] = &grid[idx(ci, cj, ck)];
                    }
                    for (int t = 0; t < 6; t++) {
                        int ia = kTetCorners[t][0], ib = kTetCorners[t][1];
                        int ic = kTetCorners[t][2], id = kTetCorners[t][3];
                        const double fa = cs[ia]->f, fb = cs[ib]->f, fc = cs[ic]->f, fd = cs[id]->f;
                        int mask = (fa < 0 ? 1:0) | (fb < 0 ? 2:0) | (fc < 0 ? 4:0) | (fd < 0 ? 8:0);
                        if (mask == 0 || mask == 15) continue; // entirely outside / inside

                        double fv[4] = { fa, fb, fc, fd };
                        std::array<double,3> pv[4] = { cp[ia], cp[ib], cp[ic], cp[id] };
                        const GridSample* sv[4] = { cs[ia], cs[ib], cs[ic], cs[id] };

                        int nInside = ((mask&1)?1:0) + ((mask&2)?1:0) + ((mask&4)?1:0) + ((mask&8)?1:0);

                        if (nInside == 1 || nInside == 3) {
                            // 1 triangle: isolate the "odd one out" corner (the
                            // single inside vertex if nInside==1, or the single
                            // outside vertex if nInside==3) and cut its 3 edges
                            // to the other three corners.
                            int lone = -1;
                            bool loneInside = (nInside == 1);
                            for (int c = 0; c < 4; c++) {
                                bool isInside = (fv[c] < 0);
                                if (isInside == loneInside) { lone = c; break; }
                            }
                            int others[3]; int oi = 0;
                            for (int c = 0; c < 4; c++) if (c != lone) others[oi++] = c;
                            SurfaceVertex v0 = cutEdge(pv[lone], *sv[lone], pv[others[0]], *sv[others[0]]);
                            SurfaceVertex v1 = cutEdge(pv[lone], *sv[lone], pv[others[1]], *sv[others[1]]);
                            SurfaceVertex v2 = cutEdge(pv[lone], *sv[lone], pv[others[2]], *sv[others[2]]);
                            emitTri(v0, v1, v2);
                        } else {
                            // nInside == 2: the cut is a planar quadrilateral.
                            int insideC[2], outsideC[2]; int ii=0, oo=0;
                            for (int c = 0; c < 4; c++) {
                                if (fv[c] < 0) insideC[ii++] = c; else outsideC[oo++] = c;
                            }
                            int a0 = insideC[0], a1 = insideC[1];
                            int b0 = outsideC[0], b1 = outsideC[1];
                            // Quad in cyclic order: (a0-b0) (a0-b1) (a1-b1) (a1-b0)
                            SurfaceVertex p_a0b0 = cutEdge(pv[a0], *sv[a0], pv[b0], *sv[b0]);
                            SurfaceVertex p_a0b1 = cutEdge(pv[a0], *sv[a0], pv[b1], *sv[b1]);
                            SurfaceVertex p_a1b1 = cutEdge(pv[a1], *sv[a1], pv[b1], *sv[b1]);
                            SurfaceVertex p_a1b0 = cutEdge(pv[a1], *sv[a1], pv[b0], *sv[b0]);
                            emitTri(p_a0b0, p_a0b1, p_a1b1);
                            emitTri(p_a0b0, p_a1b1, p_a1b0);
                        }
                    }
                }
            }
        }

        auto t1 = std::chrono::steady_clock::now();
        result.buildMillis = std::chrono::duration<double, std::milli>(t1 - t0).count();
        return result;
    }
};
