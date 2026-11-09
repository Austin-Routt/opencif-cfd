"""Geometry helpers and the per-gap flow of the five-gap slab, shared by Python CFD's post-processor and by the
comparison with COMSOL's exports (the same integration on the same points for both codes)."""
import numpy as np


def gap_edges(G):
    pitch = G['L_post'] + G['L_gap_x']
    return [(G['L_post'] + k * pitch, G['L_post'] + k * pitch + G['L_gap_x']) for k in range(G['n_gaps'])]


def in_gap(X, G):
    m = np.zeros(np.shape(X), bool)
    for x0, x1 in gap_edges(G):
        m |= (X >= x0 - 1e-9) & (X <= x1 + 1e-9)
    return m


def in_solid(X, Y, G):
    """Points strictly inside the wall between the channels (outside the gap openings)."""
    wc2, gap = G['w_c'] / 2, G['gap']
    inside = np.zeros(np.shape(X), bool)
    for x0, x1 in gap_edges(G):
        inside |= (X > x0 + 1e-9) & (X < x1 - 1e-9)
    return (Y > wc2 + 1e-9) & (Y < wc2 + gap - 1e-9) & ~inside


def gap_flows(xs, zs, V, G, depth_factor=1.0):
    """Flow through each gap from v on the gap mid-plane: V[ix, iz] on the tensor grid (xs, zs) in um, m/s;
    trapezoid in z and x over each opening (NaN on the post faces counts as 0, no slip); m^3/s, times depth_factor
    (2 for a half-depth plane)."""
    q = []
    for x0, x1 in gap_edges(G):
        sel = (xs >= x0 - 1e-9) & (xs <= x1 + 1e-9)
        Vk = np.nan_to_num(V[sel])
        q.append(float(np.trapezoid(np.trapezoid(Vk, zs, axis=1), xs[sel]) * 1e-12 * depth_factor))
    return q


GEOM = {'V1': dict(w_c=52., w_s=44., gap=19., L_post=27., L_gap_x=22., depth=140., n_gaps=5),
        'V3': dict(w_c=52., w_s=57., gap=19., L_post=27., L_gap_x=22., depth=140., n_gaps=5)}
