#!/usr/bin/env python3
"""Matplotlib renderings of the two structures in docs/octal-vs-hex.md:

  1. docs/fig-gaussian-moat.png   the Gaussian moat (2D): nested origin-islands by
     step size, with the prime-free moat ring that bounds the sqrt(8) island.
  2. docs/fig-gaussian-lattice.png the Gaussian prime "crystal" (2D full plane):
     the 8-fold-symmetric lattice the moat walks on.
  3. docs/fig-octal-crystal-3d.png the octal prime lattice (3D point cloud):
     X = n mod 8, Y = (n/8) mod 8, Z = n / 64  -- docs/lattice-3d.md, in 3D.

Needs the venv: `. .venv/bin/activate` (numpy + matplotlib). Run: python tools/lattice_plots.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(__file__), "..", "docs")


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i :: i] = False
    return s  # boolean mask: s[k] True iff k prime


def gaussian_prime_mask(R, isprime):
    """For a,b in [-R,R], True iff a+bi is a Gaussian prime."""
    g = np.arange(-R, R + 1)
    A, B = np.meshgrid(g, g, indexing="ij")
    norm = A * A + B * B
    mask = np.zeros_like(A, dtype=bool)
    # both coords nonzero: prime iff norm is a rational prime
    both = (A != 0) & (B != 0)
    mask[both] = isprime[norm[both]]
    # on an axis: prime iff |coord| is a rational prime == 3 (mod 4)
    axis = (A == 0) ^ (B == 0)
    coord = np.abs(A + B)  # the nonzero one (other is 0)
    ax_ok = axis & (coord <= len(isprime) - 1)
    sel = np.where(ax_ok)
    vals = coord[sel]
    mask[sel] = isprime[vals] & (vals % 4 == 3)
    return mask  # indexable as mask[a+R, b+R]


def origin_component(mask, R, ksq):
    """BFS the Gaussian primes reachable from 1+i with step dx^2+dy^2 <= ksq."""
    k = int(ksq**0.5)
    steps = [
        (dx, dy)
        for dx in range(-k, k + 1)
        for dy in range(-k, k + 1)
        if 0 < dx * dx + dy * dy <= ksq
    ]
    seen = np.zeros_like(mask, dtype=bool)
    start = (1 + R, 1 + R)  # 1+i
    if not mask[start]:
        return seen
    stack = [start]
    seen[start] = True
    H = mask.shape[0]
    while stack:
        ax, ay = stack.pop()
        for dx, dy in steps:
            nx, ny = ax + dx, ay + dy
            if 0 <= nx < H and 0 <= ny < H and mask[nx, ny] and not seen[nx, ny]:
                seen[nx, ny] = True
                stack.append((nx, ny))
    return seen


# --------------------------------------------------------------------------- #
# 1. Gaussian moat                                                            #
# --------------------------------------------------------------------------- #
def plot_moat():
    R = 110
    isprime = primes_upto(2 * R * R + 1)
    mask = gaussian_prime_mask(R, isprime)
    # the disk
    g = np.arange(-R, R + 1)
    A, B = np.meshgrid(g, g, indexing="ij")
    indisk = (A * A + B * B) <= R * R
    primes = mask & indisk

    # nested origin-components by step size
    layers = [(2, "√2", "#1b9e77"), (4, "2", "#d95f02"), (8, "√8", "#7570b3")]
    comp = {ksq: origin_component(mask, R, ksq) & indisk for ksq, _, _ in layers}

    fig, ax = plt.subplots(figsize=(8, 8), dpi=140)
    # all Gaussian primes (the moat is the grey ring beyond the islands)
    ax.scatter(A[primes], B[primes], s=3, c="#cfd3da", linewidths=0,
               label="Gaussian primes (the moat is the empty ring)")
    # color each prime by the *smallest* step size that reaches it from 1+i
    assigned = np.zeros_like(mask, dtype=bool)
    for ksq, lbl, col in layers:
        new = comp[ksq] & ~assigned
        assigned |= comp[ksq]
        farthest = np.sqrt((A[new] ** 2 + B[new] ** 2).max()) if new.any() else 0
        ax.scatter(A[new], B[new], s=7, c=col, linewidths=0,
                   label=f"reachable with step ≤ {lbl}  ({new.sum():,} pts, |z|≤{farthest:.0f})")
    ax.add_patch(plt.Circle((0, 0), R, fill=False, ec="#aaaaaa", ls=":", lw=0.8))
    ax.plot(1, 1, "k*", ms=11, label="start 1+i")
    ax.set_aspect("equal")
    ax.set_xlim(-R, R); ax.set_ylim(-R, R)
    ax.set_title("The Gaussian moat — the origin island is sealed off by a prime-free ring\n"
                 "(√8 island: 2,996 primes, farthest |z|≈93.5 — the walk cannot escape)")
    ax.set_xlabel("Re(z) = a"); ax.set_ylabel("Im(z) = b")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    p = os.path.join(OUT, "fig-gaussian-moat.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


# --------------------------------------------------------------------------- #
# 2. Gaussian prime lattice (the crystal the moat lives on)                   #
# --------------------------------------------------------------------------- #
def plot_gaussian_lattice():
    R = 80
    isprime = primes_upto(2 * R * R + 1)
    mask = gaussian_prime_mask(R, isprime)
    g = np.arange(-R, R + 1)
    A, B = np.meshgrid(g, g, indexing="ij")
    indisk = (A * A + B * B) <= R * R
    primes = mask & indisk
    norm = A * A + B * B

    fig, ax = plt.subplots(figsize=(8, 8), dpi=140)
    sc = ax.scatter(A[primes], B[primes], c=np.sqrt(norm[primes]), s=10,
                    cmap="twilight", linewidths=0)
    ax.set_aspect("equal")
    ax.set_xlim(-R, R); ax.set_ylim(-R, R)
    ax.set_title("Gaussian prime lattice ℤ[i] — the crystal the moat walks on\n"
                 "(8-fold dihedral symmetry: a→±a, b→±b, a↔b)")
    ax.set_xlabel("Re(z) = a"); ax.set_ylabel("Im(z) = b")
    fig.colorbar(sc, ax=ax, label="|z| = √(a²+b²)", shrink=0.8)
    fig.tight_layout()
    p = os.path.join(OUT, "fig-gaussian-lattice.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


# --------------------------------------------------------------------------- #
# 3. Octal prime crystal in 3D (docs/lattice-3d.md)                           #
# --------------------------------------------------------------------------- #
def plot_octal_crystal_3d():
    N = 8 ** 4  # 4096 -> z runs 0..63
    isprime = primes_upto(N)
    n = np.nonzero(isprime)[0]
    x = n % 8
    y = (n // 8) % 8
    z = n // 64
    corridor = n % 8  # colour by octal corridor (1,3,5,7)

    fig = plt.figure(figsize=(8, 9), dpi=140)
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(x, y, z, c=corridor, cmap="turbo", s=14, depthshade=True)
    ax.set_xlabel("n mod 8  (units digit)")
    ax.set_ylabel("(n/8) mod 8  (next digit)")
    ax.set_zlabel("n / 64  (order of magnitude)")
    ax.set_title("Octal prime lattice (3D) — primes only land on corridors 1,3,5,7\n"
                 "the empty planes are the forbidden (even / div-by-2) wedges")
    ax.view_init(elev=18, azim=-60)
    fig.colorbar(sc, ax=ax, label="octal corridor (n mod 8)", shrink=0.5, pad=0.1)
    fig.tight_layout()
    p = os.path.join(OUT, "fig-octal-crystal-3d.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    plot_moat()
    plot_gaussian_lattice()
    plot_octal_crystal_3d()
