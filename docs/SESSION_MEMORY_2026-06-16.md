# prime-octal — local session memory (claudianDB fallback, 2026-06-16)

> Written locally because claudianDB (the `claudiandb` MCP server) was temporarily
> unavailable. Mirror these facts into the knowledge graph when it's back up.

Repo: `~/Documents/prime-octal`, branch **`gaussian-lattice-moat`** (pushed to
`github.com/twigglits/prime-octal`). Single RTX 5090, 31 GB host RAM, nvcc 12 (compute_80 JIT).

---

## HEADLINE (latest): the true √26 Gaussian moat, computed for the first time

Commit **`1097c6d`**. The column-sweep engine ran `--moat-sweep 26 1050000` to completion:

- **k²=26 → BOUNDED, farthest |z| = 1,015,638.7651** — matches Tsuchimura's exact published
  value (1,015,638.765) **to the digit**. Full validation of the engine at the √26 scale.
- **Peak memory 246 MB**; peak union-find records **730,549** for a component of **~14.5 billion
  primes** (Tsuchimura octant 14,542,615,005). Reference-counted recycling: live records rose to
  730k by ~36 % of the sweep, then fell to 230k as the component finalized behind the window.
- Wall time **12 h 41 m**, 100 % CPU single-thread. Full run log: `docs/moat_sqrt26_run.log`.
- The full-disk bitmap method would have needed **~322 GB** at R=1.05M; the sweep used 246 MB
  (**~1300× less**). The O(R²) memory wall (which capped the old method at R≈400k) is GONE.

## The column-sweep engine — `src/moat_sweep.cu` (commit `69579a5`, refined `2a13447`,`71a8f7e`)

- CLI `--moat-sweep KSQ R [blockW]` (KSQ = integer k²); `--selftest`. In Makefile + `make test`.
- O(k·R) memory: with steps ≤ k (|dx| ≤ ⌊√k²⌋), a left-to-right column sweep needs only a k+1
  column window. Sweep-line **Hoshen–Kopelman union-find**: cluster roots carry aggregates
  (max-norm, origin-flag, in-window count) so departed cells keep contributing via their cluster
  and never need relabeling. Records reference-counted + recycled. GPU computes primality one
  column-block at a time (Miller–Rabin on the norm). Half-plane reduction (sweep a≥0) via the
  (a,b)→(−a,b) reflection symmetry. Seed pinned at (1,1) as a stable origin handle.
- Bugs found & fixed during build: negative-column segfault (guard a−da≥0); recycled-id origin
  handle (pin the seed record). `find` is path-halving (perf-neutral — bottleneck is cache misses
  across the union-find arrays, not find logic; single-thread ~100 % CPU).
- Validated: small moats (k²=2→11.7047, 4→45.3100, 8→93.4719); **√20 → 133,679.0655** (123 MB,
  13 min); **√26 → 1,015,638.7651** (246 MB, 12h41m). Both match Tsuchimura.

## EARLIER CORRECTION (important): old "√26 = 133679" was actually √20 (commit `15cc052`)

The earlier full-disk run `--moat-gpu 5.099 150000` reported "√26, |z|=133679, 547,583,245
primes" — that was **wrong, it was the √20 moat**. Two causes: (1) float bug, `5.099²=25.9998<26`
excluded the norm-26 steps (5,1),(1,5) → computed the k²≤25 plateau = √20's value; (2) R=150000
≪ the true √26 radius (10⁶). Caught only via EXTERNAL check against Tsuchimura's table (internal
GPU-vs-CPU validation couldn't catch it — both shared the same float K²). FIX: `moat`/`moat_gpu`/
`emoat`/`emoat_gpu` now take integer KSQ (`dx²+dy²≤KSQ` exactly). Verified by re-running:
k²=20 and k²=25 both give 547,583,245 / 133679.0655 (= √20), k²=26 escapes R=150000.
Lesson: external ground-truth verification is essential.

## OPEN frontier: √36

Now **memory-feasible** (the wall is gone) but needs ~5000× the √26 compute (~weeks
single-threaded). Tsuchimura has only the **upper bound |ξ(√36)| < 80,015,782**; exact unknown.
Time-feasibility would need parallelizing the union-find (hard) or a faster connectivity method.

## Reference: Tsuchimura's exact moat table (METR04-13, U. Tokyo)

| k | farthest \|ξ(k)\| | octant component size |
|---|---|---|
| √20 | 133,679.065 | 273,791,623 |
| √26 | 1,015,638.765 | 14,542,615,005 |
| √32 | 2,823,054.542 | 103,711,268,594 |
| √34 | < 24,289,452 (bound) | Finite |
| √36 | < 80,015,782 (bound) | Finite |

---

## The rest of the lattice investigation (UNAFFECTED by the moat bug — read prime bitmaps directly)

These were saved to claudianDB earlier this session (should already be in the graph):

1. **octal vs hex / base wheels** (docs/octal-vs-hex.md, docs/base_compare.py): representation
   buys only a constant-factor prime wheel, nothing asymptotic. Locality = ord_p(b). Octal is the
   smallest base with a local (digit-sum) rule for 7. Base 15/21 strictly dominate the mod-210
   wheel; primorial bases (210) collapse the wheel to one digit; binary wins for raw speed.

2. **Singular series = the strongest result** (commit `c3abdf9`, `--sigma`): the prime "repulsion"
   in Z[i]/Z[w] is EXACTLY the Hardy–Littlewood singular series. corr(obs,𝔖)=**0.9999** (Z[i]),
   **1.0000** (Z[w]) across small offsets. Z[i]: 30/60 offsets admissible (1+i parity obstruction
   kills the odd-norm half); Z[w]: 63/63 admissible (no norm-2 ideal). The gap-statistic "8.6% vs
   5.6% repulsion" is just the two rings' singular series differing — definitional, not anomalous.

3. **Eisenstein Z[w] port** (`77abdbd`): hexagonal lattice percolates further per step; Hecke
   equidistribution flat with 12-fold symmetry. `--ehecke --emoat --emoat-gpu`.

4. **Retracted**: the "Z[i] more angularly anisotropic" pcf claim was a sector-count artifact
   (8-vs-12 bins). The angular form-factor "CUE fit" was rejected — data is mildly sub-Poisson but
   NOT the CUE log-law (off by 10³); universal scaling collapse across rings IS real though.
   A u-fold folding artifact in the first form-factor attempt was caught and fixed.

5. **Overarching honest conclusion**: the project does NOT find a more efficient/generalized way
   to discover primes — and rigorously shows why not (constant-factor wheel ceiling; lattice
   structure = exactly the known Hardy–Littlewood/Hecke theory; no anomaly). Value = reproductions
   + negative results + one genuinely novel empirical bit (Eisenstein form-factor universality) +
   the column-sweep engine that reaches √26 (and makes √36 memory-feasible).

Commit chain (this session): cef3baf 77abdbd 1e0c8c2 c38fff5 93787dc c3abdf9 59db126 cd2b9c4
a7ff1e3 15cc052 b34564d bbf10cd 69579a5 2a13447 71a8f7e 1097c6d.
