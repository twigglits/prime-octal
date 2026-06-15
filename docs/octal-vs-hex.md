# Octal vs hexadecimal: can a base "predict" primes?

Short answer: **no base predicts primes** — primality is base-invariant. A base only
*helps* by how cheaply its digit rules detect small prime factors. This note makes the
"relative location of primes in octal vs hex" question precise, answers it, and finds the
base that is *strictly better* than octal for the same job.

## The right measure: locality = multiplicative order

For a base `b` and a prime `p` not dividing `b`, a divisibility rule for `p` always
exists. Its **span** — how far-apart digit positions must interact — is exactly

    ord_p(b)   (the multiplicative order of b mod p)

- `p | b`        → rule reads **trailing digits** (span 0, cheapest)
- `ord_p(b) = 1` → `b ≡ 1 (mod p)` → pure **digit sum** (e.g. casting out nines)
- `ord_p(b) = 2` → `b ≡ -1 (mod p)` → **alternating digit sum**
- larger order   → a longer weighted-window rule (the README's mod-5 octal trick, span 4)

So "relative location of primes in base b" = the order spectrum of b over the small primes.
A small order means the prime's multiples sit at a short, repeating stride in base b —
a *local* rule.

## Octal vs hex, head to head

Free (span ≤ 2) digit-rule primes and the resulting candidate wheel, measured at N = 2×10⁷:

| base       | cheap-rule primes        | wheel       | keeps  | concentration      |
|------------|--------------------------|-------------|--------|--------------------|
| 8 (octal)  | 2, 3, 7 (+5 at span 4)   | 42 (210 w/5)| 28.6%  | 3.50× (4.375× w/5) |
| 16 (hex)   | 2, 3, 5, 17              | 510         | 25.1%  | 3.98×              |

Locality spectrum:

| p  | ord₈ | ord₁₆ | more local |
|----|------|-------|------------|
| 3  | 2    | **1** | hex        |
| 5  | 4    | **1** | hex        |
| 7  | **1**| 3     | **octal**  |
| 17 | 8    | **2** | hex        |

### Why this project lives in base 8

**7 is local in octal and only octal.** `8 ≡ 1 (mod 7)`, so the octal digit sum tests
divisibility by 7 directly. Base 8 is the *smallest* base whose plain digit sum captures
the prime 7. That single fact is what gives the tight mod-210 predictor (keep 48/210 =
22.86%, 4.375× concentration).

Hex's free rules win on 3, 5, 17 — but 17 is a poor trade (you'd rather capture 7), and it
gets 3, 5 only marginally more locally than octal already manages. **Hex is not an upgrade
over octal; it is a different, slightly worse wheel.**

## Pushing further: the actually-optimal base for the mod-210 wheel

If the goal is the {2, 3, 5, 7} wheel, octal is *not* optimal — its rule for 5 costs span 4.
Searching all bases for the cheapest {2,3,5,7} coverage (max span over the four primes):

| base | ord(2,3,5,7) | max span | note |
|------|--------------|----------|------|
| **15** | [1, 0, 0, 1] | **1** | 3,5 trailing (15=3·5); 2,7 digit-sum (15−1=14=2·7) |
| **21** | [1, 0, 1, 0] | **1** | 3,7 trailing (21=3·7); 2,5 digit-sum (21−1=20=2²·5) |
| 6    | [0, 0, 1, 2] | 2 | |
| 14   | [0, 2, 2, 0] | 2 | |
| 8 (octal) | [0, 2, 4, 1] | 4 | the span-4 rule for 5 is the weak link |

**Base 15 and base 21 capture all of 2, 3, 5, 7 with span ≤ 1** — every one of the wheel
primes is testable by a trailing digit or a single digit sum, no weighted window needed.
The reason base 15 is so clean: `15 = 3·5` and `15 − 1 = 14 = 2·7`, so the factorizations of
`b` and `b−1` between them cover exactly {2,3,5,7}.

Going one prime further (the mod-2310 wheel, adding 11): **base 21, 55, and 34** capture the
first five primes {2,3,5,7,11} with local rules — e.g. `55 = 5·11`, `55−1 = 54 = 2·3³`,
`55+1 = 56 = 2³·7`.

### The general law

A base `b` captures a set of small primes maximally locally when those primes divide
`b`, `b−1`, or `b+1`. The best "prime-sieve bases" are therefore those where `b(b−1)(b+1)` —
three consecutive integers — is rich in small prime factors. This is why **primorial-adjacent
bases** (15 = 3·5, 21 = 3·7, 30 = 2·3·5, 55 = 5·11) outperform powers of two for cheap
primality candidate-generation.

## Takeaways

1. No base predicts primes; bases differ only in the *cost* (span = ord_p(b)) of their
   divisibility rules.
2. Octal's distinction is that it is the smallest base with a local (digit-sum) rule for 7.
3. **Base 15 and base 21 strictly dominate octal** for the mod-210 wheel — all four wheel
   primes become span ≤ 1 rules.
4. The optimal prime-sieve bases are those where `b(b−1)(b+1)` packs the most small primes —
   primorial-adjacent bases, not powers of two.

Reproduce: `python3 docs/base_compare.py` (wheels + locality + optimal-base search).

---

# Appendix: primes as a lattice, not a number line (`src/lattice.cu`)

The base/digit-rule story closes a *negative* result: representation buys only a constant
factor. The genuinely open structure in primes is **geometric**, not digit-based — it lives
in the Gaussian integers ℤ[i], where a rational prime *p* splits (`p ≡ 1 mod 4`, `p = x²+y²`),
stays inert (`p ≡ 3 mod 4`), or ramifies (`p = 2`). `src/lattice.cu` probes that geometry.

A Gaussian prime `a+bi` is: a rational prime `≡ 3 mod 4` on an axis, **or** any point whose
norm `a²+b²` is a rational prime. Built on the existing `is_prime_u64` (Miller–Rabin).

### Probe 1 — Hecke angular equidistribution (proven; pipeline sanity gate)

`./bin/lattice --hecke 4000 18` bins the angles `atan2(b,a)` of all Gaussian primes in the
first-quadrant disk. Hecke (1920) proved these **equidistribute** on (0, π/2). Measured:
flat to **reduced χ² = 0.12** at R = 4000 (1,030,827 primes), with the exact `a↔b`
conjugation mirror symmetry. This validates the lattice pipeline the way `π(8^k)` validates
the octal one.

### Probe 2 — the Gaussian moat problem (OPEN)

> Can you walk from the origin to infinity stepping only on Gaussian primes with bounded
> step size *K*?

Conjectured **no** (a prime-free "moat" always eventually blocks you), but unproven.
`./bin/lattice --moat K R` does a BFS over the prime lattice from `1+i`. Measured component
of the origin (BFS exhausted → genuinely bounded, a real moat):

| step K | √ form | component primes | farthest \|z\| |
|--------|--------|------------------|----------------|
| 1.5    | √2     | 100              | 11.70          |
| 2.0    | 2      | 720              | 45.31          |
| 2.9    | √8     | 2 996            | 93.47          |
| 3.2    | √10    | 249 508          | 1 024.35       |

Every tested step size leaves the origin in a finite island — the moat conjecture holds
in range. The √10 jump shows how fast the reachable island grows just before the next moat.

### Probe 2b — reaching the √26 record (GPU disk-sieve)

The Miller–Rabin BFS above dies at √26 (the component is ~10⁹ primes; it times out by
R = 20 000). `--moat-gpu K R` precomputes the Gaussian-prime bitmap on the GPU
(`mark_gaussian` kernel over a rational odd-sieve to R²) so the BFS is pure O(1) lookups,
folded into one quadrant by the 8-fold lattice symmetry. Result:

| step K | √ form | quadrant component | farthest \|z\| | wall time | verdict |
|--------|--------|--------------------|----------------|-----------|---------|
| 5.099  | **√26** | 547 583 245       | **133 679.07** | 4m38s     | **BOUNDED** |

This reproduces the **Gethner–Wagon–Wick √26-moat** (1998) — historically the largest known
moat — on a single RTX 5090: ~2.2 billion Gaussian primes (full plane), origin island
sealed off by a prime-free ring at radius ≈ 133 679. Validation: `make test` asserts the
GPU bitmap's farthest distance matches the CPU Miller–Rabin BFS bit-for-bit for K = 1.5, 2,
√8 before any large run.

### The scaling wall (and the architecture past it)

The next moat is **√36**, which needs ~5000× the √26 compute. The full-disk bitmap can't
get there: its rational-sieve (R² bits) and Gaussian bitmap ((R+1)² bits) blow past the
32 GB GPU at R ≈ 400 000, while √36 lives near R ~ 10⁶ (125 GB bitmap). The fix is to drop
the dense O(R²) bitmaps for a **column sweep**: because steps are bounded by K, BFS
reachability only ever connects cells within K columns, so you retain a moving window of
~K columns of visited state (O(K·R) memory, < 1 GB even at R = 10⁶) and sieve each column's
norms on the fly. That plus distributed compute (the component itself is ~10¹⁰ nodes) is
what the √36 frontier requires — a genuine HPC job, not a single-node afternoon.

Reproduce the record: `./bin/lattice --moat-gpu 5.099 150000` (needs ~7 GB RAM, ~5 min).

---

# Appendix B: the Eisenstein integers ℤ[ω] — a hexagonal cross-check

Following the council's steer (skip the √36 single-workstation grind; the higher-value
move is the less-computed *other* lattice), `src/lattice.cu` also surveys the Eisenstein
integers ℤ[ω], ω = e^{2πi/3}. A point `a+bω` embeds at `z = (a−b/2, b√3/2)`, so its squared
Euclidean norm equals the algebraic norm `N(a,b) = a²−ab+b²` — distances stay exact
integers, and the unit-1 neighbours are the 6 vectors of the hexagonal grid.

**Classification.** `a+bω` is an Eisenstein prime iff `N` is a rational prime (split
`p ≡ 1 mod 3`, or the ramified 3), or it is an associate of an inert prime `p ≡ 2 mod 3`
(sitting on a lattice axis `b=0`, `a=0`, or `a=b`). Verified in `--selftest` against
hand-computed truth (2 inert, 3 ramified, 7 split, etc.).

**Hecke equidistribution.** `--ehecke 4000` gives a flat angle histogram (reduced
χ² = 0.034 over 6.19M primes), exhibiting the lattice's **12-fold dihedral symmetry**
(the 12 bins pair up exactly). The proven-theorem sanity gate holds in ℤ[ω] too.

**Moat comparison (the novel bit).** Same BFS, hexagonal step set, GPU bitmap validated
bit-for-bit against the CPU walk. The hexagonal lattice **percolates markedly further per
step** than the square one — its denser neighbourhood (6 unit-neighbours vs 4) bridges
prime-free rings that stop the Gaussian walk:

| step K | Gaussian ℤ[i] farthest \|z\| | Eisenstein ℤ[ω] farthest \|z\| |
|--------|------------------------------|--------------------------------|
| √2     | 11.70                        | 4.36                           |
| 2      | 45.31                        | **87.85**                      |
| √8     | 93.47                        | 87.85                          |
| 3      | 93.47                        | **2 252.53** (1.58 M primes)   |

The Eisenstein component's K=2→K=3 explosion (88 → 2 252) mirrors the Gaussian √8→√10 jump
but is far more violent — a concrete, reproducible structural difference between the two
rings rather than another record radius. Both moats remain bounded at every tested K,
consistent with the conjecture (no infinite prime walk) in both lattices.

Reproduce: `./bin/lattice --ehecke R [bins]`, `./bin/lattice --emoat K R`,
`./bin/lattice --emoat-gpu K R`.

# Appendix C: nearest-neighbour gap statistics (ℤ[i] vs ℤ[ω])

`--gaps i|w R` finds every prime's nearest prime neighbour in a full-plane disk and
compares the observed mean NN distance to the random (Poisson) model, mean = 1/(2√λ) at
density λ. Both lattices show the primes sitting **farther apart than random** — a mild
repulsion/regularity — and, strikingly, by nearly the same factor at every scale:

| R | ℤ[i] obs/Poisson | ℤ[ω] obs/Poisson |
|------|------|------|
| 1 000 | 1.181 | 1.175 |
| 4 000 | 1.157 | 1.150 |
| 16 000 | 1.139 | 1.131 |

Two clean, reproducible facts:

1. **The excess drifts toward 1 as radius grows** (density thins as ~1/ln r): the prime
   gas looks progressively more Poisson the sparser it gets.
2. **The two lattices coincide to <1%** at every radius — the repulsion is essentially
   *lattice-independent*, a property of the prime distribution, not the hexagonal-vs-square
   geometry. (Gaussian runs ~0.6–0.8% above Eisenstein throughout.)

The NN-distance distributions are discrete spikes at √(norm-form values) — √2, √5, √8… for
ℤ[i]; 1, √3, 2, √7… for ℤ[ω] — the geometric fingerprint of each lattice.

Reproduce: `./bin/lattice --gaps i 8000`, `./bin/lattice --gaps w 8000`.

## Appendix D: the random-subset control (and what it overturned)

The Appendix-C comparison was against *continuous* Poisson, which is the wrong null on a
lattice: no two grid points sit closer than the minimum spacing, so *any* sparse lattice
subset beats continuous-Poisson. The clean control is a **random Bernoulli subset of the
same lattice at the same density** (deterministic splitmix64 hash per cell) — it shares the
primes' hard-core floor exactly. `--gaps` now reports both. The naive ~1.14 excess splits
cleanly in two:

| R | lattice | prime/random *(real signal)* | random/Poisson *(hard-core)* | prime/Poisson *(naive)* |
|------|------|------|------|------|
| 1 000 | ℤ[i] | 1.100 | 1.074 | 1.181 |
| 16 000 | ℤ[i] | **1.082** | 1.052 | 1.139 |
| 1 000 | ℤ[ω] | 1.062 | 1.107 | 1.175 |
| 16 000 | ℤ[ω] | **1.054** | 1.073 | 1.131 |

Two results, both of which **correct Appendix C**:

1. **Roughly half the excess was a lattice artifact.** The hard-core (random/Poisson) carries
   ~5–8% on its own; only the residual prime/random ≈ 5–9% is genuine prime correlation.
2. **The controlled signal is nearly scale-invariant.** prime/random barely moves with radius
   (ℤ[i] 1.10→1.08, ℤ[ω] 1.062→1.054), whereas the naive prime/Poisson drifted 1.18→1.14 —
   the radius-drift was mostly hard-core washing out, *not* a property of the primes.
3. **The lattices are NOT equivalent — Appendix C's "lattice-independent to <1%" was an
   artifact of hard-core contamination.** Properly controlled, **Gaussian primes are
   distinctly more repulsive/regular (~9%) than Eisenstein primes (~5.5%)** — a real,
   reproducible structural difference between the two rings.

This is exactly why the control matters: comparing to continuous Poisson made the two lattices
look identical; comparing to a matched random lattice subset reveals they are not.

**Error bars (16 random seeds, R=8000).** The point sets are huge (>10⁷), so the seed-to-seed
spread is tiny and the repulsion is overwhelmingly significant:

| lattice | prime/random | excess | significance |
|---------|--------------|--------|--------------|
| ℤ[i] | 1.0859 ± 0.0002 | **8.59 %** | ~520 σ |
| ℤ[ω] | 1.0557 ± 0.0002 | **5.57 %** | ~330 σ |

The ~3-point ℤ[i] vs ℤ[ω] gap dwarfs the error bars — it is **not** sampling noise.

**But robust ≠ anomalous.** The correct theoretical null is not a random subset at all but the
**Hardy–Littlewood singular series 𝔖** for each ring (an Euler product over prime ideals, with
ring-specific admissibility). Because 𝔖 differs between ℤ[i] and ℤ[ω] *by construction*, the
measured prime/random ratio is essentially an empirical estimate of 𝔖 — so the 8.6 % vs 5.6 %
gap is most likely the two rings' singular series differing, i.e. *definitional, not a new
phenomenon*. Confirming that needs the 𝔖 overlay (`data/𝔖`, expected ≈ 1); until then the
honest statement is: **a real, high-significance difference whose arithmetic explanation is the
expected singular series, not a surprise.**

Reproduce: `./bin/lattice --gaps i|w R [seeds]` (prints prime, random-control with seed spread).

## Appendix E: angular pair-correlation

`--pcf i|w R` extends the NN moment to the full short-range two-point function: pair counts
at separation `d` for primes vs the matched random subset, resolved by separation angle.

- **Short-range depletion confirmed.** `g_prime/g_random (d<3) ≈ 0.845` in *both* lattices —
  ~15% fewer close prime pairs than random-at-same-density, consistent with the NN result.

**RETRACTED: the "ℤ[i] is more angularly anisotropic" claim was a binning artifact.** An
earlier draft reported an angular-ratio CV of "ℤ[i] = 0.091 vs ℤ[ω] = 0.051" — but those used
*different sector counts* (8 vs 12; the code had `AB = eis?12:8`). With **equal** sector counts
the CV is dominated by how sector boundaries alias against the discrete short-range lattice
vectors, and it neither converges nor preserves the lattice ordering:

| sectors | ℤ[i] CV | ℤ[ω] CV | "more anisotropic" |
|---------|---------|---------|--------------------|
| 6  | 0.015 | 0.035 | ℤ[ω] |
| 8  | 0.091 | 0.025 | ℤ[i] |
| 12 | 0.031 | 0.068 | ℤ[ω] |
| 24 | 0.184 | 0.088 | ℤ[i] |
| 36 | 0.405 | 0.169 | ℤ[i] |

The sign of the comparison **flips with the (arbitrary) sector count** — there is no robust
angular-anisotropy difference here; the CV measures sector-vs-lattice alignment, not physics.
`--pcf` now sweeps the sector count so this is visible rather than hidden. A shell-matched
normalization (prime/random at each exact admissible distance, i.e. against the singular series)
is the only sound way to ask the angular question, and is left as future work.

Reproduce: `./bin/lattice --pcf i|w R [sectors]` (prints g(d<3) and the CV-vs-sector sweep).

### Parked

The √36 column-sweep architecture stays shelved. The principled next step (per the council) is
to compute the Hardy–Littlewood singular series 𝔖 per ring and report the residual `data/𝔖`
for both the gap and pair-correlation statistics — that turns "primes repel" into "primes obey
𝔖, here is the residual." No record-radius grind.

### Why this is the *right* place to look

Unlike base tricks, the lattice carries real arithmetic structure (Hecke characters, class
field theory) and a real open problem. It still doesn't yield a fast *generator* — primes
remain ~1/ln n dense in every direction — but it is the frontier where a CUDA survey can
measure something not already settled by the wheel.

Reproduce: `make test` (includes `lattice --selftest`); `./bin/lattice --hecke R [bins]`,
`./bin/lattice --moat K R`.
