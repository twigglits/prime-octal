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

## Appendix F: the repulsion IS the singular series (resolution)

The council's principled null was the **Hardy–Littlewood singular series** 𝔖, not a random
subset. `--sigma i|w R` computes, for every small offset δ, the measured prime/random pair
ratio and overlays it on 𝔖(δ) computed independently from the Euler product over prime ideals
`𝔖(δ) = ∏_𝔭 (1 − ω_𝔭/N𝔭)/(1 − 1/N𝔭)²`, ω_𝔭 = 1 if 𝔭|δ else 2. Result at R=4000:

| lattice | admissible offsets | corr(obs, 𝔖) | mean(obs/𝔖) |
|---------|--------------------|--------------|-------------|
| ℤ[i]  | 30 of 60 (half forbidden) | **0.9999** | 1.005 |
| ℤ[ω]  | **63 of 63 (none forbidden)** | **1.0000** | 1.006 |

The measured pair ratios reproduce 𝔖 offset-by-offset to ~0.5% (a uniform finite-size scale
factor, constant across offsets — it does not affect the correlation). This **explains every
earlier observation**:

1. **The "repulsion" is exactly 𝔖, not a new phenomenon.** Per offset the ratio is *not*
   uniformly <1 — admissible offsets are *enhanced* (ℤ[i]: 1.68 at norm 2/4/8, 2.24 at norm
   10/20; ℤ[ω]: up to 1.83). The nearest-neighbour "primes farther apart than random" was the
   aggregate of forbidden offsets plus enhanced survivors — all encoded in 𝔖.
2. **The ℤ[i] > ℤ[ω] gap is the parity obstruction, and nothing more.** ℤ[i] has a norm-2
   prime ideal (1+i), so every odd-norm offset is forbidden (𝔖 = 0 — one of the pair is always
   even): *half* of ℤ[i]'s offsets are dead. ℤ[ω] has **no norm-2 ideal**, so **all 63 offsets
   are admissible** — zero hard obstruction. That structural difference (which the gap-statistic
   only saw as "8.6% vs 5.6%") is the whole story.
3. **Confirms "robust ≠ anomalous."** The repulsion is a real, ~500σ measurement *and* fully
   predicted by classical analytic number theory. No anomaly survives.

So the lattice-prime two-point structure in both rings is the Hardy–Littlewood singular series,
verified empirically at GPU scale (corr ≈ 1.0000). The genuinely *novel* slice, per the domain
council member, would be the **Eisenstein** sector-variance / Hecke-angle form factor (the ℤ[ω]
analogue of Kurlberg–Rudnick / Rudnick–Waxman) — empirically unstudied — but that is a new line,
not a loose end.

Reproduce: `./bin/lattice --sigma i|w R [seeds]`.

## Appendix G: Hecke-angle sector variance (exploratory)

The Rudnick–Waxman / Kurlberg–Rudnick object: fold prime angles into the rotational
fundamental domain `[0, 2π/u)` (u = unit order, 4 for ℤ[i], 6 for ℤ[ω]) and measure the
variance of prime counts across K equal angular sectors — the index of dispersion
`Ψ = Var(N_k)/mean(N_k)` (Poisson ⇒ 1). A matched-density random lattice subset is the
baseline. `--formfactor i|w R` sweeps K. At R=6000:

| K | ℤ[i] Ψ_prime | Ψ_random | ratio | ℤ[ω] Ψ_prime | Ψ_random | ratio |
|------|------|------|------|------|------|------|
| 30   | 0.64 | 1.28 | 0.50 | 0.95 | 1.72 | 0.55 |
| 120  | 1.15 | 1.25 | 0.91 | 1.79 | 1.62 | 1.11 |
| 720  | 1.73 | 1.23 | 1.41 | 2.92 | 1.59 | 1.84 |
| 2520 | 2.48 | 1.27 | 1.95 | 3.89 | 1.53 | 2.54 |

**Qualitative finding (robust):** in *both* rings the prime/random ratio rises monotonically
through 1 (crossover near K≈120) — prime angles are **more uniform than a random lattice
subset at coarse angular scales (rigidity)** and **more clustered at fine scales (excess
variance)**. The u=4 and u=6 lattices behave the same way, suggesting the transition is
universal across the two rings.

**Caveats (why this is exploratory, not a clean claim):**
1. *The random baseline is not Poisson.* Ψ_random ≈ 1.25 (ℤ[i]) / 1.6 (ℤ[ω]), roughly
   K-independent — it carries the lattice's intrinsic angular granularity. That is *why* it is
   the right control, but it means "Ψ_prime vs 1" is not the signal; "Ψ_prime vs Ψ_random" is.
2. *ℤ[i] vs ℤ[ω] at equal K is apples-to-oranges* (the same mistake the pcf made). The
   fundamental domains differ (π/2 vs π/3), so equal K means different sector *width*. The
   quantitative cross-ring gap in this table is therefore **not** meaningful; a fair comparison
   needs matched angular width or matched primes-per-sector. The shared *qualitative* transition
   is what survives.
3. The fine-K excess could carry residual discrete-lattice angular structure not fully removed
   by the random control. A definitive form factor would compare against the Hecke L-function
   prediction (matched resolution, thin norm shells) — genuinely future work.

So: a real, reproducible scale-dependent angular form factor with the same shape in both rings,
reported with its caveats. Unlike the singular series (Appendix F, definitive), this one is an
exploratory measurement, not a closed result.

Reproduce: `./bin/lattice --formfactor i|w R`.

## Appendix H: matched-resolution form factor (and the folding artifact)

A council pass (skip→method) replaced Appendix G's flawed comparison with three fixes:
(1) a **thin norm shell** instead of the full disk (full disk averages over scales);
(2) **x-axis = λ = mean primes per sector**, which absorbs the differing fundamental-domain
widths (π/2 vs π/3) so ℤ[i] and ℤ[ω] are directly comparable; (3) the control is a
matched-density Bernoulli subset of **all lattice integers in the shell**, sharing the prime
set's angular "ray comb" so binning aliasing cancels. `--ffmatched i|w R`.

**An artifact caught mid-flight.** The first matched run gave Ψ_prime saturating at ≈3.9 (u=4)
and ≈5.8 (u=6) — suspiciously ≈ u. The cause: folding angles with `fmod` into `[0,2π/u)` piles
all *u* unit-associates of each prime ideal onto the *same* angle, inflating the dispersion
index by exactly u. Fix: **restrict to one fundamental sector** (each ideal counted once)
rather than fold. The selftest now gates on coarse Ψ_prime < 1.5 to catch any regression.

**Result (R=5000, corrected).** In-sector prime counts match to 0.02% (1,134,662 vs 1,134,924),
so the two rings are compared at genuinely equal λ:

| λ (primes/sector) | ℤ[i] Ψ_prime | ℤ[ω] Ψ_prime |
|-------------------|--------------|--------------|
| 2251 | 0.485 | 0.565 |
| 225  | 0.827 | 0.795 |
| 41   | 0.983 | 0.952 |
| 7.9  | 0.982 | 0.973 |
| 1.6  | 0.979 | 0.975 |

Two findings:

1. **Number-variance suppression (rigidity → Poisson).** Ψ_prime is strongly **sub-Poisson**
   (≈0.5) at coarse angular resolution (large λ) and relaxes to Poisson (Ψ→1) as λ→O(1) — the
   Rudnick–Waxman / CUE picture: prime angles are *more ordered than random* at large scales,
   washing out to shot noise at the finest scales.
2. **Universality across the rings.** The ℤ[i] (u=4) and ℤ[ω] (u=6) curves **collapse** onto
   one another at matched λ (within a few %). The Eisenstein angular form factor — empirically
   unstudied — behaves like the Gaussian one. That collapse is the novel result.

**Honest limits.** The shell is a factor-4 norm window (not strictly dyadic); the control is
Bernoulli (random) not a deterministic same-support set, so the small Ψ_prime−Ψ_control
residual is suggestive, not definitive; and no error bars / CUE-curve overlay yet. The robust,
defensible claims are the *shape* (sub-Poisson→Poisson) and the *cross-ring collapse* — not a
quantitative CUE fit. Catching the u-fold folding artifact is itself the main lesson:
fundamental-domain folding silently multiplies the variance by the unit order.

Reproduce: `./bin/lattice --ffmatched i|w R`.

## Appendix I: universality collapse, and an honest CUE verdict

The definitive test of a form factor is not one curve but **Rudnick–Waxman universality**: the
statistic should depend only on the scaling variable `δ = log K / log X` (X = norm scale), so
runs at *different R* and *different lattices* must collapse onto one master curve. `--ffmatched`
emits δ per row; sweeping R ∈ {3000, 6000, 10000} for both rings and pooling all points
(`docs/formfactor_collapse.py`):

**1. The collapse is real.** Ψ_prime(δ) from u=4 and u=6, across all three R, falls on a single
master curve — **residual RMS = 0.043** about a linear ramp `Ψ ≈ 2.03 δ − 0.20` that saturates
at Poisson (Ψ→1) near δ ≈ 0.59. The per-δ spread across lattice and R is ≤ 0.015 except at the
coarsest (few-sector, noisiest) point. So the angular form factor is a function of δ alone,
**identical for the Gaussian and the (previously unstudied) Eisenstein lattice.** That cross-ring,
cross-scale collapse is the real result.

**2. But it is NOT the CUE eigenvalue statistic** — tested and rejected. Converting to number
variance `V = Ψ·λ` and comparing:

| λ | V_data | V_Poisson | V_CUE | data/Poisson | data/CUE |
|------|--------|-----------|-------|--------------|----------|
| 8316 | 4784 | 8316 | 1.26 | 0.57 | 3795 |
| 3173 | 1542 | 3173 | 1.16 | 0.49 | 1326 |
| 863  | 437  | 863  | 1.03 | 0.51 | 424  |

The data is mildly **sub-Poisson (~0.5×)**, but the CUE number variance (the `(1/π²)log` law) is
3–4 *orders of magnitude* smaller. So at the accessible scales (R ≤ 10⁴) the primes are in a
weak-rigidity regime, **not** the CUE-rigid regime. A genuine CUE identification — if it emerges
at all — lives in the asymptotic X→∞ window with the exact RW smooth normalization, which this
finite-R sector-count setup does not reach. Claiming a "CUE fit" here would be wrong; the honest
statement is: **a universal sub-Poisson form factor (same in both rings), quantitatively far from
CUE at these scales.**

The ramp-saturating-at-Poisson *shape* is reminiscent of the U(N) form factor min(τ,1), but
shape resemblance is not a fit — the magnitudes disagree by 10³.

Reproduce: sweep `./bin/lattice --ffmatched i|w R` over several R; `python3 docs/formfactor_collapse.py`.

### Parked

Open directions, none a loose end: the √36 column-sweep moat; and — if the CUE regime is the
goal — reaching the asymptotic RW window (much larger X, the exact smooth-weighted normalization,
narrow dyadic shells), which is a research program, not a finite-R measurement. The two-point
story is closed (the singular series, Appendix F); the angular story is a confirmed universal
sub-Poisson form factor that is demonstrably not CUE at reachable scales.

### Why this is the *right* place to look

Unlike base tricks, the lattice carries real arithmetic structure (Hecke characters, class
field theory) and a real open problem. It still doesn't yield a fast *generator* — primes
remain ~1/ln n dense in every direction — but it is the frontier where a CUDA survey can
measure something not already settled by the wheel.

Reproduce: `make test` (includes `lattice --selftest`); `./bin/lattice --hecke R [bins]`,
`./bin/lattice --moat K R`.
