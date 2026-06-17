# prime-octal

Experimental mathematics: hunting for patterns in **prime number emergence through the
octal (base-8) number system**, computed end-to-end in **CUDA**.

Every primality-related computation here is phrased strictly in terms of base-8 digits.
The GPU sieves and analyzes all integers below 8^K (up to 8^12 ≈ 6.9 × 10^10) and bins
every statistic by *octal band* — the number of octal digits — so "emergence" means
"per order of magnitude in base 8".

## Octal wheel vs hex wheel — the comparative experiment

**[`docs/octal-vs-hex.md`](docs/octal-vs-hex.md)** answers the project's headline
question directly: *what is the geometric pattern of prime locations in octal, how does
it compare to hexadecimal, and what does the delta show?* It ships with a **GPU-free CPU
companion** ([`src/cpu_survey.cpp`](src/cpu_survey.cpp)) so the comparison runs anywhere.

The pattern is a **modular wheel**: a number sits on spoke `n mod b`, and a prime can
only land on a spoke coprime to the base — the "prime corridors." Octal shows 4 corridors
(spokes 1,3,5,7); hex shows 8. The **delta is pure resolution**: because 16 = 2·8, every
octal corridor splits 50/50 into two hex corridors on one extra bit — *the hex wheel is
the octal wheel at 2× zoom, with no new prime information.* Both stick at φ(b)/b = 50%
admissible density, because **a power-of-2 base only ever screens the prime 2.** The real
geometric leverage is in **primorial** bases (6, 30, 210, …) — which is exactly what the
octal *multi-digit* predictor below reconstructs (base 210). Figures:
[rays](docs/fig-rays.svg) · [dyadic refinement](docs/fig-refinement.svg) ·
[Chebyshev tilt](docs/fig-chebyshev.svg). Verified to N = 10¹⁰ (π(N) = 455,052,511 ✓).

**[`docs/lattice-3d.md`](docs/lattice-3d.md)** lifts the wheel into 3D: the digit on one
axis, **order of magnitude up the z-axis**, primes as a point cloud
([octal](docs/cloud-octal.ply) / [hex](docs/cloud-hex.ply) `.ply`, spin them in
MeshLab/Blender). The finding: the *composites* form a periodic crystal of divisibility
diagonals (verified period 210); **primes are the aperiodic vacancies in it**. Octal→hex
doubles every plane and lattice slope (`slope_hex = 2·slope_oct mod p`) — same ×2 delta;
a primorial width (30) snaps the crystal into a clean tile while powers of two can only
shear it. Built with `make lattice`.

**[`docs/the-crystal-formula.md`](docs/the-crystal-formula.md)** answers the natural
follow-up: *if the composite crystal is exactly predictable, write the equation and the
holes are the primes.* It does — `isPrime(n) = [gcd(n, ∏_{p≤√n} p) = 1]`, derived and
verified — but proves the honest catch: a **fixed** crystal's holes are primes **+**
impostors (smallest `11²=121`, precision → `4.375/ln n`); the holes equal the primes only
when the crystal **grows to √n**, which is the Sieve of Eratosthenes itself. No fixed
crystal can do it (**primality is not periodic** — proof included). Built with
`make crystal`.

## The octal machinery

Base 8 admits exact divisibility identities, analogous to decimal digit-sum tricks:

| identity | consequence |
|---|---|
| 8 ≡ 1 (mod 7) | n ≡ *octal digit sum* (mod 7) |
| 8 ≡ −1 (mod 3) | n ≡ *alternating octal digit sum* (mod 3) |
| 8^i mod 5 cycles 1,3,4,2 | weighted octal digit sum ≡ n (mod 5) |
| last octal digit | carries n mod 8 (even digit ⟺ even number) |

Stacking the four rules gives the **octal predictor**: a number survives iff its octal
digits prove it indivisible by 2, 3, 5 and 7. Survival is equivalent to gcd(n, 210) = 1
— but computed *only* by digit manipulation in base 8. The predictor keeps 48/210 =
22.857% of integers and misses no prime above 7 (the rule primes 2, 3, 5, 7 are each
rejected by the very rule they generate).

## Findings so far (RTX 5090, N = 8^12)

Full reports: [`docs/sample-report-8^10.txt`](docs/sample-report-8^10.txt),
[`docs/sample-report-8^12.txt`](docs/sample-report-8^12.txt).

1. **Prediction quality is exactly PNT-shaped.** Per band, the fraction of octal-rule
   survivors that are prime tracks (210/48)/ln n to four decimal places by band 12
   (measured 18.050% vs theoretical 18.050%), with 100% recall. The octal digits buy a
   clean 4.375× concentration of primes — and provably nothing more as n grows.
2. **Chebyshev bias is visible in the last octal digit.** Primes split across last
   octal digits 1, 3, 5, 7 (≡ p mod 8). Digit 1 is special: *every odd square ends in
   octal 1*, making it the quadratic-residue class. As theory predicts, ...1 trails its
   three peers in every band ≥ 4 — a cumulative deficit of ≈ 9.0 × 10³ primes at 8^12
   (718,592,908 vs an average 718,601,869 for the others).
3. **No octal palindromic prime has an even digit count.** Every even-length base-8
   palindrome is divisible by 9 (= 8+1). The survey's even bands are all exactly zero;
   odd-band counts run 4, 13, 47, 311, 1944, 12811 (bands 1, 3, 5, 7, 9, 11) —
   15,130 octal palindromic primes below 8^12. Smallest multi-digit: 73 = 0o111.
4. **Leading octal digit of primes is *not* Benford.** It is nearly flat and flattening
   (band 12: 14.89% → 13.91% across digits 1..7, vs Benford-8's 33.3% → 6.4%), exactly
   what a 1/ln n density implies.
5. **Prime gaps avoid 0 mod 8.** Even gaps land on residues 2, 4, 6 (mod 8) ≈ 26–27%
   each, but on 0 (mod 8) only ≈ 19.5% — the octal shadow of the gap-divisible-by-6
   preference. Most common gap: 6, then 12, then 18 (at 8^12 depth).
6. **External validation.** π(8^k) matches published values for *all* k ≤ 12, and the
   maximal gaps found match the known records (282 after 436,273,009 below 8^10; 464
   after 42,652,618,343 below 8^12).
7. **Octal digit frequency in primes is skewed odd.** Aggregate digit shares at 8^12:
   odd digits ≈ 13.7% each, even digits ≈ 11.6%, zero ≈ 10.4% — the forced-odd last
   digit and never-zero leading digit propagate through the representation.

## Build and run

Requirements: an NVIDIA GPU, `nvcc` (CUDA 12.x), GNU make, C++17 host compiler.

```sh
make            # builds bin/prime_octal and bin/test_prime_octal
make test       # full test suite (~19M checks, GPU vs CPU reference)
make run        # default survey: all n < 8^10, report + CSVs in results/
```

Direct use:

```sh
./bin/prime_octal --octal-digits 12 --out results/K12   # 68.7B numbers, ~5 s
./bin/prime_octal --judge 0o3777                        # predict one octal numeral
./bin/prime_octal --no-postpass ...                     # skip CPU gap/palindrome pass
```

K = 12 needs ≈ 4 GiB of GPU memory and ≈ 4 GiB of host RAM for the bitmap.

**No GPU?** The CPU companion needs only a C++17 host compiler and Python 3:

```sh
make cpu                       # builds bin/cpu_survey
make cpu-run N=10000000000     # octal-vs-hex survey to 1e10 (~25 s, report + CSVs)
make figures                   # renders docs/fig-*.svg from results/digits.csv
```

> **Blackwell + older nvcc note:** the Makefile embeds `compute_80` PTX only and lets
> the driver JIT it for newer GPUs (e.g. an RTX 5090 with nvcc 12.0). First run after a
> build pays a few seconds of JIT; the driver caches it.

### `--judge` example

```
$ ./bin/prime_octal --judge 0o3777
octal   3777 / decimal 2047
  vs 2: last digit 7            -> odd, passes
  vs 3: alternating digit sum +4 ≡ 1 (mod 3)  -> passes
  vs 7: digit sum 24 ≡ 3 (mod 7)              -> passes
  vs 5: weighted sum (w=1,3,4,2) ≡ 2 (mod 5)  -> passes
prediction from octal digits: CANDIDATE
ground truth (Miller-Rabin):  not prime        (2047 = 23 × 89)
```

## How it computes

- **Fused GPU kernel** (`src/sieve.cuh`): each block owns a 262,144-integer segment,
  sieves its odds in a shared-memory bitmap (block-cooperative marking for small
  primes, prime-per-thread for large), then walks the segment accumulating all octal
  statistics in shared counters before one atomic flush. 8^12 sieves in ~420 ms.
- **Octal primitives** (`src/octal_core.h`): `__host__ __device__` digit functions used
  identically by kernels, CPU reference checks, and the CLI.
- **CPU post-pass** (`src/post.h`): one linear scan of the bitmap for gaps, twins, and
  palindromes (inherently sequential bookkeeping).
- **Ground truth** (`src/primality.h`): deterministic 12-base Miller-Rabin for u64.

## Correctness

Built test-first (red → green). `make test` runs ≈ 18.8M checks:

- octal digit functions vs independent string-based references (host *and* device),
- the four divisibility identities verified exhaustively to 300,000,
- predictor ⟺ gcd(n, 210) = 1,
- GPU bitmap compared bit-for-bit against a CPU sieve up to 8^8 (16.7M),
- every statistics counter compared exactly against a CPU recomputation (K ≤ 8),
- π(8^k) against published values; Miller-Rabin against the sieve and against known
  strong-pseudoprime traps (e.g. 3825123056546413051),
- post-pass gaps/twins/palindromes against brute force.

Every survey run additionally hard-fails (exit 3) if its π(8^k) disagrees with the
published table.

## Layout

```
src/octal_core.h   base-8 digit primitives + stats layout (host/device shared)
src/sieve.cuh      fused segmented sieve + octal statistics kernel, pipeline driver
src/post.h         CPU post-pass: gaps, twins, octal palindromic primes
src/primality.h    deterministic Miller-Rabin (u64)
src/main.cu        CLI: survey reports, CSVs, --judge
src/test_main.cu   the full test suite (make test)
src/cpu_survey.cpp GPU-free octal-wheel vs hex-wheel comparative survey
tools/visualize.py renders the prime-wheel figures (SVG, stdlib only)
docs/              sample reports (8^10, 8^12, CPU 1e10) + octal-vs-hex.md + figures
results/           generated reports + CSVs (gitignored)
```

## CUDA vs Rust

The project brief was CUDA if possible, Rust as fallback. CUDA sufficed — the RTX 5090
sieves and pattern-bins 68.7 billion integers in under half a second of kernel time, so
no Rust port was needed.

## Future directions

- Stream segments without a global bitmap to push past 8^12 (memory-free statistics).
- Track the mod-64 (last two octal digits) Chebyshev races at depth.
- Richer octal-digit features (n-grams, positional correlations) as predictors beyond
  the wheel — measuring how much primality signal base-8 digits really carry.
- Octal repunit/repdigit prime hunts.
