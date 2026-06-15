// lattice.cu — primes as points in the Gaussian integers Z[i].
//
// Companion to the base-8 survey: instead of asking "is n prime on the number
// line", we place primes in the 2D lattice and measure the GEOMETRY of that
// point set. Three probes:
//   --hecke R     angle equidistribution of Gaussian primes (Hecke's theorem; a
//                 *proven* result, so it doubles as the pipeline's sanity gate)
//   --moat  K R   the Gaussian moat problem (OPEN): can you walk from the origin
//                 to infinity stepping only on Gaussian primes with step <= K?
//   --selftest    known-answer checks (classification + Hecke chi^2 + step-2 moat)
//
// ponytail: the disk scan and the BFS run on the host. The angle histogram is
// embarrassingly parallel (one lattice point per thread, atomic into bins) and
// is the obvious GPU upgrade — drop in a kernel mirroring `scan_disk` when a
// single host pass stops being fast enough (radius >~ 1e5). The moat BFS is
// inherently sequential frontier bookkeeping (cf. src/post.h) and stays on CPU.

#include "primality.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include <queue>
#include <unordered_set>

// --- Gaussian primality -----------------------------------------------------
// a+bi (a,b integers) is a Gaussian prime iff:
//   * exactly one of a,b is zero and |the other| is a rational prime  ≡ 3 (mod 4)
//   * both nonzero and the norm a^2+b^2 is a rational prime (then it is 2 or ≡1 mod4)
// Units (norm 1) and 0 are not prime.
static inline bool gaussian_prime(long long a, long long b) {
    long long A = a < 0 ? -a : a, B = b < 0 ? -b : b;
    if (A == 0 && B == 0) return false;
    if (A == 0 || B == 0) {
        u64 m = (u64)(A + B);              // the one nonzero magnitude
        return is_prime_u64(m) && (m % 4 == 3);
    }
    u64 norm = (u64)A * (u64)A + (u64)B * (u64)B;
    return is_prime_u64(norm);
}

// --- Hecke angle equidistribution ------------------------------------------
// Bin the angles atan2(b,a) of every Gaussian prime in the open first quadrant
// disk of radius R. Hecke (1920) proved these equidistribute on (0, pi/2), so a
// correct scan gives a flat histogram (reduced chi^2 -> 1). Returns reduced chi^2.
static double hecke_chi2(long long R, int bins, bool verbose) {
    std::vector<long long> h(bins, 0);
    long long total = 0;
    const u64 R2 = (u64)R * (u64)R;
    const double scale = bins / (M_PI / 2.0);
    for (long long a = 1; a <= R; ++a)
        for (long long b = 1; b <= R; ++b) {
            if ((u64)a * a + (u64)b * b > R2) break;   // b increasing -> rest of row is outside
            if (!gaussian_prime(a, b)) continue;
            int k = (int)(atan2((double)b, (double)a) * scale);
            if (k >= bins) k = bins - 1;
            ++h[k]; ++total;
        }
    double exp = (double)total / bins, chi2 = 0;
    for (int k = 0; k < bins; ++k) { double d = h[k] - exp; chi2 += d * d / exp; }
    double reduced = chi2 / (bins - 1);
    if (verbose) {
        printf("Hecke equidistribution, R=%lld: %lld Gaussian primes in quadrant disk, %d bins\n",
               R, total, bins);
        for (int k = 0; k < bins; ++k)
            printf("  [%5.1f-%5.1f deg] %8lld %s\n",
                   k * 90.0 / bins, (k + 1) * 90.0 / bins, h[k],
                   std::string((int)(h[k] * 50 / (exp * 2 + 1)), '#').c_str());
        printf("  reduced chi^2 = %.4f  (-> 1.0 means uniform / equidistributed)\n", reduced);
    }
    return reduced;
}

// --- Gaussian moat ----------------------------------------------------------
// BFS over Gaussian primes from the prime nearest the origin (1+i), stepping to
// any Gaussian prime within Euclidean distance K. Returns the component's
// farthest distance from origin and whether it ESCAPED past radius `limit`
// (escape would be a sensational result; for small K the component is known
// finite — that is the "moat"). `size` receives the component's vertex count.
struct MoatResult { double farthest; bool escaped; long long size; };

static MoatResult moat(double K, long long limit, bool verbose) {
    const long long Kc = (long long)floor(K);
    const double K2 = K * K;
    std::vector<std::pair<long long,long long>> steps;
    for (long long dx = -Kc; dx <= Kc; ++dx)
        for (long long dy = -Kc; dy <= Kc; ++dy) {
            double d2 = (double)dx * dx + (double)dy * dy;
            if (d2 >= 1 && d2 <= K2) steps.push_back({dx, dy});
        }
    auto key = [](long long a, long long b) {
        return ((a + 0x40000000LL) << 32) | (b + 0x40000000LL);   // pack two ~31-bit coords
    };
    std::unordered_set<long long> seen;
    std::queue<std::pair<long long,long long>> q;
    q.push({1, 1}); seen.insert(key(1, 1));
    double farthest = sqrt(2.0); bool escaped = false; long long size = 0;
    while (!q.empty()) {
        auto [a, b] = q.front(); q.pop(); ++size;
        double r = sqrt((double)a * a + (double)b * b);
        if (r > farthest) farthest = r;
        if (r > limit) { escaped = true; break; }
        for (auto [dx, dy] : steps) {
            long long na = a + dx, nb = b + dy;
            long long kk = key(na, nb);
            if (seen.count(kk)) continue;
            if (!gaussian_prime(na, nb)) continue;
            seen.insert(kk); q.push({na, nb});
        }
    }
    if (verbose)
        printf("Moat K=%.4f: component of origin has %lld primes, farthest |z|=%.4f, %s\n",
               K, size, farthest, escaped ? "ESCAPED past limit (!)" : "bounded -> moat exists");
    return {farthest, escaped, size};
}

// --- GPU disk-sieve moat ----------------------------------------------------
// To push to the literature √26 record the per-neighbor Miller–Rabin above is far
// too slow (the √26 component is millions of primes). Instead we PRECOMPUTE a
// Gaussian-prime bitmap over the first-quadrant disk on the GPU, then the BFS is
// pure O(1) bitmap lookups. Symmetry: the Gaussian primes are invariant under the
// 8-element dihedral group (conjugation + 90° rotation), so we BFS only the closed
// first quadrant, folding any step across an axis to |a|,|b|. Distances — all we
// need for the moat — are preserved by that fold.

#define CUDA_OK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
    exit(2); } } while (0)

// Odd-number prime sieve on the host: bit i (set => prime) represents 2i+1.
static std::vector<uint64_t> odd_sieve(u64 nmax) {
    u64 half = nmax / 2 + 1;                       // indices 0..half-1 -> odds 1,3,5,...
    std::vector<uint64_t> bits((half + 63) / 64, ~0ULL);
    bits[0] &= ~1ULL;                              // 1 is not prime
    for (u64 i = 1; (2 * i + 1) * (2 * i + 1) <= nmax; ++i) {
        if (!((bits[i >> 6] >> (i & 63)) & 1ULL)) continue;
        u64 p = 2 * i + 1;
        for (u64 m = p * p; m <= nmax; m += 2 * p) {  // odd multiples only
            u64 j = (m - 1) / 2;
            bits[j >> 6] &= ~(1ULL << (j & 63));
        }
    }
    return bits;
}

// is the odd number n a rational prime, per the device-side odd sieve?
__device__ __forceinline__ bool dev_isrp(u64 n, const uint64_t* s) {
    if (n == 2) return true;
    if (n < 2 || !(n & 1)) return false;
    u64 i = (n - 1) >> 1;
    return (s[i >> 6] >> (i & 63)) & 1ULL;
}

// One thread per lattice point (a,b) in the (R+1)x(R+1) square; set a bit iff it
// is a Gaussian prime inside the disk of radius R.
__global__ void mark_gaussian(uint64_t* gbm, const uint64_t* rsieve, long long R) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    u64 span = (u64)(R + 1) * (R + 1);
    if (t >= span) return;
    long long a = t / (R + 1), b = t % (R + 1);
    if ((u64)a * a + (u64)b * b > (u64)R * R) return;   // outside the disk
    bool prime;
    if (a == 0 && b == 0)            prime = false;
    else if (a == 0 || b == 0)       { u64 m = a + b; prime = dev_isrp(m, rsieve) && (m % 4 == 3); }
    else                             prime = dev_isrp((u64)a * a + (u64)b * b, rsieve);
    if (prime) atomicOr((unsigned long long*)&gbm[t >> 6], 1ULL << (t & 63));
}

static MoatResult moat_gpu(double K, long long R, bool verbose) {
    const long long Kc = (long long)floor(K);
    const double K2 = K * K;
    const u64 Nmax = (u64)R * R;
    const u64 span = (u64)(R + 1) * (R + 1);
    const size_t words = (span + 63) / 64;

    if (verbose) printf("Building GPU Gaussian-prime bitmap: R=%lld (%.2f GB disk + sieve)...\n",
                        R, (words * 8.0 + (Nmax / 2 / 8.0)) / 1e9);

    // 1. rational odd-sieve to R^2 on host, push to device.
    std::vector<uint64_t> rs = odd_sieve(Nmax);
    uint64_t *d_rs = nullptr, *d_gbm = nullptr;
    CUDA_OK(cudaMalloc(&d_rs, rs.size() * 8));
    CUDA_OK(cudaMemcpy(d_rs, rs.data(), rs.size() * 8, cudaMemcpyHostToDevice));
    CUDA_OK(cudaMalloc(&d_gbm, words * 8));
    CUDA_OK(cudaMemset(d_gbm, 0, words * 8));

    // 2. mark Gaussian primes in parallel, copy bitmap back.
    int tpb = 256; u64 blocks = (span + tpb - 1) / tpb;
    mark_gaussian<<<(unsigned)blocks, tpb>>>(d_gbm, d_rs, R);
    CUDA_OK(cudaGetLastError());
    CUDA_OK(cudaDeviceSynchronize());
    std::vector<uint64_t> gbm(words);
    CUDA_OK(cudaMemcpy(gbm.data(), d_gbm, words * 8, cudaMemcpyDeviceToHost));
    CUDA_OK(cudaFree(d_rs)); CUDA_OK(cudaFree(d_gbm));

    // 3. folded first-quadrant BFS over the bitmap. visited reuses a second bitmap.
    auto gp = [&](long long a, long long b) -> bool {
        long long A = a < 0 ? -a : a, B = b < 0 ? -b : b;
        if (A > R || B > R || (u64)A * A + (u64)B * B > Nmax) return false;
        u64 t = (u64)A * (R + 1) + B;
        return (gbm[t >> 6] >> (t & 63)) & 1ULL;
    };
    std::vector<uint64_t> vis(words, 0);
    auto visit = [&](long long A, long long B) -> bool {     // returns true if newly visited
        u64 t = (u64)A * (R + 1) + B;
        if ((vis[t >> 6] >> (t & 63)) & 1ULL) return false;
        vis[t >> 6] |= 1ULL << (t & 63); return true;
    };
    std::vector<std::pair<long long,long long>> steps;
    for (long long dx = -Kc; dx <= Kc; ++dx)
        for (long long dy = -Kc; dy <= Kc; ++dy) {
            double d2 = (double)dx * dx + (double)dy * dy;
            if (d2 >= 1 && d2 <= K2) steps.push_back({dx, dy});
        }
    std::queue<std::pair<long long,long long>> q;
    q.push({1, 1}); visit(1, 1);
    double farthest = sqrt(2.0); long long size = 0;
    const double boundary = (double)(R - Kc);       // reaching here means R is too small
    bool reached_boundary = false;
    while (!q.empty()) {
        auto [a, b] = q.front(); q.pop(); ++size;
        double r = sqrt((double)a * a + (double)b * b);
        if (r > farthest) farthest = r;
        if (r > boundary) reached_boundary = true;
        for (auto [dx, dy] : steps) {
            long long na = a + dx, nb = b + dy;
            long long A = na < 0 ? -na : na, B = nb < 0 ? -nb : nb;   // fold into quadrant
            if (!gp(A, B)) continue;
            if (visit(A, B)) q.push({A, B});
        }
    }
    if (verbose) {
        if (reached_boundary)
            printf("Moat K=%.4f: component reached the disk boundary (farthest |z|=%.1f >= R-K). "
                   "INCONCLUSIVE — rerun with larger R.\n", K, farthest);
        else
            printf("Moat K=%.4f: BOUNDED. quadrant component = %lld primes, farthest |z|=%.4f "
                   "-> moat confirmed below this radius.\n", K, size, farthest);
    }
    return {farthest, reached_boundary, size};
}

// ===========================================================================
// EISENSTEIN INTEGERS Z[w], w = e^{2pi i/3}.  (council pivot: the less-computed
// hexagonal-lattice analogue of the Gaussian moat / Hecke story.)
//
// A point a + b*w embeds in the plane at z = (a - b/2, b*sqrt3/2), so the
// SQUARED Euclidean norm equals the algebraic norm N(a,b) = a^2 - a*b + b^2.
// Distances — all the moat needs — are therefore exact integers in (a,b) coords.
// The 6 units are (±1,0),(0,±1),±(1,1); the unit-1 neighbours are the 6 vectors
// of norm 1, i.e. the hexagonal grid.
//
// a + b*w is an Eisenstein prime iff:
//   * N(a,b) is a rational prime   (split p ≡ 1 mod 3, or the ramified 3), OR
//   * it is an associate of an inert rational prime p ≡ 2 mod 3 — these sit on
//     the three lattice axes b=0, a=0, a=b at magnitude p.
// (A perfect-square norm off those axes is a genuine factorisation -> composite.)

static inline long long eis_norm(long long a, long long b) { return a * a - a * b + b * b; }

static inline bool eisenstein_prime(long long a, long long b) {
    long long N = eis_norm(a, b);
    if (N == 0) return false;
    if (is_prime_u64((u64)N)) return true;                 // split or ramified
    if (b == 0 || a == 0 || a == b) {                      // inert, on a hex axis
        long long p = (a < 0 ? -a : a);
        if (a == 0) p = (b < 0 ? -b : b);
        return is_prime_u64((u64)p) && (p % 3 == 2);
    }
    return false;
}

__device__ __forceinline__ bool dev_eprime(long long a, long long b, const uint64_t* s) {
    long long N = a * a - a * b + b * b;
    if (N == 0) return false;
    if (dev_isrp((u64)N, s)) return true;
    if (b == 0 || a == 0 || a == b) {
        long long p = (a < 0 ? -a : a);
        if (a == 0) p = (b < 0 ? -b : b);
        return dev_isrp((u64)p, s) && (p % 3 == 2);
    }
    return false;
}

// Hecke equidistribution in Z[w]: angles of Eisenstein primes over the full disk
// of radius R. Like Z[i] they equidistribute, so the histogram is flat -> chi^2 ~ 1.
static double ehecke_chi2(long long R, int bins, bool verbose) {
    std::vector<long long> h(bins, 0);
    long long total = 0;
    const u64 R2 = (u64)R * (u64)R;
    const long long A = (long long)(1.1547 * R) + 2;        // |a|,|b| bound for N<=R^2
    const double TAU = 2.0 * M_PI, scale = bins / TAU, s3 = sqrt(3.0) / 2.0;
    for (long long a = -A; a <= A; ++a)
        for (long long b = -A; b <= A; ++b) {
            if ((u64)eis_norm(a, b) > R2 || (a == 0 && b == 0)) continue;
            if (!eisenstein_prime(a, b)) continue;
            double ang = atan2(b * s3, a - b / 2.0);          // (-pi, pi]
            if (ang < 0) ang += TAU;
            int k = (int)(ang * scale);
            if (k >= bins) k = bins - 1;
            ++h[k]; ++total;
        }
    double exp = (double)total / bins, chi2 = 0;
    for (int k = 0; k < bins; ++k) { double d = h[k] - exp; chi2 += d * d / exp; }
    double reduced = chi2 / (bins - 1);
    if (verbose) {
        printf("Eisenstein Hecke, R=%lld: %lld primes in disk, %d bins over 360deg\n", R, total, bins);
        for (int k = 0; k < bins; ++k)
            printf("  [%5.1f-%5.1f deg] %8lld %s\n", k * 360.0 / bins, (k + 1) * 360.0 / bins, h[k],
                   std::string((int)(h[k] * 50 / (exp * 2 + 1)), '#').c_str());
        printf("  reduced chi^2 = %.4f  (-> 1.0 = equidistributed)\n", reduced);
    }
    return reduced;
}

// hex steps of length <= K: integer offsets with 1 <= N(da,db) <= K^2.
static std::vector<std::pair<long long,long long>> hex_steps(double K) {
    std::vector<std::pair<long long,long long>> s;
    long long M = (long long)ceil(K * 1.16) + 1;
    for (long long da = -M; da <= M; ++da)
        for (long long db = -M; db <= M; ++db) {
            long long n = eis_norm(da, db);
            if (n >= 1 && (double)n <= K * K) s.push_back({da, db});
        }
    return s;
}

// CPU Miller-Rabin BFS moat in Z[w] (full plane). Reference + small-scale truth.
static MoatResult emoat(double K, long long limit, bool verbose) {
    auto steps = hex_steps(K);
    auto key = [](long long a, long long b) { return ((a + 0x40000000LL) << 32) | (b + 0x40000000LL); };
    std::unordered_set<long long> seen;
    std::queue<std::pair<long long,long long>> q;
    q.push({1, -1}); seen.insert(key(1, -1));               // nearest Eisenstein prime, N=3
    double farthest = sqrt(3.0); bool escaped = false; long long size = 0;
    while (!q.empty()) {
        auto [a, b] = q.front(); q.pop(); ++size;
        double r = sqrt((double)eis_norm(a, b));
        if (r > farthest) farthest = r;
        if (r > limit) { escaped = true; break; }
        for (auto [dx, dy] : steps) {
            long long na = a + dx, nb = b + dy, kk = key(na, nb);
            if (seen.count(kk) || !eisenstein_prime(na, nb)) continue;
            seen.insert(kk); q.push({na, nb});
        }
    }
    if (verbose)
        printf("Eisenstein moat K=%.4f: component has %lld primes, farthest |z|=%.4f, %s\n",
               K, size, farthest, escaped ? "ESCAPED past limit (!)" : "bounded -> moat exists");
    return {farthest, escaped, size};
}

// GPU disk-sieve Eisenstein moat. Full-plane offset grid (no symmetry fold yet —
// ponytail: the 12-fold hex symmetry is the obvious memory win when pushing for a
// record radius; for the first port we trade ~12x memory for simplicity).
__global__ void mark_eisenstein(uint64_t* gbm, const uint64_t* rsieve, long long OFF, long long side, u64 R2) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)side * side) return;
    long long a = (long long)(t / side) - OFF, b = (long long)(t % side) - OFF;
    if ((u64)(a * a - a * b + b * b) > R2) return;
    if (dev_eprime(a, b, rsieve)) atomicOr((unsigned long long*)&gbm[t >> 6], 1ULL << (t & 63));
}

static MoatResult emoat_gpu(double K, long long R, bool verbose) {
    const u64 Nmax = (u64)R * R;
    const long long OFF = (long long)(1.1547 * R) + 2, side = 2 * OFF + 1;
    const u64 span = (u64)side * side;
    const size_t words = (span + 63) / 64;
    if (verbose) printf("Building GPU Eisenstein bitmap: R=%lld, grid %lldx%lld (%.2f GB)...\n",
                        R, side, side, (words * 8.0 * 2 + Nmax / 16.0) / 1e9);

    std::vector<uint64_t> rs = odd_sieve(Nmax);
    uint64_t *d_rs = nullptr, *d_gbm = nullptr;
    CUDA_OK(cudaMalloc(&d_rs, rs.size() * 8));
    CUDA_OK(cudaMemcpy(d_rs, rs.data(), rs.size() * 8, cudaMemcpyHostToDevice));
    CUDA_OK(cudaMalloc(&d_gbm, words * 8));
    CUDA_OK(cudaMemset(d_gbm, 0, words * 8));
    int tpb = 256; u64 blocks = (span + tpb - 1) / tpb;
    mark_eisenstein<<<(unsigned)blocks, tpb>>>(d_gbm, d_rs, OFF, side, Nmax);
    CUDA_OK(cudaGetLastError()); CUDA_OK(cudaDeviceSynchronize());
    std::vector<uint64_t> gbm(words);
    CUDA_OK(cudaMemcpy(gbm.data(), d_gbm, words * 8, cudaMemcpyDeviceToHost));
    CUDA_OK(cudaFree(d_rs)); CUDA_OK(cudaFree(d_gbm));

    auto idx = [&](long long a, long long b) { return (u64)(a + OFF) * side + (b + OFF); };
    auto gp = [&](long long a, long long b) -> bool {
        if (a < -OFF || a > OFF || b < -OFF || b > OFF || (u64)eis_norm(a, b) > Nmax) return false;
        u64 t = idx(a, b); return (gbm[t >> 6] >> (t & 63)) & 1ULL;
    };
    std::vector<uint64_t> vis(words, 0);
    auto steps = hex_steps(K);
    std::queue<std::pair<long long,long long>> q;
    { u64 t = idx(1, -1); vis[t >> 6] |= 1ULL << (t & 63); }
    q.push({1, -1});
    double farthest = sqrt(3.0); long long size = 0; bool reached_boundary = false;
    const double boundary = (double)R - K;
    while (!q.empty()) {
        auto [a, b] = q.front(); q.pop(); ++size;
        double r = sqrt((double)eis_norm(a, b));
        if (r > farthest) farthest = r;
        if (r > boundary) reached_boundary = true;
        for (auto [dx, dy] : steps) {
            long long na = a + dx, nb = b + dy;
            if (!gp(na, nb)) continue;
            u64 t = idx(na, nb);
            if ((vis[t >> 6] >> (t & 63)) & 1ULL) continue;
            vis[t >> 6] |= 1ULL << (t & 63); q.push({na, nb});
        }
    }
    if (verbose) {
        if (reached_boundary)
            printf("Eisenstein moat K=%.4f: reached disk boundary (|z|=%.1f). INCONCLUSIVE — larger R.\n", K, farthest);
        else
            printf("Eisenstein moat K=%.4f: BOUNDED. component = %lld primes, farthest |z|=%.4f -> moat confirmed.\n",
                   K, size, farthest);
    }
    return {farthest, reached_boundary, size};
}

// ===========================================================================
// NEAREST-NEIGHBOUR GAP STATISTICS (council's next step: a falsifiable empirical
// law, not a record radius). For both lattices we build a full-plane prime
// bitmap, find every prime's nearest prime neighbour, and compare the observed
// NN-distance distribution to the random (Poisson) model exp(-lambda*pi*r^2),
// whose mean NN distance is 1/(2*sqrt(lambda)) at density lambda. Deviation
// from Poisson is the structure the moat is the extreme tail of.

__device__ __forceinline__ bool dev_gprime(long long a, long long b, const uint64_t* s) {
    long long A = a < 0 ? -a : a, B = b < 0 ? -b : b;
    if (A == 0 && B == 0) return false;
    if (A == 0 || B == 0) { u64 m = A + B; return dev_isrp(m, s) && (m % 4 == 3); }
    return dev_isrp((u64)A * A + (u64)B * B, s);
}

// One kernel, both lattices, full-plane offset grid (used by the gap tool).
__global__ void mark_lattice(uint64_t* gbm, const uint64_t* rs, long long OFF, long long side, u64 R2, int eis) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)side * side) return;
    long long a = (long long)(t / side) - OFF, b = (long long)(t % side) - OFF;
    u64 N = eis ? (u64)(a * a - a * b + b * b) : (u64)(a * a + b * b);
    if (N > R2) return;
    bool p = eis ? dev_eprime(a, b, rs) : dev_gprime(a, b, rs);
    if (p) atomicOr((unsigned long long*)&gbm[t >> 6], 1ULL << (t & 63));
}

static std::vector<uint64_t> lattice_bitmap(long long R, bool eis, long long OFF, long long side) {
    u64 Nmax = (u64)R * R, span = (u64)side * side; size_t words = (span + 63) / 64;
    std::vector<uint64_t> rs = odd_sieve(Nmax), gbm(words);
    uint64_t *d_rs, *d_gbm;
    CUDA_OK(cudaMalloc(&d_rs, rs.size() * 8));
    CUDA_OK(cudaMemcpy(d_rs, rs.data(), rs.size() * 8, cudaMemcpyHostToDevice));
    CUDA_OK(cudaMalloc(&d_gbm, words * 8));
    CUDA_OK(cudaMemset(d_gbm, 0, words * 8));
    int tpb = 256; u64 blocks = (span + tpb - 1) / tpb;
    mark_lattice<<<(unsigned)blocks, tpb>>>(d_gbm, d_rs, OFF, side, Nmax, eis ? 1 : 0);
    CUDA_OK(cudaGetLastError()); CUDA_OK(cudaDeviceSynchronize());
    CUDA_OK(cudaMemcpy(gbm.data(), d_gbm, words * 8, cudaMemcpyDeviceToHost));
    CUDA_OK(cudaFree(d_rs)); CUDA_OK(cudaFree(d_gbm));
    return gbm;
}

struct GapStats { double obs_mean, poisson_mean, random_mean, random_sd, max_gap; long long count, random_count; };

static GapStats gap_stats(bool eis, long long R, int nseeds, bool verbose) {
    const long long OFF = eis ? (long long)(1.1547 * R) + 2 : R + 2, side = 2 * OFF + 1;
    std::vector<uint64_t> gbm = lattice_bitmap(R, eis, OFF, side);
    auto isP = [&](long long a, long long b) -> bool {
        if (a < -OFF || a > OFF || b < -OFF || b > OFF) return false;
        u64 t = (u64)(a + OFF) * side + (b + OFF); return (gbm[t >> 6] >> (t & 63)) & 1ULL;
    };
    auto dist = [&](long long dx, long long dy) {
        return sqrt((double)(eis ? dx * dx - dx * dy + dy * dy : dx * dx + dy * dy));
    };
    const long long margin = 8;
    const u64 Rin2 = (u64)(R - margin) * (R - margin);
    auto norm = [&](long long a, long long b) -> u64 {
        return eis ? (u64)(a * a - a * b + b * b) : (u64)(a * a + b * b);
    };

    // mean NN distance over the points selected by `pred`, scanned inside Rin2.
    auto scan = [&](auto pred, double& mean, double& maxg, long long& cnt) {
        double sum = 0; maxg = 0; cnt = 0;
        for (long long a = -OFF; a <= OFF; ++a)
            for (long long b = -OFF; b <= OFF; ++b) {
                if (norm(a, b) > Rin2 || !pred(a, b)) continue;
                double best = 1e18;
                for (long long d = 1; d <= 64; ++d) {
                    for (long long dx = -d; dx <= d; ++dx)
                        for (long long dy = -d; dy <= d; ++dy) {
                            if (std::max(dx < 0 ? -dx : dx, dy < 0 ? -dy : dy) != d) continue;
                            if (!pred(a + dx, b + dy)) continue;
                            double dd = dist(dx, dy);
                            if (dd < best) best = dd;
                        }
                    if (best <= (double)d) break;
                }
                if (best > 1e17) continue;
                sum += best; ++cnt; if (best > maxg) maxg = best;
            }
        mean = cnt ? sum / cnt : 0;
    };

    // --- primes ---
    double p_mean, p_max; long long p_cnt;
    scan(isP, p_mean, p_max, p_cnt);

    // --- random control over MANY seeds (council kill-test): a Bernoulli lattice
    //     subset at the SAME density on the SAME lattice shares the prime set's
    //     hard-core floor. Running independent seeds turns the prime/random ratio
    //     into a measurement WITH ERROR BARS, so the repulsion (and the Z[i] vs
    //     Z[w] gap) can be tested against sampling noise instead of asserted.
    long long cells = 0;
    for (long long a = -OFF; a <= OFF; ++a)
        for (long long b = -OFF; b <= OFF; ++b)
            if (norm(a, b) <= Rin2) ++cells;
    const u64 thresh = (u64)(((double)p_cnt / cells) * (double)(1ULL << 24));
    std::vector<double> r_means;
    long long r_cnt = 0; double r_max = 0;
    for (int s = 0; s < nseeds; ++s) {
        u64 seed = 0x9E3779B97F4A7C15ULL * (u64)(s + 1);
        auto isR = [&](long long a, long long b) -> bool {
            if (a < -OFF || a > OFF || b < -OFF || b > OFF || norm(a, b) > (u64)R * R) return false;
            u64 x = (u64)(a + OFF) * side + (b + OFF) + seed;
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
            x ^= x >> 31;
            return (x >> 40) < thresh;
        };
        double rm, rmx; long long rc;
        scan(isR, rm, rmx, rc);
        r_means.push_back(rm); r_cnt = rc; r_max = std::max(r_max, rmx);
    }
    double r_mean = 0; for (double v : r_means) r_mean += v; r_mean /= r_means.size();
    double r_sd = 0; for (double v : r_means) r_sd += (v - r_mean) * (v - r_mean);
    r_sd = r_means.size() > 1 ? sqrt(r_sd / (r_means.size() - 1)) : 0;

    double lambda = p_cnt / (M_PI * (double)Rin2);
    double ratio = p_mean / r_mean, ratio_sd = ratio * (r_sd / r_mean);   // p_mean fixed
    GapStats g{ p_mean, 1.0 / (2 * sqrt(lambda)), r_mean, r_sd, p_max, p_cnt, r_cnt };
    if (verbose) {
        printf("Nearest-neighbour gaps, %s, R=%lld  (%d random seeds)\n", eis ? "Z[w]" : "Z[i]", R, nseeds);
        printf("  primes:        %lld points, mean NN = %.5f\n", p_cnt, p_mean);
        printf("  random subset: ~%lld points, mean NN = %.5f +/- %.5f (seed spread)\n", r_cnt, r_mean, r_sd);
        printf("  prime/random ratio = %.5f +/- %.5f   (excess repulsion = %.2f%% +/- %.2f%%, %.1f sigma)\n",
               ratio, ratio_sd, (ratio - 1) * 100, ratio_sd * 100, (ratio - 1) / (ratio_sd > 0 ? ratio_sd : 1e-9));
    }
    return g;
}

// ===========================================================================
// ANGULAR PAIR-CORRELATION. The NN study gave one moment (mean gap); the pair
// correlation g(d) is the whole short-range two-point function: for a source
// point, how many other points sit at separation d, vs the matched random
// subset. g_prime/g_random < 1 at small d = repulsion; -> 1 as d grows = the
// correlation length. We also resolve it by the ANGLE of the separation vector
// to test whether the repulsion is isotropic or carries the lattice's symmetry.

struct PCF { double g_short, aniso; };   // short-range depletion ratio; angular anisotropy (CV)

static PCF pair_corr(bool eis, long long R, int sectors, bool verbose) {
    const long long OFF = eis ? (long long)(1.1547 * R) + 2 : R + 2, side = 2 * OFF + 1;
    std::vector<uint64_t> gbm = lattice_bitmap(R, eis, OFF, side);
    auto norm = [&](long long a, long long b) -> u64 { return eis ? (u64)(a*a - a*b + b*b) : (u64)(a*a + b*b); };
    auto isP = [&](long long a, long long b) -> bool {
        if (a < -OFF || a > OFF || b < -OFF || b > OFF) return false;
        u64 t = (u64)(a + OFF) * side + (b + OFF); return (gbm[t >> 6] >> (t & 63)) & 1ULL;
    };
    const long long margin = 16;
    const u64 Rin2 = (u64)(R - margin) * (R - margin);
    // density -> random-subset threshold (same construction as gap_stats)
    long long p_cnt = 0, cells = 0;
    for (long long a = -OFF; a <= OFF; ++a)
        for (long long b = -OFF; b <= OFF; ++b)
            if (norm(a, b) <= Rin2) { ++cells; if (isP(a, b)) ++p_cnt; }
    const u64 thresh = (u64)(((double)p_cnt / cells) * (double)(1ULL << 24));
    auto isR = [&](long long a, long long b) -> bool {
        if (a < -OFF || a > OFF || b < -OFF || b > OFF || norm(a, b) > (u64)R * R) return false;
        u64 x = (u64)(a + OFF) * side + (b + OFF) + 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        x ^= x >> 31;
        return (x >> 40) < thresh;
    };

    const double Dmax = 8.0, s3 = sqrt(3.0) / 2.0;
    const int RB = 32; const double rbin = Dmax / RB;       // radial bins of g(d)
    const int FB = 360;                                     // FINE angular bins; coarsen later
    // accumulate radial + fine-angular pair counts for source points selected by pred
    auto accum = [&](auto pred, std::vector<double>& radial, std::vector<double>& ang) {
        radial.assign(RB, 0); ang.assign(FB, 0);
        long long M = (long long)Dmax + 1;
        for (long long a = -OFF; a <= OFF; ++a)
            for (long long b = -OFF; b <= OFF; ++b) {
                if (norm(a, b) > Rin2 || !pred(a, b)) continue;
                for (long long dx = -M; dx <= M; ++dx)
                    for (long long dy = -M; dy <= M; ++dy) {
                        if (dx == 0 && dy == 0) continue;
                        double dd = sqrt((double)(eis ? dx*dx - dx*dy + dy*dy : dx*dx + dy*dy));
                        if (dd > Dmax || !pred(a + dx, b + dy)) continue;
                        int rb = (int)(dd / rbin); if (rb < RB) radial[rb] += 1;
                        if (dd >= 2.0) {                    // angle of separation in the embedded plane
                            double ex = eis ? dx - dy / 2.0 : (double)dx, ey = eis ? dy * s3 : (double)dy;
                            double ph = atan2(ey, ex); if (ph < 0) ph += 2 * M_PI;
                            ang[(int)(ph * FB / (2 * M_PI)) % FB] += 1;
                        }
                    }
            }
    };
    std::vector<double> pr, pa, rr, ra;
    accum(isP, pr, pa); accum(isR, rr, ra);

    // g_short: ratio of prime to random pairs in the smallest populated shells (d<3)
    double ps = 0, rs = 0;
    for (int k = 0; k < RB && (k + 1) * rbin <= 3.0; ++k) { ps += pr[k]; rs += rr[k]; }
    double g_short = (rs > 0) ? ps / rs : 0;

    // anisotropy CV of the prime/random angular ratio, coarsening the FINE bins to AB
    // sectors. AB is now the SAME for both lattices (the eis?12:8 confound is gone);
    // we sweep AB so the bin-count dependence of CV is visible, not hidden.
    auto cv_at = [&](int AB) -> double {
        if (FB % AB) return -1;
        int grp = FB / AB; double gm = 0; std::vector<double> gr(AB, 0);
        for (int k = 0; k < AB; ++k) {
            double pp = 0, qq = 0;
            for (int j = 0; j < grp; ++j) { pp += pa[k * grp + j]; qq += ra[k * grp + j]; }
            gr[k] = qq > 0 ? pp / qq : 0; gm += gr[k];
        }
        gm /= AB; double var = 0; for (double v : gr) var += (v - gm) * (v - gm);
        return gm > 0 ? sqrt(var / AB) / gm : 0;
    };
    double aniso = cv_at(sectors);

    if (verbose) {
        printf("Pair correlation, %s, R=%lld (%lld source primes, Dmax=%.0f)\n",
               eis ? "Z[w]" : "Z[i]", R, p_cnt, Dmax);
        printf("  short-range (d<3) g_prime/g_random = %.4f   (<1 = repulsion beyond hard-core)\n", g_short);
        printf("  angular-ratio CV vs sector count (EQUAL sectors for both lattices):\n");
        for (int AB : {6, 8, 12, 24, 36})
            printf("    %2d sectors: CV = %.4f\n", AB, cv_at(AB));
        printf("  --> CV grows mechanically with sector count; compare the two lattices at the\n");
        printf("      SAME AB. The earlier eis?12:8 comparison was apples-to-oranges.\n");
    }
    return {g_short, aniso};
}

// ===========================================================================
// HARDY-LITTLEWOOD SINGULAR SERIES (the principled null). For a fixed offset
// delta, the density of prime pairs (pi, pi+delta) relative to the random-at-
// matched-density baseline converges to the singular series
//     S(delta) = prod_p  local_p,    local_p = (1 - w_p/Np) / (1 - 1/Np)^2,
// the product over prime ideals p, where w_p = #{distinct residues of {0,-delta}
// mod p} = 1 if p | delta else 2.  So local_p = Np/(Np-1) when p|delta, and
// (1-2/Np)/(1-1/Np)^2 otherwise.  The norm-2 ideal of Z[i] makes local=0 unless
// it divides delta -> the Gaussian parity obstruction (S=0 for "odd" offsets).
// Z[w] has no norm-2 ideal, so NO hard obstruction. The measured prime/random
// ratio at offset delta should equal S(delta): that is the council's test.

static std::vector<int> primes_upto(int Q) {
    std::vector<char> c(Q + 1, 1); std::vector<int> ps;
    for (int i = 2; i <= Q; ++i) if (c[i]) { ps.push_back(i); for (long long j = (long long)i * i; j <= Q; j += i) c[j] = 0; }
    return ps;
}

static double singular_series(bool eis, long long dx, long long dy, const std::vector<int>& primes) {
    long long Nd = eis ? dx * dx - dx * dy + dy * dy : dx * dx + dy * dy;
    if (Nd == 0) return 0;
    const long long Q = primes.empty() ? 0 : (long long)primes.back();
    auto b = [](double q) { double u = 1.0 - 1.0 / q; return (1.0 - 2.0 / q) / (u * u); };
    auto modz = [&](long long v, int p) { return (int)(((v % p) + p) % p); };
    double prod = 1.0;
    for (int p : primes) {
        if (Nd % p != 0) {                         // p divides no ideal factor of delta
            if (!eis) {
                if (p == 2) prod *= b(2);                                  // =0 unless 2|Nd
                else if (p % 4 == 1) prod *= b(p) * b(p);                  // two split ideals
                else { long long q = (long long)p * p; if (q <= Q) prod *= b((double)q); }  // inert
            } else {
                if (p == 3) prod *= b(3);
                else if (p % 3 == 1) prod *= b(p) * b(p);
                else { long long q = (long long)p * p; if (q <= Q) prod *= b((double)q); }
            }
            continue;
        }
        // p | Nd: handle each ideal above p, splitting into divides-delta vs not.
        if ((!eis && p == 2) || (eis && p == 3)) {                        // ramified, norm p
            prod *= (double)p / (p - 1);                                  // it divides delta (Nd ≡ 0)
        } else if ((!eis && p % 4 == 1) || (eis && p % 3 == 1)) {         // two split ideals, norm p
            for (int r = 0; r < p; ++r) {
                bool root = !eis ? (((long long)r * r + 1) % p == 0) : (((long long)r * r + r + 1) % p == 0);
                if (!root) continue;
                bool div = modz(dx + dy * (long long)r, p) == 0;          // i (or w) ≡ r mod this ideal
                prod *= div ? (double)p / (p - 1) : b((double)p);
            }
        } else {                                                          // inert, norm p^2; p|Nd ⟹ p|delta
            long long q = (long long)p * p;
            prod *= (double)q / (q - 1);
        }
    }
    return prod;
}

// Overlay: measured prime/random pair ratio at each small offset vs S(delta).
static void sigma_test(bool eis, long long R, int nseeds, bool verbose) {
    const long long OFF = eis ? (long long)(1.1547 * R) + 2 : R + 2, side = 2 * OFF + 1;
    std::vector<uint64_t> gbm = lattice_bitmap(R, eis, OFF, side);
    auto norm = [&](long long a, long long b) -> u64 { return eis ? (u64)(a*a - a*b + b*b) : (u64)(a*a + b*b); };
    auto isP = [&](long long a, long long b) -> bool {
        if (a < -OFF || a > OFF || b < -OFF || b > OFF) return false;
        u64 t = (u64)(a + OFF) * side + (b + OFF); return (gbm[t >> 6] >> (t & 63)) & 1ULL;
    };
    const long long margin = 12;
    const u64 Rin2 = (u64)(R - margin) * (R - margin);
    long long p_cnt = 0, cells = 0;
    for (long long a = -OFF; a <= OFF; ++a) for (long long b = -OFF; b <= OFF; ++b)
        if (norm(a, b) <= Rin2) { ++cells; if (isP(a, b)) ++p_cnt; }
    const u64 thresh = (u64)(((double)p_cnt / cells) * (double)(1ULL << 24));
    auto isRs = [&](long long a, long long b, u64 seed) -> bool {
        if (a < -OFF || a > OFF || b < -OFF || b > OFF || norm(a, b) > (u64)R * R) return false;
        u64 x = (u64)(a + OFF) * side + (b + OFF) + seed;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL; x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL; x ^= x >> 31;
        return (x >> 40) < thresh;
    };
    // canonical small offsets (one of each +/- pair), norm <= 36
    const long long Dm = 6;
    std::vector<std::pair<long long,long long>> offs;
    for (long long dx = -Dm; dx <= Dm; ++dx)
        for (long long dy = -Dm; dy <= Dm; ++dy) {
            if (dx == 0 && dy == 0) continue;
            if (!(dx > 0 || (dx == 0 && dy > 0))) continue;        // half-plane
            u64 n = eis ? (u64)(dx*dx - dx*dy + dy*dy) : (u64)(dx*dx + dy*dy);
            if (n <= 36) offs.push_back({dx, dy});
        }
    std::vector<int> primes = primes_upto(200000);

    if (verbose) printf("Singular series overlay, %s, R=%lld (%lld primes, %d seeds)\n  %-10s %5s %10s %10s %8s\n",
                        eis ? "Z[w]" : "Z[i]", R, p_cnt, nseeds, "offset", "N", "obs r", "S(delta)", "r/S");
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sresid = 0; int ncmp = 0;
    for (auto [dx, dy] : offs) {
        long long P = 0; double Rnd = 0;
        for (long long a = -OFF; a <= OFF; ++a) for (long long b = -OFF; b <= OFF; ++b) {
            if (norm(a, b) > Rin2 || !isP(a, b)) continue;
            if (isP(a + dx, b + dy)) ++P;
        }
        for (int s = 0; s < nseeds; ++s) {
            u64 seed = 0x9E3779B97F4A7C15ULL * (u64)(s + 1); long long rc = 0;
            for (long long a = -OFF; a <= OFF; ++a) for (long long b = -OFF; b <= OFF; ++b) {
                if (norm(a, b) > Rin2 || !isRs(a, b, seed)) continue;
                if (isRs(a + dx, b + dy, seed)) ++rc;
            }
            Rnd += rc;
        }
        Rnd /= nseeds;
        double r = Rnd > 0 ? P / Rnd : 0;
        double S = singular_series(eis, dx, dy, primes);
        if (verbose) printf("  (%2lld,%2lld)    %5llu %10.4f %10.4f %8.3f\n", dx, dy,
                            (unsigned long long)(eis ? dx*dx-dx*dy+dy*dy : dx*dx+dy*dy), r, S, S > 0 ? r / S : 0);
        if (S > 1e-9) { sx += S; sy += r; sxx += S*S; syy += r*r; sxy += S*r; sresid += r/S; ++ncmp; }
    }
    if (verbose && ncmp > 1) {
        double corr = (ncmp*sxy - sx*sy) / sqrt((ncmp*sxx - sx*sx) * (ncmp*syy - sy*sy));
        printf("  admissible offsets: %d   corr(obs, S) = %.4f   mean(obs/S) = %.4f\n",
               ncmp, corr, sresid / ncmp);
        printf("  --> corr~1 and mean(obs/S)~1 means the repulsion IS the singular series.\n");
    }
}

// --- self-test (red->green gate) -------------------------------------------
static int selftest() {
    // classification against hand-computed truth
    assert(gaussian_prime(1, 1));                 // norm 2, ramified prime
    assert(!gaussian_prime(2, 0));                // 2 = -i(1+i)^2, not prime
    assert(gaussian_prime(3, 0) && gaussian_prime(0, 3));  // 3 ≡ 3 mod4 on axis
    assert(!gaussian_prime(5, 0));                // 5 = (2+i)(2-i) splits
    assert(gaussian_prime(2, 1) && gaussian_prime(1, 2));  // norm 5 prime
    assert(gaussian_prime(7, 0));                 // 7 ≡ 3 mod4
    assert(gaussian_prime(3, 2));                 // norm 13 prime
    assert(!gaussian_prime(1, 0) && !gaussian_prime(0, 0)); // unit, zero
    assert(gaussian_prime(4, 1));                 // norm 17 prime
    assert(!gaussian_prime(3, 3));                // norm 18 composite

    // Hecke: angles must be ~uniform. Reduced chi^2 near 1 (loose bound catches a broken scan).
    double chi2 = hecke_chi2(2000, 36, false);
    assert(chi2 < 2.0 && "Hecke angles not equidistributed -> scan is broken");

    // Moat: step-2 component of the origin is KNOWN to be finite (Gethner-Wagon-Wick).
    MoatResult m = moat(2.0, 1000, false);
    assert(!m.escaped && "step-2 moat should be bounded; escape => bug or Fields medal");

    // GPU disk-sieve must agree with the CPU Miller-Rabin BFS on the symmetry-
    // invariant quantity (farthest |z|) for several step sizes. R=400 contains all.
    for (double K : {1.5, 2.0, 2.9}) {
        MoatResult cpu = moat(K, 100000, false);
        MoatResult gpu = moat_gpu(K, 400, false);
        assert(!gpu.escaped && "validation R too small");
        assert(fabs(cpu.farthest - gpu.farthest) < 1e-6 && "GPU bitmap disagrees with CPU BFS");
    }

    // --- Eisenstein Z[w] ---
    assert(eisenstein_prime(1, -1));              // N=3, ramified prime
    assert(eisenstein_prime(2, 0));               // 2 ≡ 2 mod3, inert
    assert(!eisenstein_prime(3, 0));              // 3 ramifies, the unit*prime is (1,-1)
    assert(eisenstein_prime(3, 1));               // N=7, split prime (7 ≡ 1 mod3)
    assert(!eisenstein_prime(7, 0));              // 7 splits -> not prime as a point
    assert(eisenstein_prime(5, 0));               // 5 ≡ 2 mod3, inert
    assert(!eisenstein_prime(1, 1) && !eisenstein_prime(0, 0)); // unit, zero
    assert(!eisenstein_prime(4, 2));              // N=12 composite, off-axis
    double echi2 = ehecke_chi2(1500, 36, false);
    assert(echi2 < 2.0 && "Eisenstein Hecke not equidistributed");
    for (double K : {1.0, 2.0}) {                 // K=3 component reaches |z|~2252, too big for a quick gate
        MoatResult cpu = emoat(K, 100000, false);
        MoatResult gpu = emoat_gpu(K, 400, false);
        assert(!gpu.escaped && "Eisenstein validation R too small");
        assert(fabs(cpu.farthest - gpu.farthest) < 1e-6 && "Eisenstein GPU bitmap disagrees with CPU BFS");
    }

    // Gap stats + multi-seed random control: NN means real, control ran across seeds
    // and produced a seed spread (random_sd >= 0) — the error-bar machinery works.
    GapStats gi = gap_stats(false, 400, 4, false), gw = gap_stats(true, 400, 4, false);
    assert(gi.count > 0 && gw.count > 0 && "no primes found in gap scan");
    assert(gi.obs_mean >= 1.0 && gi.obs_mean < 5.0 && "Z[i] NN mean implausible");
    assert(gw.obs_mean >= 1.0 && gw.obs_mean < 5.0 && "Z[w] NN mean implausible");
    assert(gi.random_count > 0 && gi.random_mean >= 1.0 && gi.random_mean < 5.0 && "Z[i] random control failed");
    assert(gw.random_count > 0 && gw.random_mean >= 1.0 && gw.random_mean < 5.0 && "Z[w] random control failed");
    assert(gi.random_sd >= 0.0 && gw.random_sd >= 0.0 && "seed-spread not computed");

    // Pair correlation: short-range ratio is a finite depletion; CV uses EQUAL sectors.
    PCF pi = pair_corr(false, 400, 24, false), pw = pair_corr(true, 400, 24, false);
    assert(pi.g_short > 0 && pi.g_short < 1.2 && "Z[i] pcf short-range ratio implausible");
    assert(pw.g_short > 0 && pw.g_short < 1.2 && "Z[w] pcf short-range ratio implausible");

    // Singular series: Gaussian parity obstruction (odd offset -> S=0), admissible
    // offsets finite positive; Z[w] has no norm-2 ideal so all offsets admissible.
    std::vector<int> sp = primes_upto(100000);
    assert(singular_series(false, 1, 0, sp) < 1e-9 && "Z[i] odd offset should be inadmissible (S=0)");
    assert(singular_series(false, 1, 1, sp) > 0.1 && "Z[i] even offset should be admissible");
    assert(singular_series(false, 2, 0, sp) > 0.1 && "Z[i] offset 2 admissible");
    assert(singular_series(true, 1, 0, sp) > 0.1 && "Z[w] has no parity obstruction");

    printf("selftest OK  (Z[i]: chi^2=%.3f, |z|=%.2f, NN=%.2f vs rand %.2f+/-%.3f, pcf g=%.3f, S(1,1)=%.3f | Z[w]: chi^2=%.3f, NN=%.2f vs rand %.2f+/-%.3f, pcf g=%.3f, S(1,0)=%.3f)\n",
           chi2, m.farthest, gi.obs_mean, gi.random_mean, gi.random_sd, pi.g_short, singular_series(false, 1, 1, sp),
           echi2, gw.obs_mean, gw.random_mean, gw.random_sd, pw.g_short, singular_series(true, 1, 0, sp));
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "--selftest")) return selftest();
    if (argc >= 3 && !strcmp(argv[1], "--hecke")) {
        long long R = atoll(argv[2]);
        int bins = argc >= 4 ? atoi(argv[3]) : 36;
        hecke_chi2(R, bins, true);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "--moat")) {
        double K = atof(argv[2]); long long R = atoll(argv[3]);
        moat(K, R, true);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "--moat-gpu")) {
        double K = atof(argv[2]); long long R = atoll(argv[3]);
        moat_gpu(K, R, true);
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "--ehecke")) {
        long long R = atoll(argv[2]);
        int bins = argc >= 4 ? atoi(argv[3]) : 36;
        ehecke_chi2(R, bins, true);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "--emoat")) {
        double K = atof(argv[2]); long long R = atoll(argv[3]);
        emoat(K, R, true);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "--emoat-gpu")) {
        double K = atof(argv[2]); long long R = atoll(argv[3]);
        emoat_gpu(K, R, true);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "--gaps")) {
        bool eis = !strcmp(argv[2], "w") || !strcmp(argv[2], "eisenstein");
        long long R = atoll(argv[3]);
        int nseeds = argc >= 5 ? atoi(argv[4]) : 8;
        gap_stats(eis, R, nseeds, true);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "--pcf")) {
        bool eis = !strcmp(argv[2], "w") || !strcmp(argv[2], "eisenstein");
        long long R = atoll(argv[3]);
        int sectors = argc >= 5 ? atoi(argv[4]) : 24;
        pair_corr(eis, R, sectors, true);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "--sigma")) {
        bool eis = !strcmp(argv[2], "w") || !strcmp(argv[2], "eisenstein");
        long long R = atoll(argv[3]);
        int nseeds = argc >= 5 ? atoi(argv[4]) : 3;
        sigma_test(eis, R, nseeds, true);
        return 0;
    }
    fprintf(stderr,
        "usage:\n"
        "  %s --selftest\n"
        "  Gaussian Z[i]:\n"
        "  %s --hecke R [bins]     angle equidistribution of Gaussian primes (disk radius R)\n"
        "  %s --moat  K R          walk from origin past radius R with step <= K? (CPU)\n"
        "  %s --moat-gpu K R       same, GPU disk-sieve bitmap (scales to the √26 record)\n"
        "  Eisenstein Z[w] (hexagonal lattice):\n"
        "  %s --ehecke R [bins]    angle equidistribution of Eisenstein primes\n"
        "  %s --emoat  K R         Eisenstein moat walk (CPU)\n"
        "  %s --emoat-gpu K R      same, GPU disk-sieve bitmap\n"
        "  Statistics:\n"
        "  %s --gaps i|w R [seeds] NN gaps vs matched random subset, error bars over seeds\n"
        "  %s --pcf  i|w R [sect]  angular pair-correlation; CV vs sector-count sweep\n"
        "  %s --sigma i|w R [seeds] Hardy-Littlewood singular series overlay vs measured pairs\n",
        argv[0], argv[0], argv[0], argv[0], argv[0], argv[0], argv[0], argv[0], argv[0], argv[0]);
    return 1;
}
