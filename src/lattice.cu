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

    printf("selftest OK  (classification + Hecke chi^2=%.3f + step-2 moat |z|=%.2f + GPU==CPU farthest for K=1.5,2,2.9)\n",
           chi2, m.farthest);
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
    fprintf(stderr,
        "usage:\n"
        "  %s --selftest\n"
        "  %s --hecke R [bins]     angle equidistribution of Gaussian primes (disk radius R)\n"
        "  %s --moat  K R          can primes be walked from origin past radius R with step <= K?\n"
        "  %s --moat-gpu K R       same, GPU disk-sieve bitmap (scales to the √26 record)\n",
        argv[0], argv[0], argv[0], argv[0]);
    return 1;
}
