// sieve.cuh — GPU pipeline: fused segmented sieve of Eratosthenes + octal band statistics.
//
// Bitmap convention: one bit per ODD number, bit i <-> n = 2i + 1, bit set = composite.
// n = 1 is force-marked composite; bits covering n >= N (tail of last segment) are
// force-marked composite so downstream scans need no bounds logic.
//
// Each block owns one segment of SEG_SPAN integers: it sieves the segment's odds
// in a shared-memory bitmap, then walks the segment once more to accumulate the
// per-band octal statistics into shared u32 counters, and finally flushes both
// the bitmap (to global memory, for the CPU post-pass) and the counters (via
// global u64 atomics).
#pragma once
#include <chrono>
#include <cstring>
#include <string>
#include <vector>

#include "octal_core.h"

constexpr int SEG_BITS_LOG2 = 17;
constexpr int SEG_BITS = 1 << SEG_BITS_LOG2;  // odd numbers per segment
constexpr int SEG_WORDS = SEG_BITS / 32;      // 4096 u32 = 16 KB shared
constexpr u64 SEG_SPAN = 2ULL * SEG_BITS;     // 262144 integers per segment
constexpr int SIEVE_BLOCK = 256;
constexpr u32 COOP_PRIME_LIMIT = 256;  // primes below this are marked block-cooperatively

struct PipelineResult {
    Stats stats{};            // includes host-side adjustment for the even prime 2
    std::vector<u32> bitmap;  // odds bitmap (empty unless want_bitmap)
    u64 N = 0;
    int K = 0;                // N = 8^K
    float kernel_ms = 0.0f;
    double wall_s = 0.0;
};

#if defined(__CUDACC__)

__device__ __forceinline__ u64 first_odd_multiple_geq(u32 p, u64 lo) {
    u64 m = (lo + p - 1) / p * p;
    if (!(m & 1)) m += p;  // p is odd, so m or m+p is odd
    return m;
}

__global__ void sieve_stats_kernel(const u32 *__restrict__ base_primes, int n_primes,
                                   int n_coop, u64 N, u32 *__restrict__ g_bitmap,
                                   u64 *__restrict__ g_stats) {
    __shared__ u32 s_bm[SEG_WORDS];
    __shared__ u32 s_st[STATS_FLAT];
    const u64 base = (u64)blockIdx.x * SEG_SPAN;
    const u64 seg_hi = base + SEG_SPAN;

    for (int i = threadIdx.x; i < SEG_WORDS; i += blockDim.x) s_bm[i] = 0;
    for (int i = threadIdx.x; i < STATS_FLAT; i += blockDim.x) s_st[i] = 0;
    __syncthreads();

    // Small primes hit the segment many times: the whole block cooperates on each,
    // threads taking every blockDim-th multiple. (Uniform loop, no divergence.)
    for (int j = 0; j < n_coop; ++j) {
        u32 p = base_primes[j];
        u64 pp = (u64)p * p;
        if (pp >= seg_hi) break;
        u64 m = (pp >= base) ? pp : first_odd_multiple_geq(p, base);
        u64 b0 = (m - base) >> 1;  // consecutive odd multiples of p are p bits apart
        for (u64 bit = b0 + (u64)threadIdx.x * p; bit < SEG_BITS;
             bit += (u64)blockDim.x * p)
            atomicOr(&s_bm[bit >> 5], 1u << ((int)bit & 31));
    }
    // Larger primes hit the segment a few times each: one prime per thread.
    for (int j = n_coop + threadIdx.x; j < n_primes; j += blockDim.x) {
        u32 p = base_primes[j];
        u64 pp = (u64)p * p;
        if (pp >= seg_hi) break;  // primes are sorted ascending
        u64 m = (pp >= base) ? pp : first_odd_multiple_geq(p, base);
        for (u64 bit = (m - base) >> 1; bit < SEG_BITS; bit += p)
            atomicOr(&s_bm[bit >> 5], 1u << ((int)bit & 31));
    }
    if (base == 0 && threadIdx.x == 0) atomicOr(&s_bm[0], 1u);  // n = 1 is not prime
    if (seg_hi > N) {  // force-mark the tail (n >= N) composite
        u64 tail0 = (N > base) ? ((N - base) >> 1) : 0;
        for (u64 bit = tail0 + threadIdx.x; bit < SEG_BITS; bit += blockDim.x)
            atomicOr(&s_bm[bit >> 5], 1u << ((int)bit & 31));
    }
    __syncthreads();

    // Octal statistics for every odd n in the segment.
    for (int bit = threadIdx.x; bit < SEG_BITS; bit += blockDim.x) {
        u64 n = base + 2ULL * bit + 1;
        if (n >= N) break;
        if (n == 1) continue;
        bool prime = !((s_bm[bit >> 5] >> (bit & 31)) & 1);
        int band = oct_num_digits(n);
        u32 *bs = &s_st[(band - 1) * BAND_STRIDE];
        atomicAdd(&bs[IDX_ODD_TOTAL], 1u);
        bool cand = oct_prime_candidate(n);
        if (cand) atomicAdd(&bs[IDX_CAND], 1u);
        if (prime) {
            atomicAdd(&bs[IDX_PRIMES], 1u);
            if (cand) atomicAdd(&bs[IDX_CANDP], 1u);
            atomicAdd(&bs[IDX_LAST + (int)(n & 7)], 1u);
            atomicAdd(&bs[IDX_LEAD + oct_leading_digit(n)], 1u);
            atomicAdd(&bs[IDX_M64 + (int)(n & 63)], 1u);
            u64 t = n;
            do { atomicAdd(&bs[IDX_DF + (int)(t & 7)], 1u); t >>= 3; } while (t);
        }
    }
    for (int i = threadIdx.x; i < SEG_WORDS; i += blockDim.x)
        g_bitmap[(u64)blockIdx.x * SEG_WORDS + i] = s_bm[i];
    __syncthreads();
    for (int i = threadIdx.x; i < STATS_FLAT; i += blockDim.x)
        if (s_st[i])
            atomicAdd((unsigned long long *)&g_stats[i], (unsigned long long)s_st[i]);
}

// Sieve and analyze [1, 8^K). K in [1, 12]. Returns false with a message in err on failure.
inline bool run_pipeline(int K, bool want_bitmap, PipelineResult &out, std::string &err) {
    if (K < 1 || K > 12) {
        err = "K must be in [1, 12]";
        return false;
    }
    auto t0 = std::chrono::steady_clock::now();
    const u64 N = pow8(K);
    out.N = N;
    out.K = K;
    out.bitmap.clear();
    memset(&out.stats, 0, sizeof(Stats));

    // Odd base primes up to floor(sqrt(N)) via a small CPU sieve (2 is implicit in
    // the odds-only bitmap).
    u64 S = 1;
    while ((S + 1) * (S + 1) <= N) ++S;  // N <= 2^36, no overflow
    std::vector<u32> base_primes;
    {
        std::vector<u8> comp(S + 1, 0);
        for (u64 p = 3; p * p <= S; p += 2)
            if (!comp[p])
                for (u64 m = p * p; m <= S; m += 2 * p) comp[m] = 1;
        for (u64 p = 3; p <= S; p += 2)
            if (!comp[p]) base_primes.push_back((u32)p);
    }
    int n_coop = 0;
    while (n_coop < (int)base_primes.size() && base_primes[n_coop] < COOP_PRIME_LIMIT)
        ++n_coop;

    const u64 nsegs = (N + SEG_SPAN - 1) / SEG_SPAN;
    const u64 nwords = nsegs * SEG_WORDS;

    u32 *d_primes = nullptr, *d_bitmap = nullptr;
    u64 *d_stats = nullptr;
    cudaEvent_t ev0 = nullptr, ev1 = nullptr;
    auto fail = [&](const std::string &what, cudaError_t e) {
        err = what + ": " + cudaGetErrorString(e);
        if (d_primes) cudaFree(d_primes);
        if (d_bitmap) cudaFree(d_bitmap);
        if (d_stats) cudaFree(d_stats);
        if (ev0) cudaEventDestroy(ev0);
        if (ev1) cudaEventDestroy(ev1);
        return false;
    };
    cudaError_t e;
    size_t primes_bytes = base_primes.size() * sizeof(u32);
    if ((e = cudaMalloc(&d_primes, primes_bytes ? primes_bytes : 4)) != cudaSuccess)
        return fail("cudaMalloc(base primes)", e);
    if (primes_bytes &&
        (e = cudaMemcpy(d_primes, base_primes.data(), primes_bytes,
                        cudaMemcpyHostToDevice)) != cudaSuccess)
        return fail("cudaMemcpy(base primes)", e);
    if ((e = cudaMalloc(&d_bitmap, nwords * sizeof(u32))) != cudaSuccess)
        return fail("cudaMalloc(bitmap)", e);
    if ((e = cudaMalloc(&d_stats, STATS_FLAT * sizeof(u64))) != cudaSuccess)
        return fail("cudaMalloc(stats)", e);
    if ((e = cudaMemset(d_stats, 0, STATS_FLAT * sizeof(u64))) != cudaSuccess)
        return fail("cudaMemset(stats)", e);

    cudaEventCreate(&ev0);
    cudaEventCreate(&ev1);
    cudaEventRecord(ev0);
    sieve_stats_kernel<<<(unsigned)nsegs, SIEVE_BLOCK>>>(
        d_primes, (int)base_primes.size(), n_coop, N, d_bitmap, d_stats);
    cudaEventRecord(ev1);
    if ((e = cudaDeviceSynchronize()) != cudaSuccess) return fail("kernel", e);
    if ((e = cudaGetLastError()) != cudaSuccess) return fail("kernel launch", e);
    cudaEventElapsedTime(&out.kernel_ms, ev0, ev1);

    if ((e = cudaMemcpy(out.stats.flat(), d_stats, STATS_FLAT * sizeof(u64),
                        cudaMemcpyDeviceToHost)) != cudaSuccess)
        return fail("cudaMemcpy(stats)", e);
    if (N > 2) {  // account for the even prime 2 (band 1, octal digit 2)
        u64 *b1 = out.stats.c[0];
        b1[IDX_PRIMES]++;
        b1[IDX_LAST + 2]++;
        b1[IDX_LEAD + 2]++;
        b1[IDX_DF + 2]++;
        b1[IDX_M64 + 2]++;
    }
    if (want_bitmap) {
        try {
            out.bitmap.resize(nwords);
        } catch (const std::bad_alloc &) {
            return fail("host bitmap allocation", cudaErrorMemoryAllocation);
        }
        if ((e = cudaMemcpy(out.bitmap.data(), d_bitmap, nwords * sizeof(u32),
                            cudaMemcpyDeviceToHost)) != cudaSuccess)
            return fail("cudaMemcpy(bitmap)", e);
    }
    cudaFree(d_primes);
    cudaFree(d_bitmap);
    cudaFree(d_stats);
    cudaEventDestroy(ev0);
    cudaEventDestroy(ev1);
    out.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

#else
// Translation units not compiled by nvcc only get the data structures above.
bool run_pipeline(int K, bool want_bitmap, PipelineResult &out, std::string &err);
#endif
