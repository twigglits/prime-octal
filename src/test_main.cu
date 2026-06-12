// test_main.cu — full test suite for prime-octal.
//
// Reference implementations here are deliberately independent of the production
// code paths: octal digit references go through printf("%llo") strings or
// div/mod-8 loops, primality ground truth comes from a classic CPU sieve.
//
// Build & run: make test
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <string>
#include <vector>

#include "octal_core.h"
#include "post.h"
#include "primality.h"
#include "sieve.cuh"

static long g_checks = 0, g_fails = 0;
static int g_section_fail_prints = 0;

#define SECTION(name)                            \
    do {                                         \
        printf("== %s\n", name);                 \
        fflush(stdout);                          \
        g_section_fail_prints = 0;               \
    } while (0)

#define CHECK(cond, fmt, ...)                                                  \
    do {                                                                       \
        g_checks++;                                                            \
        if (!(cond)) {                                                         \
            g_fails++;                                                         \
            if (++g_section_fail_prints <= 10)                                 \
                printf("  FAIL %s:%d  " fmt "\n", __FILE__, __LINE__,          \
                       ##__VA_ARGS__);                                         \
            if (g_section_fail_prints == 11)                                   \
                printf("  ... further failures in this section suppressed\n"); \
        }                                                                      \
    } while (0)

#define CUCHECK(call)                                                         \
    do {                                                                      \
        cudaError_t e_ = (call);                                              \
        CHECK(e_ == cudaSuccess, "CUDA error: %s", cudaGetErrorString(e_));   \
    } while (0)

// ---------------------------------------------------------------------------
// Independent reference implementations
// ---------------------------------------------------------------------------

static std::string oct_str(u64 n) {
    char b[32];
    snprintf(b, sizeof b, "%llo", n);
    return b;
}
static int ref_num_digits(u64 n) { return (int)oct_str(n).size(); }
static int ref_last_digit(u64 n) { std::string s = oct_str(n); return s[s.size() - 1] - '0'; }
static int ref_leading_digit(u64 n) { return oct_str(n)[0] - '0'; }
static int ref_digit_sum(u64 n) {
    int t = 0;
    for (char c : oct_str(n)) t += c - '0';
    return t;
}
static int ref_alt_sum(u64 n) {  // d0 - d1 + d2 - ... counted from least significant digit
    std::string s = oct_str(n);
    int L = (int)s.size(), t = 0;
    for (int i = 0; i < L; ++i) {
        int d = s[L - 1 - i] - '0';
        t += (i % 2) ? -d : d;
    }
    return t;
}
static int ref_weighted_mod5(u64 n) {  // weights 8^i mod 5 = 1,3,4,2 cycling from LSD
    static const int w[4] = {1, 3, 4, 2};
    std::string s = oct_str(n);
    int L = (int)s.size(), t = 0;
    for (int i = 0; i < L; ++i) t += (s[L - 1 - i] - '0') * w[i % 4];
    return t % 5;
}
static bool ref_palindrome(u64 n) {
    std::string s = oct_str(n), r = s;
    std::reverse(r.begin(), r.end());
    return s == r;
}
// Digits via div/mod (second independent path, used by the bulk stats reference).
static int ref_digits_divmod(u64 n, int d[MAX_BANDS]) {
    int L = 0;
    do { d[L++] = (int)(n % 8); n /= 8; } while (n);
    return L;
}

static std::vector<u8> ref_sieve(u64 N) {  // comp[n] = 1 iff n composite (or 0/1)
    std::vector<u8> comp(N, 0);
    if (N > 0) comp[0] = 1;
    if (N > 1) comp[1] = 1;
    for (u64 p = 2; p * p < N; ++p)
        if (!comp[p])
            for (u64 m = p * p; m < N; m += p) comp[m] = 1;
    return comp;
}

// Mirrors the pipeline's statistics semantics exactly (see sieve.cuh):
// odd n only, n=1 skipped, then host-side adjustment for the prime 2.
static void ref_stats(u64 N, const std::vector<u8> &comp, Stats &s) {
    memset(&s, 0, sizeof s);
    for (u64 n = 3; n < N; n += 2) {
        int dg[MAX_BANDS];
        int L = ref_digits_divmod(n, dg);
        u64 *bs = s.c[L - 1];
        bs[IDX_ODD_TOTAL]++;
        int alt = 0, ds = 0, w5 = 0;
        static const int w[4] = {1, 3, 4, 2};
        for (int i = 0; i < L; ++i) {
            ds += dg[i];
            alt += (i % 2) ? -dg[i] : dg[i];
            w5 += dg[i] * w[i % 4];
        }
        bool cand = (dg[0] % 2 == 1) && (((alt % 3) + 3) % 3 != 0) && (ds % 7 != 0) &&
                    (w5 % 5 != 0);
        bool prime = !comp[n];
        if (cand) bs[IDX_CAND]++;
        if (prime) {
            bs[IDX_PRIMES]++;
            if (cand) bs[IDX_CANDP]++;
            bs[IDX_LAST + dg[0]]++;
            bs[IDX_LEAD + dg[L - 1]]++;
            bs[IDX_M64 + (n & 63)]++;
            for (int i = 0; i < L; ++i) bs[IDX_DF + dg[i]]++;
        }
    }
    if (N > 2) {  // the even prime, same adjustment the pipeline applies host-side
        u64 *b1 = s.c[0];
        b1[IDX_PRIMES]++;
        b1[IDX_LAST + 2]++;
        b1[IDX_LEAD + 2]++;
        b1[IDX_DF + 2]++;
        b1[IDX_M64 + 2]++;
    }
}

// Mirrors post_pass semantics from the full prime list.
static void ref_post(u64 N, const std::vector<u8> &comp, PostStats &out) {
    std::vector<u64> primes;
    for (u64 n = 2; n < N; ++n)
        if (!comp[n]) primes.push_back(n);
    u64 prev = 0;
    for (u64 p : primes) {
        int dg[MAX_BANDS];
        int L = ref_digits_divmod(p, dg);
        if (prev) {
            u64 gap = p - prev;
            out.gap_hist[gap < GAP_HIST_MAX ? gap : GAP_HIST_MAX]++;
            out.gap_mod8[gap & 7]++;
            if (gap == 2) { out.twins_total++; out.band_twins[L - 1]++; }
            if (gap > out.max_gap) { out.max_gap = gap; out.max_gap_after = prev; }
            if (gap > out.band_max_gap[L - 1]) {
                out.band_max_gap[L - 1] = gap;
                out.band_max_gap_after[L - 1] = prev;
            }
        }
        bool pal = true;
        for (int i = 0, j = L - 1; i < j; ++i, --j)
            if (dg[i] != dg[j]) { pal = false; break; }
        if (pal) {
            out.palin_total++;
            out.band_palin[L - 1]++;
            if ((int)out.palin_examples[L - 1].size() < PALIN_EXAMPLES_MAX)
                out.palin_examples[L - 1].push_back(p);
        }
        out.primes_seen++;
        prev = p;
    }
}

static std::vector<u64> feature_test_values() {
    std::vector<u64> v;
    for (u64 n = 0; n < (1ULL << 19); ++n) v.push_back(n);
    const u64 extra[] = {
        1ULL << 19, (1ULL << 20) + 12345, 0xFFFFFFFFULL, 1ULL << 32, (1ULL << 33) + 7,
        0x123456789ABCDEFULL, 1234567890123456789ULL,
        pow8(10) - 1, pow8(10), pow8(10) + 1,
        pow8(12) - 1, pow8(12), pow8(12) + 1,
        pow8(21) - 1, pow8(21), pow8(21) + 1,
        1ULL << 63, ~0ULL,
    };
    for (u64 n : extra) v.push_back(n);
    return v;
}

// ---------------------------------------------------------------------------
// 1. Octal digit primitives vs string references
// ---------------------------------------------------------------------------
static void test_octal_features() {
    SECTION("octal digit primitives vs string-based references");
    std::vector<u64> vals = feature_test_values();
    for (u64 n : vals) {
        CHECK(oct_num_digits(n) == ref_num_digits(n), "num_digits(%llu): got %d want %d", n,
              oct_num_digits(n), ref_num_digits(n));
        CHECK(oct_last_digit(n) == ref_last_digit(n), "last_digit(%llu): got %d want %d", n,
              oct_last_digit(n), ref_last_digit(n));
        CHECK(oct_leading_digit(n) == ref_leading_digit(n),
              "leading_digit(%llu): got %d want %d", n, oct_leading_digit(n),
              ref_leading_digit(n));
        CHECK(oct_digit_sum(n) == ref_digit_sum(n), "digit_sum(%llu): got %d want %d", n,
              oct_digit_sum(n), ref_digit_sum(n));
        CHECK(oct_alt_sum(n) == ref_alt_sum(n), "alt_sum(%llu): got %d want %d", n,
              oct_alt_sum(n), ref_alt_sum(n));
        CHECK(oct_weighted_mod5(n) == ref_weighted_mod5(n),
              "weighted_mod5(%llu): got %d want %d", n, oct_weighted_mod5(n),
              ref_weighted_mod5(n));
        CHECK(oct_is_palindrome(n) == ref_palindrome(n), "palindrome(%llu): got %d want %d",
              n, (int)oct_is_palindrome(n), (int)ref_palindrome(n));
    }
}

// ---------------------------------------------------------------------------
// 2. Base-8 divisibility identities (the math the predictor rests on)
// ---------------------------------------------------------------------------
static void test_divisibility_identities() {
    SECTION("base-8 divisibility identities");
    for (u64 n = 1; n < 300000; ++n) {
        CHECK((u64)(oct_last_digit(n) & 1) == (n & 1), "parity via last digit, n=%llu", n);
        int a = oct_alt_sum(n);
        CHECK((u64)(((a % 3) + 3) % 3) == n % 3, "mod 3 via alternating sum, n=%llu", n);
        CHECK((u64)(oct_digit_sum(n) % 7) == n % 7, "mod 7 via digit sum, n=%llu", n);
        CHECK((u64)oct_weighted_mod5(n) == n % 5, "mod 5 via weighted sum, n=%llu", n);
    }
}

// ---------------------------------------------------------------------------
// 3. Octal candidate rule == coprimality to 210
// ---------------------------------------------------------------------------
static void test_candidate_rule() {
    SECTION("octal candidate rule == gcd(n,210)==1");
    for (u64 n = 1; n < 300000; ++n) {
        bool want = std::gcd(n, (u64)210) == 1;
        CHECK(oct_prime_candidate(n) == want, "candidate(%llu): got %d want %d", n,
              (int)oct_prime_candidate(n), (int)want);
    }
    // The four primes 2,3,5,7 each fail their own octal rule by construction.
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL})
        CHECK(!oct_prime_candidate(p), "candidate(%llu) should be false", p);
}

// ---------------------------------------------------------------------------
// 4. Palindrome facts in octal
// ---------------------------------------------------------------------------
static void test_palindrome_facts() {
    SECTION("octal palindrome facts");
    std::vector<u8> comp = ref_sieve(1 << 16);
    // No 2-digit octal palindromic prime exists: "dd" in base 8 is 9d, divisible by 3.
    for (int d = 1; d <= 7; ++d) {
        u64 n = (u64)d * 9;
        CHECK(oct_is_palindrome(n), "0o%d%d should be octal palindrome", d, d);
        CHECK(comp[n] == 1, "9*%d should be composite", d);
    }
    CHECK(oct_is_palindrome(73) && !comp[73], "73 = 0o111 is a palindromic prime");
    CHECK(oct_is_palindrome(89) && !comp[89], "89 = 0o131 is a palindromic prime");
}

// ---------------------------------------------------------------------------
// 5. Device/host equivalence of the octal primitives
// ---------------------------------------------------------------------------
__global__ void feat_kernel(const u64 *in, int cnt, int *nd, int *last, int *lead, int *ds,
                            int *alt, int *w5, int *pal, int *cand) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= cnt) return;
    u64 n = in[i];
    nd[i] = oct_num_digits(n);
    last[i] = oct_last_digit(n);
    lead[i] = oct_leading_digit(n);
    ds[i] = oct_digit_sum(n);
    alt[i] = oct_alt_sum(n);
    w5[i] = oct_weighted_mod5(n);
    pal[i] = oct_is_palindrome(n) ? 1 : 0;
    cand[i] = oct_prime_candidate(n) ? 1 : 0;
}

static void test_device_equivalence() {
    SECTION("device/host equivalence of octal primitives");
    std::vector<u64> vals;
    for (u64 n = 1; n < (1ULL << 20); ++n) vals.push_back(n);
    for (u64 n : feature_test_values()) vals.push_back(n);
    int cnt = (int)vals.size();

    u64 *d_in = nullptr;
    int *d_out[8] = {};
    CUCHECK(cudaMalloc(&d_in, cnt * sizeof(u64)));
    for (int j = 0; j < 8; ++j) CUCHECK(cudaMalloc(&d_out[j], cnt * sizeof(int)));
    CUCHECK(cudaMemcpy(d_in, vals.data(), cnt * sizeof(u64), cudaMemcpyHostToDevice));
    feat_kernel<<<(cnt + 255) / 256, 256>>>(d_in, cnt, d_out[0], d_out[1], d_out[2],
                                            d_out[3], d_out[4], d_out[5], d_out[6],
                                            d_out[7]);
    CUCHECK(cudaGetLastError());
    std::vector<int> h(8ULL * cnt);
    for (int j = 0; j < 8; ++j)
        CUCHECK(cudaMemcpy(h.data() + (u64)j * cnt, d_out[j], cnt * sizeof(int),
                           cudaMemcpyDeviceToHost));
    for (int i = 0; i < cnt; ++i) {
        u64 n = vals[i];
        CHECK(h[i] == oct_num_digits(n), "gpu num_digits(%llu)", n);
        CHECK(h[cnt + i] == oct_last_digit(n), "gpu last_digit(%llu)", n);
        CHECK(h[2ULL * cnt + i] == oct_leading_digit(n), "gpu leading_digit(%llu)", n);
        CHECK(h[3ULL * cnt + i] == oct_digit_sum(n), "gpu digit_sum(%llu)", n);
        CHECK(h[4ULL * cnt + i] == oct_alt_sum(n), "gpu alt_sum(%llu)", n);
        CHECK(h[5ULL * cnt + i] == oct_weighted_mod5(n), "gpu weighted_mod5(%llu)", n);
        CHECK(h[6ULL * cnt + i] == (int)oct_is_palindrome(n), "gpu palindrome(%llu)", n);
        CHECK(h[7ULL * cnt + i] == (int)oct_prime_candidate(n), "gpu candidate(%llu)", n);
    }
    cudaFree(d_in);
    for (int j = 0; j < 8; ++j) cudaFree(d_out[j]);
}

// ---------------------------------------------------------------------------
// 6. Miller-Rabin ground truth
// ---------------------------------------------------------------------------
static void test_miller_rabin() {
    SECTION("deterministic Miller-Rabin for u64");
    std::vector<u8> comp = ref_sieve(1 << 20);
    for (u64 n = 0; n < (1ULL << 20); ++n)
        CHECK(is_prime_u64(n) == (n >= 2 && !comp[n]), "is_prime_u64(%llu)", n);
    const u64 known_primes[] = {
        2ULL, 3ULL, 999999937ULL, 1000003ULL, 2147483647ULL,            // 2^31 - 1
        2305843009213693951ULL,                                         // 2^61 - 1
        18446744073709551557ULL,                                        // largest u64 prime
    };
    for (u64 p : known_primes) CHECK(is_prime_u64(p), "%llu is prime", p);
    const u64 known_composites[] = {
        0ULL, 1ULL, 561ULL, 1105ULL, 1729ULL, 2465ULL, 2821ULL, 6601ULL, 8911ULL,
        25326001ULL,               // strong pseudoprime to bases 2,3,5
        3825123056546413051ULL,    // strong pseudoprime to all bases <= 23
        1000036000099ULL,          // 1000003 * 1000033
        ~0ULL, 1ULL << 63,
    };
    for (u64 c : known_composites) CHECK(!is_prime_u64(c), "%llu is composite", c);
}

// ---------------------------------------------------------------------------
// 7+8. GPU pipeline vs CPU reference, known pi(8^k), structural invariants
// ---------------------------------------------------------------------------
static void check_pipeline_invariants(const PipelineResult &res) {
    const Stats &s = res.stats;
    for (int b = 1; b <= MAX_BANDS; ++b) {
        u64 primes = s.get(b, IDX_PRIMES);
        if (b > res.K) {
            CHECK(primes == 0, "band %d beyond K=%d must be empty", b, res.K);
            continue;
        }
        u64 sl = 0, sd = 0, sm = 0;
        for (int d = 0; d < 8; ++d) sl += s.get(b, IDX_LAST + d);
        for (int d = 0; d < 8; ++d) sd += s.get(b, IDX_LEAD + d);
        for (int r = 0; r < 64; ++r) sm += s.get(b, IDX_M64 + r);
        CHECK(sl == primes, "K=%d band %d: last-digit sum %llu != primes %llu", res.K, b, sl,
              primes);
        CHECK(sd == primes, "K=%d band %d: lead-digit sum %llu != primes %llu", res.K, b, sd,
              primes);
        CHECK(sm == primes, "K=%d band %d: mod64 sum %llu != primes %llu", res.K, b, sm,
              primes);
        u64 want_odd = (b == 1) ? 3 : 7 * pow8(b - 1) / 2;
        CHECK(s.get(b, IDX_ODD_TOTAL) == want_odd, "K=%d band %d: odd_total %llu != %llu",
              res.K, b, s.get(b, IDX_ODD_TOTAL), want_odd);
        if (b == 1)
            CHECK(s.get(b, IDX_CANDP) == 0, "band 1: no candidate is prime (2,3,5,7 all fail)");
        else
            CHECK(s.get(b, IDX_CANDP) == primes,
                  "K=%d band %d: octal rules must have perfect recall (%llu != %llu)", res.K,
                  b, s.get(b, IDX_CANDP), primes);
    }
}

static void test_pipeline() {
    SECTION("GPU sieve+stats pipeline vs CPU reference (K=1..7)");
    for (int K = 1; K <= 7; ++K) {
        PipelineResult res;
        std::string err;
        bool ok = run_pipeline(K, true, res, err);
        CHECK(ok, "run_pipeline(K=%d) failed: %s", K, err.c_str());
        if (!ok) continue;
        u64 N = pow8(K);
        CHECK(res.N == N, "K=%d: N mismatch", K);
        std::vector<u8> comp = ref_sieve(N);

        // Bitmap: bit i <-> n = 2i+1; set iff composite, n==1, or tail (n >= N).
        u64 nbits = (u64)res.bitmap.size() * 32;
        CHECK(nbits * 2 >= N, "K=%d: bitmap too small", K);
        u64 bad = 0;
        for (u64 bit = 0; bit < nbits; ++bit) {
            u64 n = 2 * bit + 1;
            bool set = (res.bitmap[bit >> 5] >> (bit & 31)) & 1;
            bool want = (n >= N) || comp[n];
            if (set != want) bad++;
        }
        CHECK(bad == 0, "K=%d: %llu bitmap bits disagree with reference sieve", K, bad);

        // Full statistics comparison, every counter of every band.
        Stats want;
        ref_stats(N, comp, want);
        u64 diffs = 0;
        for (int i = 0; i < STATS_FLAT; ++i)
            if (res.stats.flat()[i] != want.flat()[i]) {
                diffs++;
                if (diffs <= 5)
                    printf("  stats[band %d, idx %d]: got %llu want %llu\n",
                           i / BAND_STRIDE + 1, i % BAND_STRIDE, res.stats.flat()[i],
                           want.flat()[i]);
            }
        CHECK(diffs == 0, "K=%d: %llu stats counters disagree with reference", K, diffs);

        for (int k = 1; k <= K; ++k)
            CHECK(res.stats.pi_upto_band(k) == PI_8_POW[k],
                  "K=%d: pi(8^%d) = %llu, known value %llu", K, k,
                  res.stats.pi_upto_band(k), PI_8_POW[k]);
        check_pipeline_invariants(res);
    }

    SECTION("GPU pipeline at K=8 (16.7M) vs reference and known pi");
    PipelineResult res;
    std::string err;
    bool ok = run_pipeline(8, true, res, err);
    CHECK(ok, "run_pipeline(K=8) failed: %s", err.c_str());
    if (ok) {
        u64 N = pow8(8);
        std::vector<u8> comp = ref_sieve(N);
        Stats want;
        ref_stats(N, comp, want);
        u64 diffs = 0;
        for (int i = 0; i < STATS_FLAT; ++i)
            if (res.stats.flat()[i] != want.flat()[i]) diffs++;
        CHECK(diffs == 0, "K=8: %llu stats counters disagree with reference", diffs);
        CHECK(res.stats.pi_upto_band(8) == PI_8_POW[8], "pi(8^8) = %llu, known %llu",
              res.stats.pi_upto_band(8), PI_8_POW[8]);
        check_pipeline_invariants(res);
    }
}

// ---------------------------------------------------------------------------
// 9. CPU post-pass (gaps, twins, palindromes) vs reference
// ---------------------------------------------------------------------------
static void test_post_pass() {
    SECTION("post-pass gaps/twins/palindromes vs reference (K=6,7)");
    for (int K = 6; K <= 7; ++K) {
        PipelineResult res;
        std::string err;
        bool ok = run_pipeline(K, true, res, err);
        CHECK(ok, "run_pipeline(K=%d) failed: %s", K, err.c_str());
        if (!ok) continue;
        u64 N = pow8(K);
        std::vector<u8> comp = ref_sieve(N);
        PostStats got, want;
        post_pass(res.bitmap, N, got);
        ref_post(N, comp, want);

        CHECK(got.primes_seen == want.primes_seen, "K=%d: primes_seen %llu want %llu", K,
              got.primes_seen, want.primes_seen);
        CHECK(got.primes_seen == res.stats.pi_upto_band(K),
              "K=%d: post-pass prime count vs pipeline stats", K);
        u64 gh = 0;
        for (int g = 0; g <= GAP_HIST_MAX; ++g)
            if (got.gap_hist[g] != want.gap_hist[g]) gh++;
        CHECK(gh == 0, "K=%d: %llu gap histogram slots disagree", K, gh);
        CHECK(want.gap_hist[GAP_HIST_MAX] == 0, "K=%d: no gap should overflow histogram", K);
        for (int r = 0; r < 8; ++r)
            CHECK(got.gap_mod8[r] == want.gap_mod8[r], "K=%d: gap_mod8[%d] %llu want %llu", K,
                  r, got.gap_mod8[r], want.gap_mod8[r]);
        CHECK(got.twins_total == want.twins_total, "K=%d: twins %llu want %llu", K,
              got.twins_total, want.twins_total);
        CHECK(got.max_gap == want.max_gap && got.max_gap_after == want.max_gap_after,
              "K=%d: max gap %llu after %llu, want %llu after %llu", K, got.max_gap,
              got.max_gap_after, want.max_gap, want.max_gap_after);
        CHECK(got.palin_total == want.palin_total, "K=%d: palindromic primes %llu want %llu",
              K, got.palin_total, want.palin_total);
        for (int b = 0; b < MAX_BANDS; ++b) {
            CHECK(got.band_twins[b] == want.band_twins[b], "K=%d band %d twins", K, b + 1);
            CHECK(got.band_max_gap[b] == want.band_max_gap[b] &&
                      got.band_max_gap_after[b] == want.band_max_gap_after[b],
                  "K=%d band %d max gap", K, b + 1);
            CHECK(got.band_palin[b] == want.band_palin[b], "K=%d band %d palindromes", K,
                  b + 1);
            CHECK(got.palin_examples[b] == want.palin_examples[b],
                  "K=%d band %d palindrome examples", K, b + 1);
        }
    }
}

int main() {
    test_octal_features();
    test_divisibility_identities();
    test_candidate_rule();
    test_palindrome_facts();
    test_device_equivalence();
    test_miller_rabin();
    test_pipeline();
    test_post_pass();
    printf("\n%s — %ld checks, %ld failures\n", g_fails ? "RED" : "GREEN", g_checks,
           g_fails);
    return g_fails ? 1 : 0;
}
