// octal_core.h — base-8 digit primitives shared by host code, tests, and CUDA kernels.
//
// Everything here treats an integer strictly through its octal (base-8) digits.
// The divisibility rules are exact identities of base 8:
//   8 ≡ 1 (mod 7)  → n ≡ octal digit sum (mod 7)
//   8 ≡ −1 (mod 3) → n ≡ alternating octal digit sum (mod 3)
//   8^i mod 5 cycles 1,3,4,2 → n ≡ weighted octal digit sum (mod 5)
//   last octal digit even ⟺ n even
#pragma once
#include <cstdint>

typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned char u8;

#if defined(__CUDACC__)
#define HD __host__ __device__ __forceinline__
#else
#define HD inline
#endif

// A u64 has at most 22 octal digits (ceil(64/3)).
constexpr int MAX_BANDS = 22;

HD u64 pow8(int k) { return 1ULL << (3 * k); }  // valid for k <= 21

// Known prime counts pi(8^k) for k = 1..12 (standard pi(2^n) tables, n = 3k).
constexpr u64 PI_8_POW[13] = {
    0ULL,           // pi(1)
    4ULL,           // pi(8)
    18ULL,          // pi(64)
    97ULL,          // pi(512)
    564ULL,         // pi(4096)
    3512ULL,        // pi(32768)
    23000ULL,       // pi(262144)
    155611ULL,      // pi(2097152)
    1077871ULL,     // pi(16777216)
    7603553ULL,     // pi(134217728)
    54400028ULL,    // pi(1073741824)
    393615806ULL,   // pi(8589934592)
    2874398515ULL,  // pi(68719476736)
};

// ---------------------------------------------------------------------------
// Flat per-band statistics layout (band = number of octal digits, 1-based).
// One u64 counter array per band; identical layout on host and device so the
// GPU can atomically accumulate into it and tests can compare it field-for-field.
// ---------------------------------------------------------------------------
enum : int {
    IDX_ODD_TOTAL = 0,  // odd n counted in this band (n = 1 excluded)
    IDX_PRIMES    = 1,  // primes in this band
    IDX_CAND      = 2,  // octal-rule candidates in this band (odd n only)
    IDX_CANDP     = 3,  // candidates that are actually prime
    IDX_LAST      = 4,  // +8: primes by last octal digit
    IDX_LEAD      = 12, // +8: primes by leading octal digit
    IDX_DF        = 20, // +8: octal digit frequency over all digits of all primes
    IDX_M64       = 28, // +64: primes by last two octal digits (n mod 64)
    BAND_STRIDE   = 92,
};
constexpr int STATS_FLAT = MAX_BANDS * BAND_STRIDE;

struct Stats {
    u64 c[MAX_BANDS][BAND_STRIDE];
    u64 *flat() { return &c[0][0]; }
    const u64 *flat() const { return &c[0][0]; }
    u64 get(int band, int idx) const { return c[band - 1][idx]; }  // band is 1-based
    u64 pi_upto_band(int k) const {
        u64 s = 0;
        for (int b = 1; b <= k && b <= MAX_BANDS; ++b) s += get(b, IDX_PRIMES);
        return s;
    }
};

// ---------------------------------------------------------------------------
// Octal digit primitives
// ---------------------------------------------------------------------------

HD int oct_msb(u64 n) {  // position of most significant set bit; n must be > 0
#ifdef __CUDA_ARCH__
    return 63 - __clzll((long long)n);
#else
    return 63 - __builtin_clzll(n);
#endif
}

HD int oct_num_digits(u64 n) { return n ? oct_msb(n) / 3 + 1 : 1; }

HD int oct_last_digit(u64 n) { return (int)(n & 7); }

HD int oct_leading_digit(u64 n) { return (int)(n >> 3 * (oct_num_digits(n) - 1)); }

HD int oct_digit_sum(u64 n) {
    int s = 0;
    do { s += (int)(n & 7); n >>= 3; } while (n);
    return s;
}

HD int oct_alt_sum(u64 n) {  // d0 - d1 + d2 - ... from the least significant digit
    int a = 0, sign = 1;
    do { a += sign * (int)(n & 7); sign = -sign; n >>= 3; } while (n);
    return a;
}

HD int oct_weighted_mod5(u64 n) {  // weights 8^i mod 5 = 1,3,4,2 cycling from LSD
    const int w[4] = {1, 3, 4, 2};
    int s = 0, i = 0;
    do { s += (int)(n & 7) * w[i & 3]; ++i; n >>= 3; } while (n);
    return s % 5;  // s <= 22 digits * 7 * 4 = 616, no overflow
}

HD bool oct_is_palindrome(u64 n) {
    int d[MAX_BANDS], L = 0;
    do { d[L++] = (int)(n & 7); n >>= 3; } while (n);
    for (int i = 0, j = L - 1; i < j; ++i, --j)
        if (d[i] != d[j]) return false;
    return true;
}

// The "octal predictor": n survives iff its octal digits pass all four
// base-8 divisibility tests (2 via last digit, 3 via alternating sum,
// 7 via digit sum, 5 via weighted sum). Equivalent to gcd(n, 210) == 1.
HD bool oct_prime_candidate(u64 n) {
    if (!(n & 1)) return false;                          // last octal digit even
    if (((oct_alt_sum(n) % 3) + 3) % 3 == 0) return false;  // divisible by 3
    if (oct_digit_sum(n) % 7 == 0) return false;         // divisible by 7
    if (oct_weighted_mod5(n) == 0) return false;         // divisible by 5
    return true;
}
