// post.h — CPU post-pass over the odds bitmap: prime gaps, twins, octal palindromic primes.
#pragma once
#include <vector>

#include "octal_core.h"

constexpr int GAP_HIST_MAX = 1024;  // gaps below 8^12 are far smaller than this
constexpr int PALIN_EXAMPLES_MAX = 8;

struct PostStats {
    u64 gap_hist[GAP_HIST_MAX + 1] = {};  // [gap] = count; slot GAP_HIST_MAX = overflow
    u64 gap_mod8[8] = {};                 // gap & 7 (the 2->3 gap of 1 is the only odd gap)
    u64 twins_total = 0;
    u64 max_gap = 0;
    u64 max_gap_after = 0;                // lower prime of the maximal gap
    u64 band_twins[MAX_BANDS] = {};       // twin counted in band of the upper prime
    u64 band_max_gap[MAX_BANDS] = {};     // gap attributed to band of the upper prime
    u64 band_max_gap_after[MAX_BANDS] = {};
    u64 band_palin[MAX_BANDS] = {};
    u64 palin_total = 0;
    std::vector<u64> palin_examples[MAX_BANDS];  // first few per band
    u64 primes_seen = 0;                  // cross-check against Stats
};

// Scan the odds bitmap produced by run_pipeline. The prime 2 is accounted for here
// (it is not in the bitmap); tail bits (n >= N) arrive pre-marked composite.
inline void post_pass(const std::vector<u32> &bitmap, u64 N, PostStats &out) {
    out = PostStats{};
    if (N <= 2) return;
    out.primes_seen = 1;  // the prime 2: 1-digit palindrome, no predecessor gap
    out.palin_total = 1;
    out.band_palin[0] = 1;
    out.palin_examples[0].push_back(2);

    u64 prev = 2;
    int band = 1;
    u64 band_hi = 8;
    const u64 nwords = bitmap.size();
    for (u64 w = 0; w < nwords; ++w) {
        u32 inv = ~bitmap[w];  // set bits = primes
        while (inv) {
            int b = __builtin_ctz(inv);
            inv &= inv - 1;
            u64 n = 2 * ((w << 5) + b) + 1;
            while (n >= band_hi) {
                ++band;
                band_hi = (band >= MAX_BANDS) ? ~0ULL : pow8(band);
            }
            u64 gap = n - prev;
            out.gap_hist[gap < GAP_HIST_MAX ? gap : GAP_HIST_MAX]++;
            out.gap_mod8[gap & 7]++;
            if (gap == 2) {
                out.twins_total++;
                out.band_twins[band - 1]++;
            }
            if (gap > out.max_gap) {
                out.max_gap = gap;
                out.max_gap_after = prev;
            }
            if (gap > out.band_max_gap[band - 1]) {
                out.band_max_gap[band - 1] = gap;
                out.band_max_gap_after[band - 1] = prev;
            }
            if (oct_is_palindrome(n)) {
                out.palin_total++;
                out.band_palin[band - 1]++;
                if ((int)out.palin_examples[band - 1].size() < PALIN_EXAMPLES_MAX)
                    out.palin_examples[band - 1].push_back(n);
            }
            out.primes_seen++;
            prev = n;
        }
    }
}
