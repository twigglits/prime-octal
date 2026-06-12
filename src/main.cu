// main.cu — prime-octal CLI: GPU survey of prime emergence in base 8.
//
//   prime_octal [--octal-digits K] [--out DIR] [--no-postpass]
//   prime_octal --judge OCTAL          (predict a single number from its octal digits)
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#include "octal_core.h"
#include "post.h"
#include "primality.h"
#include "sieve.cuh"

static std::string oct_s(u64 n) {
    char b[32];
    snprintf(b, sizeof b, "%llo", n);
    return b;
}

// Offset logarithmic integral Li(x) = li(x) - li(2), the PNT estimate of pi(x).
static double Li(double x) {
    if (x < 2) return 0.0;
    auto li = [](double v) {
        double L = std::log(v), term = 1.0, s = 0.0;
        for (int k = 1; k < 300; ++k) {
            term *= L / k;
            double add = term / k;
            s += add;
            if (add < 1e-14 * s && k > (int)L) break;
        }
        return 0.5772156649015328606 + std::log(L) + s;
    };
    return li(x) - li(2.0);
}

// Writes every line to stdout and to the report file.
struct Tee {
    FILE *f = nullptr;
    void p(const char *fmt, ...) {
        va_list a;
        va_start(a, fmt);
        va_list b;
        va_copy(b, a);
        vprintf(fmt, a);
        if (f) vfprintf(f, fmt, b);
        va_end(a);
        va_end(b);
    }
};

// ---------------------------------------------------------------------------
// --judge: the octal-rule predictor applied to one number, with ground truth
// ---------------------------------------------------------------------------
static int judge(const std::string &arg) {
    std::string s = arg;
    if (s.size() > 2 && s[0] == '0' && (s[1] == 'o' || s[1] == 'O')) s = s.substr(2);
    if (s.empty()) {
        fprintf(stderr, "judge: empty octal numeral\n");
        return 1;
    }
    u64 n = 0;
    for (char c : s) {
        if (c < '0' || c > '7') {
            fprintf(stderr, "judge: '%c' is not an octal digit (base 8 uses 0-7 only)\n", c);
            return 1;
        }
        if (n > (~0ULL - (u64)(c - '0')) / 8) {
            fprintf(stderr, "judge: %s does not fit in 64 bits\n", arg.c_str());
            return 1;
        }
        n = n * 8 + (u64)(c - '0');
    }

    printf("octal   %s\n", oct_s(n).c_str());
    printf("decimal %llu\n", n);
    printf("octal digits: %d (leading %d, last %d)\n\n", oct_num_digits(n),
           oct_leading_digit(n), oct_last_digit(n));

    int last = oct_last_digit(n);
    int alt = oct_alt_sum(n);
    int alt3 = ((alt % 3) + 3) % 3;
    int ds = oct_digit_sum(n);
    int w5 = oct_weighted_mod5(n);
    printf("octal digit rules (each an exact base-8 divisibility identity):\n");
    printf("  vs 2: last digit %d           -> %s\n", last,
           (last & 1) ? "odd, passes" : "even => divisible by 2");
    printf("  vs 3: alternating digit sum %+d ≡ %d (mod 3)  -> %s\n", alt, alt3,
           alt3 ? "passes" : "divisible by 3");
    printf("  vs 7: digit sum %d ≡ %d (mod 7)       -> %s\n", ds, ds % 7,
           (ds % 7) ? "passes" : "divisible by 7");
    printf("  vs 5: weighted sum (w=1,3,4,2) ≡ %d (mod 5) -> %s\n", w5,
           w5 ? "passes" : "divisible by 5");

    bool cand = oct_prime_candidate(n);
    bool prime = is_prime_u64(n);
    printf("\nprediction from octal digits: %s\n",
           cand ? "CANDIDATE — survives all four octal digit rules"
                : "NOT PRIME — an octal digit rule eliminates it");
    if (cand && n >= 4096) {  // the asymptotic density is meaningless for tiny n
        double density = (210.0 / 48.0) / std::log((double)n);
        printf("  (about %.1f%% of octal-rule survivors of this magnitude are prime)\n",
               100.0 * density);
    }
    printf("ground truth (Miller-Rabin):  %s\n", prime ? "PRIME" : "not prime");
    if (!cand && prime)
        printf("note: %llu is one of the four rule primes 2,3,5,7 — each is rejected by\n"
               "the very divisibility rule it generates.\n",
               n);
    return 0;
}

// ---------------------------------------------------------------------------
// Report + CSV emission
// ---------------------------------------------------------------------------
static FILE *open_csv(const std::string &dir, const char *name, const char *header) {
    std::string path = dir + "/" + name;
    FILE *f = fopen(path.c_str(), "w");
    if (f) fprintf(f, "%s\n", header);
    return f;
}

int main(int argc, char **argv) {
    int K = 10;
    std::string out_dir = "results";
    bool postpass = true;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "--octal-digits" || a == "-k") && i + 1 < argc) {
            K = atoi(argv[++i]);
        } else if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (a == "--no-postpass") {
            postpass = false;
        } else if (a == "--judge" && i + 1 < argc) {
            return judge(argv[++i]);
        } else if (a == "--help" || a == "-h") {
            printf("usage: %s [--octal-digits K] [--out DIR] [--no-postpass]\n"
                   "       %s --judge OCTAL\n\n"
                   "Surveys all n < 8^K (K in [1,12], default 10) on the GPU and reports\n"
                   "prime-emergence patterns in the octal number system.\n",
                   argv[0], argv[0]);
            return 0;
        } else {
            fprintf(stderr, "unknown argument: %s (try --help)\n", argv[i]);
            return 1;
        }
    }

    PipelineResult res;
    std::string err;
    if (!run_pipeline(K, postpass, res, err)) {
        fprintf(stderr, "error: %s\n", err.c_str());
        return 1;
    }
    PostStats post;
    if (postpass) post_pass(res.bitmap, res.N, post);

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    std::string report_path = out_dir + "/report-8^" + std::to_string(K) + ".txt";
    Tee t;
    t.f = fopen(report_path.c_str(), "w");

    cudaDeviceProp prop;
    memset(&prop, 0, sizeof prop);
    cudaGetDeviceProperties(&prop, 0);

    const Stats &S = res.stats;
    u64 PI = S.pi_upto_band(K);

    t.p("prime-octal — prime numbers through the base-8 lens\n");
    t.p("====================================================\n");
    t.p("N = 8^%d = %llu (octal %s)\n", K, res.N, oct_s(res.N).c_str());
    t.p("GPU: %s | sieve+stats kernel %.1f ms | pipeline %.2f s\n", prop.name,
        res.kernel_ms, res.wall_s);
    t.p("pi(8^%d) = %llu\n\n", K, PI);

    // -- 1. emergence per band --------------------------------------------------
    t.p("1. PRIME EMERGENCE PER OCTAL BAND   (band b = numbers with b octal digits)\n");
    t.p("%4s %26s %14s %14s %8s %14s %s\n", "b", "range (octal)", "primes",
        "Li estimate", "ratio", "pi(8^b)", "known");
    for (int b = 1; b <= K; ++b) {
        u64 lo = (b == 1) ? 1 : pow8(b - 1), hi = pow8(b);
        double expect = Li((double)hi) - Li((double)lo);
        u64 primes = S.get(b, IDX_PRIMES);
        u64 cum = S.pi_upto_band(b);
        const char *known = (b <= 12) ? (cum == PI_8_POW[b] ? "✓" : "✗ MISMATCH") : "-";
        char range[40];
        snprintf(range, sizeof range, "%s..%s", oct_s(lo).c_str(), oct_s(hi - 1).c_str());
        t.p("%4d %26s %14llu %14.0f %8.4f %14llu   %s\n", b, range, primes, expect,
            expect > 0 ? primes / expect : 0.0, cum, known);
    }
    t.p("\n");

    // -- 2. the octal predictor ---------------------------------------------------
    t.p("2. OCTAL PREDICTOR — candidates from base-8 digit rules only\n");
    t.p("   rules: last digit odd | alt digit sum ≢ 0 (mod 3) | digit sum ≢ 0 (mod 7)\n");
    t.p("          | weighted digit sum (1,3,4,2) ≢ 0 (mod 5)\n");
    t.p("%4s %16s %16s %8s %14s %10s %10s %8s\n", "b", "odd numbers", "candidates",
        "cand%", "primes", "precision", "theory", "recall");
    for (int b = 1; b <= K; ++b) {
        u64 odd = S.get(b, IDX_ODD_TOTAL), cand = S.get(b, IDX_CAND);
        u64 primes = S.get(b, IDX_PRIMES), candp = S.get(b, IDX_CANDP);
        u64 lo = (b == 1) ? 1 : pow8(b - 1), hi = pow8(b);
        double expect = Li((double)hi) - Li((double)lo);
        t.p("%4d %16llu %16llu %7.2f%% %14llu ", b, odd, cand,
            odd ? 100.0 * cand / odd : 0.0, primes);
        if (cand)
            t.p("%9.3f%% %9.3f%% ", 100.0 * candp / cand, 100.0 * expect / cand);
        else
            t.p("%10s %10s ", "-", "-");
        t.p("%7.2f%%\n", primes ? 100.0 * candp / primes : 0.0);
    }
    t.p("   every prime > 7 survives the rules (recall 100%% beyond band 1); the four\n");
    t.p("   rule primes 2,3,5,7 are each rejected by their own rule. The rules keep\n");
    t.p("   48/210 = 22.857%% of integers, a 4.375x concentration of primes.\n\n");

    // -- 3. last octal digit ------------------------------------------------------
    t.p("3. LAST OCTAL DIGIT OF PRIMES (p mod 8) — per band\n");
    t.p("   every odd square ends in octal 1, so digit 1 is the quadratic-residue\n");
    t.p("   class: Chebyshev bias predicts it trails its peers.\n");
    t.p("%4s %14s %14s %14s %14s %16s\n", "b", "...1", "...3", "...5", "...7",
        "bias D(b) = avg(3,5,7) - c1");
    for (int b = 2; b <= K; ++b) {
        u64 c1 = S.get(b, IDX_LAST + 1), c3 = S.get(b, IDX_LAST + 3);
        u64 c5 = S.get(b, IDX_LAST + 5), c7 = S.get(b, IDX_LAST + 7);
        t.p("%4d %14llu %14llu %14llu %14llu %+16.1f\n", b, c1, c3, c5, c7,
            (c3 + c5 + c7) / 3.0 - (double)c1);
    }
    {
        u64 c1 = 0, c3 = 0, c5 = 0, c7 = 0;
        for (int b = 1; b <= K; ++b) {
            c1 += S.get(b, IDX_LAST + 1);
            c3 += S.get(b, IDX_LAST + 3);
            c5 += S.get(b, IDX_LAST + 5);
            c7 += S.get(b, IDX_LAST + 7);
        }
        t.p("%4s %14llu %14llu %14llu %14llu %+16.1f\n", "all", c1, c3, c5, c7,
            (c3 + c5 + c7) / 3.0 - (double)c1);
    }
    t.p("\n");

    // -- 4. leading octal digit ---------------------------------------------------
    t.p("4. LEADING OCTAL DIGIT OF PRIMES vs Benford(base 8)\n");
    t.p("%6s", "b");
    for (int d = 1; d <= 7; ++d) t.p(" %9s%d", "lead ", d);
    t.p("\n");
    for (int b = 2; b <= K; ++b) {
        u64 primes = S.get(b, IDX_PRIMES);
        t.p("%6d", b);
        for (int d = 1; d <= 7; ++d)
            t.p(" %9.3f%%", primes ? 100.0 * S.get(b, IDX_LEAD + d) / primes : 0.0);
        t.p("\n");
    }
    t.p("%6s", "B8");
    for (int d = 1; d <= 7; ++d)
        t.p(" %9.3f%%", 100.0 * std::log1p(1.0 / d) / std::log(8.0));
    t.p("   <- Benford base-8 reference\n\n");

    // -- 5. octal digit frequency ------------------------------------------------
    t.p("5. OCTAL DIGIT FREQUENCY ACROSS ALL DIGITS OF ALL PRIMES\n");
    {
        u64 df[8] = {}, tot = 0;
        for (int b = 1; b <= K; ++b)
            for (int d = 0; d < 8; ++d) { df[d] += S.get(b, IDX_DF + d); tot += S.get(b, IDX_DF + d); }
        t.p("%8s", "digit");
        for (int d = 0; d < 8; ++d) t.p(" %9d", d);
        t.p("\n%8s", "share");
        for (int d = 0; d < 8; ++d) t.p(" %8.3f%%", tot ? 100.0 * df[d] / tot : 0.0);
        t.p("\n\n");
    }

    // -- 6. last two octal digits -------------------------------------------------
    t.p("6. LAST TWO OCTAL DIGITS (p mod 64) — extremes over all primes\n");
    {
        u64 m64[64] = {};
        for (int b = 1; b <= K; ++b)
            for (int r = 0; r < 64; ++r) m64[r] += S.get(b, IDX_M64 + r);
        std::vector<int> idx;
        for (int r = 1; r < 64; r += 2) idx.push_back(r);  // odd residues only
        std::sort(idx.begin(), idx.end(), [&](int a, int b2) { return m64[a] > m64[b2]; });
        t.p("   most common: ");
        for (int i = 0; i < 5; ++i)
            t.p("..%s (%llu)  ", oct_s(idx[i]).c_str(), m64[idx[i]]);
        t.p("\n   least common:");
        for (int i = (int)idx.size() - 5; i < (int)idx.size(); ++i)
            t.p(" ..%s (%llu) ", oct_s(idx[i]).c_str(), m64[idx[i]]);
        t.p("\n\n");
    }

    // -- 7. gaps -------------------------------------------------------------------
    if (postpass) {
        t.p("7. PRIME GAPS THROUGH THE OCTAL LENS\n");
        u64 evengaps = 0;
        for (int r = 0; r < 8; ++r) evengaps += post.gap_mod8[r];
        evengaps -= post.gap_mod8[1];  // the single 2->3 gap of 1
        t.p("   gap mod 8 (even gaps):");
        for (int r = 0; r < 8; r += 2)
            t.p("  %d: %.3f%%", r, evengaps ? 100.0 * post.gap_mod8[r] / evengaps : 0.0);
        t.p("\n   most common gaps: ");
        {
            std::vector<int> gi;
            for (int g = 1; g < GAP_HIST_MAX; ++g)
                if (post.gap_hist[g]) gi.push_back(g);
            std::sort(gi.begin(), gi.end(),
                      [&](int a, int b2) { return post.gap_hist[a] > post.gap_hist[b2]; });
            for (int i = 0; i < (int)gi.size() && i < 8; ++i)
                t.p("%d (%llu)  ", gi[i], post.gap_hist[gi[i]]);
        }
        t.p("\n   twin primes: %llu | max gap: %llu after %llu (octal %s)\n", post.twins_total,
            post.max_gap, post.max_gap_after, oct_s(post.max_gap_after).c_str());
        t.p("%6s %14s %18s %26s\n", "b", "twins", "max gap", "after prime (octal)");
        for (int b = 1; b <= K; ++b)
            t.p("%6d %14llu %18llu %26s\n", b, post.band_twins[b - 1],
                post.band_max_gap[b - 1], oct_s(post.band_max_gap_after[b - 1]).c_str());
        t.p("\n");

        // -- 8. palindromes ---------------------------------------------------------
        t.p("8. OCTAL PALINDROMIC PRIMES (digits read the same both ways in base 8)\n");
        t.p("   total: %llu  (no 2-digit octal palindromic prime exists: 0oDD = 9*D)\n",
            post.palin_total);
        t.p("%6s %12s   %s\n", "b", "count", "first examples (octal)");
        for (int b = 1; b <= K; ++b) {
            t.p("%6d %12llu   ", b, post.band_palin[b - 1]);
            for (u64 v : post.palin_examples[b - 1]) t.p("%s ", oct_s(v).c_str());
            t.p("\n");
        }
        t.p("\n");
    }

    t.p("9. THE BASE-8 IDENTITIES BEHIND THE PREDICTOR\n");
    t.p("   8 ≡ 1 (mod 7): n ≡ its octal digit sum (mod 7)\n");
    t.p("   8 ≡ -1 (mod 3): n ≡ its alternating octal digit sum (mod 3)\n");
    t.p("   8^i mod 5 cycles 1,3,4,2: weighted octal digit sum gives n mod 5\n");
    t.p("   last octal digit carries n mod 8; odd squares all end in octal 1\n");

    // -- CSVs ----------------------------------------------------------------------
    if (FILE *f = open_csv(out_dir, "pi_emergence.csv",
                           "band,range_lo,range_hi,primes,li_estimate,pi_cumulative")) {
        for (int b = 1; b <= K; ++b) {
            u64 lo = (b == 1) ? 1 : pow8(b - 1), hi = pow8(b);
            fprintf(f, "%d,%llu,%llu,%llu,%.3f,%llu\n", b, lo, hi, S.get(b, IDX_PRIMES),
                    Li((double)hi) - Li((double)lo), S.pi_upto_band(b));
        }
        fclose(f);
    }
    if (FILE *f = open_csv(out_dir, "predictor.csv",
                           "band,odd_total,candidates,primes,cand_primes")) {
        for (int b = 1; b <= K; ++b)
            fprintf(f, "%d,%llu,%llu,%llu,%llu\n", b, S.get(b, IDX_ODD_TOTAL),
                    S.get(b, IDX_CAND), S.get(b, IDX_PRIMES), S.get(b, IDX_CANDP));
        fclose(f);
    }
    if (FILE *f = open_csv(out_dir, "last_digit.csv", "band,d1,d3,d5,d7")) {
        for (int b = 1; b <= K; ++b)
            fprintf(f, "%d,%llu,%llu,%llu,%llu\n", b, S.get(b, IDX_LAST + 1),
                    S.get(b, IDX_LAST + 3), S.get(b, IDX_LAST + 5), S.get(b, IDX_LAST + 7));
        fclose(f);
    }
    if (FILE *f = open_csv(out_dir, "leading_digit.csv", "band,d1,d2,d3,d4,d5,d6,d7")) {
        for (int b = 1; b <= K; ++b) {
            fprintf(f, "%d", b);
            for (int d = 1; d <= 7; ++d) fprintf(f, ",%llu", S.get(b, IDX_LEAD + d));
            fprintf(f, "\n");
        }
        fclose(f);
    }
    if (FILE *f = open_csv(out_dir, "digit_freq.csv", "band,d0,d1,d2,d3,d4,d5,d6,d7")) {
        for (int b = 1; b <= K; ++b) {
            fprintf(f, "%d", b);
            for (int d = 0; d < 8; ++d) fprintf(f, ",%llu", S.get(b, IDX_DF + d));
            fprintf(f, "\n");
        }
        fclose(f);
    }
    if (FILE *f = open_csv(out_dir, "mod64.csv", "residue,octal,count")) {
        for (int r = 1; r < 64; r += 2) {
            u64 c = 0;
            for (int b = 1; b <= K; ++b) c += S.get(b, IDX_M64 + r);
            fprintf(f, "%d,%s,%llu\n", r, oct_s(r).c_str(), c);
        }
        fclose(f);
    }
    if (postpass) {
        if (FILE *f = open_csv(out_dir, "gap_hist.csv", "gap,count")) {
            for (int g = 1; g <= GAP_HIST_MAX; ++g)
                if (post.gap_hist[g]) fprintf(f, "%d,%llu\n", g, post.gap_hist[g]);
            fclose(f);
        }
        if (FILE *f = open_csv(out_dir, "band_gaps_twins.csv",
                               "band,twins,max_gap,max_gap_after")) {
            for (int b = 1; b <= K; ++b)
                fprintf(f, "%d,%llu,%llu,%llu\n", b, post.band_twins[b - 1],
                        post.band_max_gap[b - 1], post.band_max_gap_after[b - 1]);
            fclose(f);
        }
        if (FILE *f = open_csv(out_dir, "palindromes.csv", "band,count,examples_octal")) {
            for (int b = 1; b <= K; ++b) {
                fprintf(f, "%d,%llu,", b, post.band_palin[b - 1]);
                for (size_t i = 0; i < post.palin_examples[b - 1].size(); ++i)
                    fprintf(f, "%s%s", i ? " " : "", oct_s(post.palin_examples[b - 1][i]).c_str());
                fprintf(f, "\n");
            }
            fclose(f);
        }
    }

    if (t.f) {
        fclose(t.f);
        printf("\nreport: %s | CSVs in %s/\n", report_path.c_str(), out_dir.c_str());
    }

    // Hard self-verification against the known pi(8^k) table.
    for (int k = 1; k <= K && k <= 12; ++k)
        if (S.pi_upto_band(k) != PI_8_POW[k]) {
            fprintf(stderr, "WARNING: pi(8^%d) = %llu disagrees with known %llu\n", k,
                    S.pi_upto_band(k), PI_8_POW[k]);
            return 3;
        }
    return 0;
}
