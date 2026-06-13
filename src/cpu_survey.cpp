// cpu_survey.cpp — where do primes land in the octal wheel vs the hex wheel?
//
// Goal (project brief): find the geometric pattern of prime locations in base 8,
// compare it to base 16, and take the delta. This is a CPU companion to the CUDA
// survey so the experiment runs anywhere (no GPU required).
//
// The "geometry" is the modular wheel: lay the integers around a circle with `b`
// spokes (a clock with b hours). A number n sits on spoke (n mod b). A prime p>b
// can only sit on a spoke s with gcd(s,b)=1 — otherwise p would be divisible by a
// factor it shares with the base. Those coprime spokes are the prime corridors.
//
// We measure the corridors for base 8 and base 16 (and the 2^k family, plus the
// primorial / decimal bases for contrast) directly from a sieve.
//
// Build:  c++ -O3 -std=c++17 -o bin/cpu_survey src/cpu_survey.cpp
// Run:    ./bin/cpu_survey [N]        (default N = 1e9)

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

using std::uint64_t;
using std::uint32_t;

// ---------------------------------------------------------------------------
// odds-only sieve of Eratosthenes: bit k represents the odd number 2k+1.
// ---------------------------------------------------------------------------
struct Sieve {
    uint64_t N;
    std::vector<uint64_t> comp;  // composite flag for odd numbers; bit k -> 2k+1
    explicit Sieve(uint64_t n) : N(n) {
        uint64_t halfN = (N >> 1) + 1;             // indices 0..(N-1)/2
        comp.assign((halfN + 63) / 64, 0);
        set(0);                                     // 1 is not prime
        for (uint64_t p = 3; p * p <= N; p += 2) {
            if (test((p - 1) / 2)) continue;        // p already composite
            for (uint64_t m = p * p; m <= N; m += 2 * p) set((m - 1) / 2);
        }
    }
    inline void set(uint64_t k) { comp[k >> 6] |= (1ULL << (k & 63)); }
    inline bool test(uint64_t k) const { return comp[k >> 6] >> (k & 63) & 1; }
    // visit every prime p in [2, N]
    template <class F> void for_each_prime(F f) const {
        f(2);
        for (uint64_t k = 1; 2 * k + 1 <= N; ++k)
            if (!test(k)) f(2 * k + 1);
    }
};

static int octal_digits(uint64_t n) {            // number of base-8 digits
    int d = 0; do { ++d; n >>= 3; } while (n); return d;
}
static int gcd_i(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static int phi(int m) { int c = 0; for (int i = 1; i <= m; ++i) c += (gcd_i(i, m) == 1); return c; }

static const int MAXBAND = 24;

int main(int argc, char** argv) {
    uint64_t N = (argc > 1) ? strtoull(argv[1], nullptr, 10) : 1000000000ULL;

    fprintf(stderr, "sieving to N = %llu ...\n", (unsigned long long)N);
    Sieve s(N);

    // last-digit histograms for the bases we care about
    std::vector<uint64_t> m8(8, 0), m16(16, 0), m32(32, 0);
    std::vector<uint64_t> m6(6, 0), m10(10, 0), m30(30, 0), m210(210, 0);
    // band-resolved (octal band) octal & hex last digit, to watch the race vs scale
    static uint64_t band8[MAXBAND][8] = {};
    static uint64_t band16[MAXBAND][16] = {};
    uint64_t pi = 0;

    s.for_each_prime([&](uint64_t p) {
        ++pi;
        m8[p & 7]++; m16[p & 15]++; m32[p & 31]++;
        m6[p % 6]++; m10[p % 10]++; m30[p % 30]++; m210[p % 210]++;
        int b = octal_digits(p);
        if (b < MAXBAND) { band8[b][p & 7]++; band16[b][p & 15]++; }
    });

    // ---------------------------------------------------------------- report
    printf("prime-octal — octal wheel vs hex wheel (CPU survey)\n");
    printf("====================================================\n");
    printf("N = %llu   pi(N) = %llu\n\n", (unsigned long long)N, (unsigned long long)pi);

    auto pct = [&](uint64_t c) { return 100.0 * (double)c / (double)pi; };

    // -- octal wheel --
    printf("1. OCTAL WHEEL  (last octal digit = p mod 8)\n");
    printf("   spoke  coprime?   count        share     note\n");
    for (int d = 0; d < 8; ++d) {
        const char* note = (gcd_i(d, 8) != 1) ? "forbidden (shares 2 with base)"
                          : (d == 1)          ? "quadratic-residue class (odd^2 = 1 mod 8)"
                                              : "prime corridor";
        printf("   %d      %-8s  %12llu  %7.4f%%   %s\n", d,
               gcd_i(d, 8) == 1 ? "yes" : "no", (unsigned long long)m8[d], pct(m8[d]), note);
    }
    printf("   -> %d of 8 spokes carry primes (phi(8)=%d) = %.1f%% admissible density\n\n",
           phi(8), phi(8), 100.0 * phi(8) / 8);

    // -- hex wheel --
    printf("2. HEX WHEEL  (last hex digit = p mod 16)\n");
    printf("   spoke  coprime?   count        share     note\n");
    for (int d = 0; d < 16; ++d) {
        bool cop = gcd_i(d, 16) == 1;
        bool qr  = (d == 1 || d == 9);  // odd squares are 1 or 9 (mod 16)
        const char* note = !cop ? "forbidden (even)"
                         : qr   ? "quadratic-residue class (odd^2 in {1,9} mod 16)"
                                : "prime corridor";
        printf("   %2d     %-8s  %12llu  %7.4f%%   %s\n", d,
               cop ? "yes" : "no", (unsigned long long)m16[d], pct(m16[d]), note);
    }
    printf("   -> %d of 16 spokes carry primes (phi(16)=%d) = %.1f%% admissible density\n\n",
           phi(16), phi(16), 100.0 * phi(16) / 16);

    // -- the delta: octal -> hex refinement --
    printf("3. THE DELTA  (each octal spoke splits into two hex spokes: d and d+8)\n");
    printf("   the only new information is bit 3 (the 8's bit) of p.\n");
    printf("   octal d   count       ->  hex d      hex d+8     low/high split\n");
    for (int d : {1, 3, 5, 7}) {
        uint64_t lo = m16[d], hi = m16[d + 8];
        printf("   %d        %12llu  ->  %11llu %11llu   %6.3f%% / %6.3f%%\n",
               d, (unsigned long long)m8[d], (unsigned long long)lo, (unsigned long long)hi,
               100.0 * lo / m8[d], 100.0 * hi / m8[d]);
    }
    printf("   reading: the extra hex bit splits every octal corridor ~50/50.\n");
    printf("   the hex wheel is the octal wheel at 2x angular resolution -- no new\n");
    printf("   prime structure, just one more bit of p revealed.\n\n");

    // -- Chebyshev / quadratic-residue bias --
    printf("4. CHEBYSHEV BIAS  (which corridors run light?)\n");
    {
        // octal corridors ranked
        std::vector<int> oc = {1, 3, 5, 7};
        std::sort(oc.begin(), oc.end(), [&](int a, int b){ return m8[a] < m8[b]; });
        printf("   octal corridors fewest->most: ");
        for (int d : oc) printf("%d(%.4f%%) ", d, pct(m8[d]));
        printf("\n   -> spoke 1 is the quadratic-residue class and trails.\n");
        std::vector<int> hc = {1, 3, 5, 7, 9, 11, 13, 15};
        std::sort(hc.begin(), hc.end(), [&](int a, int b){ return m16[a] < m16[b]; });
        printf("   hex corridors  fewest->most: ");
        for (int d : hc) printf("%d(%.4f%%) ", d, pct(m16[d]));
        printf("\n   -> spokes 1 and 9 (the QR classes) trail; and octal-1 = hex-1 + hex-9,\n");
        printf("      so the octal deficit splits exactly into the two deficient hex spokes.\n\n");
    }

    // -- why powers of two plateau: admissible density across bases --
    printf("5. WHY OCTAL=HEX AT THE LAST DIGIT: admissible density phi(b)/b\n");
    printf("   a base's last digit only screens the primes that divide the base.\n");
    printf("   base b   factorization      phi(b)/b     primes screened by last digit\n");
    struct Row { int b; const char* f; const char* scr; };
    Row rows[] = {
        {2,  "2",        "2"},
        {8,  "2^3",      "2  (only!)"},
        {16, "2^4",      "2  (only!)"},
        {32, "2^5",      "2  (only!)"},
        {10, "2*5",      "2, 5"},
        {6,  "2*3",      "2, 3"},
        {30, "2*3*5",    "2, 3, 5"},
        {210,"2*3*5*7",  "2, 3, 5, 7"},
    };
    for (auto& r : rows)
        printf("   %-7d  %-16s   %6.2f%%      %s\n", r.b, r.f, 100.0 * phi(r.b) / r.b, r.scr);
    printf("   -> every power of two sticks at 50%%: it only ever screens the prime 2.\n");
    printf("      octal and hex are informationally identical at the last digit.\n");
    printf("      the real corridors tighten with PRIMORIAL bases (6,30,210,...),\n");
    printf("      which is exactly what the base-8 *multi-digit* predictor reconstructs.\n\n");

    // -- equidistribution sanity (primorial base 30) --
    printf("6. EQUIDISTRIBUTION in a primorial wheel (base 30, the 8 corridors)\n");
    printf("   ");
    for (int d = 0; d < 30; ++d) if (gcd_i(d, 30) == 1) printf("%d:%.3f%% ", d, pct(m30[d]));
    printf("\n   -> Dirichlet: primes spread evenly across the coprime spokes.\n\n");

    // ---------------------------------------------------------------- CSV out
    // long format the python viz reads: base,residue,coprime,count
    FILE* f = fopen("results/digits.csv", "w");
    if (f) {
        fprintf(f, "base,residue,coprime,count\n");
        auto dump = [&](int base, const std::vector<uint64_t>& h) {
            for (int d = 0; d < base; ++d)
                fprintf(f, "%d,%d,%d,%llu\n", base, d, gcd_i(d, base) == 1,
                        (unsigned long long)h[d]);
        };
        dump(8, m8); dump(16, m16); dump(32, m32);
        dump(6, m6); dump(10, m10); dump(30, m30);
        fclose(f);
        fprintf(stderr, "wrote results/digits.csv\n");
    }
    FILE* fb = fopen("results/bands.csv", "w");
    if (fb) {
        fprintf(fb, "band,base,residue,count\n");
        for (int b = 1; b < MAXBAND; ++b) {
            uint64_t tot = 0; for (int d = 0; d < 8; ++d) tot += band8[b][d];
            if (!tot) continue;
            for (int d = 0; d < 8; ++d)  fprintf(fb, "%d,8,%d,%llu\n",  b, d, (unsigned long long)band8[b][d]);
            for (int d = 0; d < 16; ++d) fprintf(fb, "%d,16,%d,%llu\n", b, d, (unsigned long long)band16[b][d]);
        }
        fclose(fb);
        fprintf(stderr, "wrote results/bands.csv\n");
    }
    // headline numbers for the writeup
    FILE* fs = fopen("results/summary.csv", "w");
    if (fs) {
        fprintf(fs, "key,value\n");
        fprintf(fs, "N,%llu\n", (unsigned long long)N);
        fprintf(fs, "pi,%llu\n", (unsigned long long)pi);
        fprintf(fs, "oct_admissible_pct,%.4f\n", 100.0 * phi(8) / 8);
        fprintf(fs, "hex_admissible_pct,%.4f\n", 100.0 * phi(16) / 16);
        fclose(fs);
    }
    return 0;
}
