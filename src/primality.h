// primality.h — host-side deterministic Miller-Rabin for u64 (ground truth for --judge).
#pragma once
#include "octal_core.h"

inline u64 mulmod_u64(u64 a, u64 b, u64 m) {
    return (u64)((unsigned __int128)a * b % m);
}

inline u64 powmod_u64(u64 a, u64 e, u64 m) {
    u64 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod_u64(r, a, m);
        a = mulmod_u64(a, a, m);
        e >>= 1;
    }
    return r;
}

// Deterministic for all n < 2^64 with bases {2,3,5,7,11,13,17,19,23,29,31,37}.
inline bool is_prime_u64(u64 n) {
    static const u64 bases[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return false;
    for (u64 p : bases)
        if (n % p == 0) return n == p;
    u64 d = n - 1;
    int r = 0;
    while (!(d & 1)) { d >>= 1; ++r; }
    for (u64 a : bases) {
        u64 x = powmod_u64(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int i = 0; i < r - 1; ++i) {
            x = mulmod_u64(x, x, n);
            if (x == n - 1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}
