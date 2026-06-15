// moat_sweep.cu — column-sweep Gaussian moat with O(k·R) memory (no O(R²) bitmap).
//
// The full-disk method (src/lattice.cu --moat-gpu) needs an (R+1)² bitmap and
// walls out at R≈400k on 32 GB. The √26 component ends at |z|≈1,015,638
// (Tsuchimura), far past that. This engine removes the wall: with step length
// ≤ k (so |dx| ≤ Kc = ⌊√Ksq⌋), two adjacent primes differ in x by ≤ Kc, so a
// left-to-right column sweep only ever needs a window of Kc+1 columns resident.
//
// Connectivity is a sweep-line Hoshen–Kopelman union-find: each prime is a node,
// adjacent primes are unioned, and every CLUSTER (union-find root) carries its
// aggregates — max-norm, origin-flag, count-of-cells-still-in-window. A cell that
// leaves the window keeps contributing through its cluster's aggregate (merges
// take max/or), so departed cells never need relabeling. Records are
// reference-counted and recycled, so live memory is O(window cells) = O(Kc·R).
//
// Half-plane reduction: Gaussian primes are symmetric under (a,b)→(−a,b); any
// path's a<0 segments reflect into a≥0 with steps only getting shorter, so the
// origin component restricted to a≥0 is connected and reaches the same max |z|.
// We therefore sweep a = 0..R only and seed at (1,1).
//
// Primality is computed on the GPU one column-block at a time (Miller–Rabin on
// the norm a²+b²), so the GPU memory is just the block, not the disk.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <cuda_runtime.h>

typedef unsigned long long u64;

#define CUDA_OK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
    exit(2); } } while (0)

// --- device primality (deterministic MR for n < 3.3e18) ---------------------
__device__ __forceinline__ u64 d_mulmod(u64 a, u64 b, u64 m) { return (u64)((unsigned __int128)a * b % m); }
__device__ u64 d_powmod(u64 a, u64 e, u64 m) { u64 r = 1; a %= m; while (e) { if (e & 1) r = d_mulmod(r, a, m); a = d_mulmod(a, a, m); e >>= 1; } return r; }
__device__ bool d_isprime(u64 n) {
    const u64 B[12] = {2,3,5,7,11,13,17,19,23,29,31,37};
    if (n < 2) return false;
    for (int i = 0; i < 12; ++i) { u64 p = B[i]; if (n % p == 0) return n == p; }
    u64 d = n - 1; int r = 0; while (!(d & 1)) { d >>= 1; ++r; }
    for (int i = 0; i < 12; ++i) {
        u64 x = d_powmod(B[i], d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int j = 0; j < r - 1; ++j) { x = d_mulmod(x, x, n); if (x == n - 1) { comp = false; break; } }
        if (comp) return false;
    }
    return true;
}
__device__ bool d_gp(long long a, long long b) {
    long long A = a < 0 ? -a : a, B = b < 0 ? -b : b;
    if (A == 0 && B == 0) return false;
    if (A == 0 || B == 0) { u64 m = (u64)(A + B); return d_isprime(m) && (m % 4 == 3); }
    return d_isprime((u64)A * A + (u64)B * B);
}

// prime bits for W columns [a0, a0+W); b ∈ [−R,R]; only norm ≤ R². bit (w,bi).
__global__ void mark_cols(long long a0, long long W, long long R, long long Wstride, u64* out) {
    long long Bspan = 2 * R + 1;
    u64 tot = (u64)W * Bspan;
    for (u64 idx = (u64)blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += (u64)gridDim.x * blockDim.x) {
        long long w = idx / Bspan, bi = idx % Bspan, b = bi - R, a = a0 + w;
        if ((u64)a * a + (u64)b * b > (u64)R * R) continue;
        if (d_gp(a, b)) atomicOr((unsigned long long*)&out[(u64)w * Wstride + (bi >> 6)], 1ULL << (bi & 63));
    }
}

// --- reference-counted union-find with cluster aggregates -------------------
// ref[id] = (# window cells pointing to id) + (# records c≠id with par[c]==id).
struct DSU {
    std::vector<int> par, rnk, ref, fl;
    std::vector<char> org;
    std::vector<u64> mx;
    std::vector<long long> wc;
    int alloc(bool o, u64 n) {
        int id;
        if (!fl.empty()) { id = fl.back(); fl.pop_back(); }
        else { id = par.size(); par.push_back(0); rnk.push_back(0); ref.push_back(0); org.push_back(0); mx.push_back(0); wc.push_back(0); }
        par[id] = id; rnk[id] = 0; ref[id] = 1; org[id] = o; mx[id] = n; wc[id] = 1;
        return id;
    }
    void decref(int id) {                       // iterative cascade-free
        for (;;) {
            if (--ref[id] != 0) return;
            int p = par[id];
            fl.push_back(id);
            if (p == id) return;                // root: no parent reference to drop
            id = p;
        }
    }
    std::vector<int> path, dec;
    int find(int x) {
        int r = x; while (par[r] != r) r = par[r];
        path.clear();
        for (int n = x; par[n] != r; n = par[n]) path.push_back(n);
        for (int m : path) { dec.push_back(par[m]); par[m] = r; ref[r]++; }   // reparent all first
        for (int oldp : dec) decref(oldp);                                    // then drop old-parent refs
        dec.clear();
        return r;
    }
    // returns merged root
    int unite(int x, int y) {
        int rx = find(x), ry = find(y);
        if (rx == ry) return rx;
        if (rnk[rx] < rnk[ry]) { int t = rx; rx = ry; ry = t; }
        par[ry] = rx; ref[rx]++;
        if (rnk[rx] == rnk[ry]) rnk[rx]++;
        mx[rx] = mx[rx] > mx[ry] ? mx[rx] : mx[ry];
        org[rx] = org[rx] || org[ry];
        wc[rx] += wc[ry];
        return rx;
    }
};

struct SweepResult { double farthest; bool escaped; long long maxRecords; };

static SweepResult moat_sweep(long long Ksq, long long R, long long blockW, bool verbose) {
    const long long Kc = (long long)floor(sqrt((double)Ksq));
    const long long Bspan = 2 * R + 1, Wstride = (Bspan + 63) / 64, WC = Kc + 1;
    const u64 boundsq = (u64)(R - Kc) * (R - Kc);

    // step offsets (da≥1 backward cols; plus same-col db>0 backward)
    std::vector<std::pair<long long,long long>> back;     // (da≥1, db)
    for (long long da = 1; da <= Kc; ++da)
        for (long long db = -Kc; db <= Kc; ++db)
            if (da * da + db * db <= Ksq) back.push_back({da, db});
    std::vector<long long> same;                          // db≥1 within same column
    for (long long db = 1; db <= Kc; ++db) if (db * db <= Ksq) same.push_back(db);

    // GPU block buffer + host copy
    u64 *d_block; CUDA_OK(cudaMalloc(&d_block, (size_t)blockW * Wstride * 8));
    std::vector<u64> h_block((size_t)blockW * Wstride);
    long long blockStart = -1;
    auto column_words = [&](long long a) -> const u64* {
        if (blockStart < 0 || a >= blockStart + blockW) {
            blockStart = a;
            CUDA_OK(cudaMemset(d_block, 0, (size_t)blockW * Wstride * 8));
            u64 span = (u64)blockW * Bspan; int tpb = 256;
            u64 blocks = (span + tpb - 1) / tpb; if (blocks > 65535) blocks = 65535;
            mark_cols<<<(unsigned)blocks, tpb>>>(blockStart, blockW, R, Wstride, d_block);
            CUDA_OK(cudaGetLastError()); CUDA_OK(cudaDeviceSynchronize());
            CUDA_OK(cudaMemcpy(h_block.data(), d_block, (size_t)blockW * Wstride * 8, cudaMemcpyDeviceToHost));
        }
        return &h_block[(size_t)(a - blockStart) * Wstride];
    };

    std::vector<u64> winp((size_t)WC * Wstride);          // window prime bits
    std::vector<int> cell((size_t)WC * Bspan, -1);        // window cell -> record id
    auto isP = [&](long long col, long long b) -> bool {
        if (b < -R || b > R) return false;
        long long s = col % WC; u64 bi = b + R;
        return (winp[(size_t)s * Wstride + (bi >> 6)] >> (bi & 63)) & 1ULL;
    };
    auto cref = [&](long long col, long long b) -> int& {
        return cell[(size_t)(col % WC) * Bspan + (b + R)];
    };

    DSU d;
    bool escaped = false, bounded = false, born = false;
    int seedRid = -1; double farthest = 0; long long maxRecords = 0;

    for (long long a = 0; a <= R && !escaped && !bounded; ++a) {
        long long s = a % WC;
        // a column leaves the window: the slot we are about to reuse held column a-WC.
        if (a >= WC) {
            u64* w = &winp[(size_t)s * Wstride];
            for (u64 word = 0; word < (u64)Wstride; ++word) {
                u64 bits = w[word];
                while (bits) {
                    int bit = __builtin_ctzll(bits); bits &= bits - 1;
                    long long b = (long long)(word * 64 + bit) - R;
                    int rid = cref(a - WC, b);
                    if (rid < 0) continue;
                    int root = d.find(rid);
                    d.wc[root]--;
                    if (d.wc[root] == 0 && d.org[root]) { bounded = true; farthest = sqrt((double)d.mx[root]); }
                    cref(a - WC, b) = -1;
                    d.decref(rid);
                }
            }
        }
        if (bounded) break;
        // load column a prime bits into the window slot
        std::memcpy(&winp[(size_t)s * Wstride], column_words(a), (size_t)Wstride * 8);
        // add column a's primes: create records and union with window neighbours
        const u64* w = &winp[(size_t)s * Wstride];
        for (u64 word = 0; word < (u64)Wstride; ++word) {
            u64 bits = w[word];
            while (bits) {
                int bit = __builtin_ctzll(bits); bits &= bits - 1;
                long long b = (long long)(word * 64 + bit) - R;
                bool seed = (a == 1 && b == 1);
                int rid = d.alloc(seed, (u64)a * a + (u64)b * b);
                cref(a, b) = rid;
                if (seed) { born = true; seedRid = rid; d.ref[rid]++; }   // pin seed: stable origin handle
                for (auto& o : back) {                     // backward columns (a-da >= 0 guaranteed below)
                    long long ca = a - o.first; if (ca < 0) continue;
                    long long nb = b + o.second; int nr = isP(ca, nb) ? cref(ca, nb) : -1;
                    if (nr >= 0) d.unite(rid, nr);
                }
                for (long long db : same) {                // same column, already-added (b-db)
                    int nr = (isP(a, b - db)) ? cref(a, b - db) : -1;
                    if (nr >= 0) d.unite(rid, nr);
                }
            }
        }
        if (born) {
            int root = d.find(seedRid);                    // pinned -> always valid
            if (d.mx[root] >= boundsq) { escaped = true; farthest = sqrt((double)d.mx[root]); }
        }
        long long live = (long long)d.par.size() - (long long)d.fl.size();
        if (live > maxRecords) maxRecords = live;
    }
    if (!escaped && !bounded && born) {                    // swept to R without finalizing
        int root = d.find(seedRid); escaped = true; farthest = sqrt((double)d.mx[root]);
    }
    CUDA_OK(cudaFree(d_block));

    if (verbose) {
        printf("Column-sweep moat k=sqrt(%lld)=%.4f, R=%lld, window=%lld cols, peak records=%lld (%.1f MB)\n",
               Ksq, sqrt((double)Ksq), R, WC, maxRecords, maxRecords * 28.0 / 1e6 + cell.size() * 4.0 / 1e6);
        if (escaped)
            printf("  INCONCLUSIVE: origin component reached |z|=%.4f (>= R-Kc) — escapes this disk, rerun larger R\n", farthest);
        else
            printf("  BOUNDED: origin component sealed off, farthest |z|=%.4f -> moat confirmed\n", farthest);
    }
    return {farthest, escaped, maxRecords};
}

static int selftest() {
    // Validate against the full-disk BFS / Tsuchimura values, with O(k·R) memory.
    struct { long long Ksq, R; double far; } cases[] = {
        {2, 200, 11.7047}, {4, 200, 45.3100}, {8, 300, 93.4719},
    };
    for (auto& c : cases) {
        SweepResult r = moat_sweep(c.Ksq, c.R, 256, false);
        assert(!r.escaped && "small moat should be bounded");
        assert(fabs(r.farthest - c.far) < 1e-3 && "column-sweep farthest disagrees with known value");
    }
    // memory really is O(k·R): peak records must be tiny vs the component size.
    SweepResult big = moat_sweep(8, 2000, 256, false);
    assert(big.maxRecords < 200000 && "memory not bounded — recycling broken");
    printf("selftest OK  (k^2=2,4,8 farthest match known BFS values; peak records bounded)\n");
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "--selftest")) return selftest();
    if (argc >= 4 && !strcmp(argv[1], "--moat-sweep")) {
        long long Ksq = atoll(argv[2]), R = atoll(argv[3]);
        long long blockW = argc >= 5 ? atoll(argv[4]) : 256;
        moat_sweep(Ksq, R, blockW, true);
        return 0;
    }
    fprintf(stderr,
        "usage:\n"
        "  %s --selftest\n"
        "  %s --moat-sweep KSQ R [blockW]   column-sweep Gaussian moat, O(k*R) memory\n"
        "      KSQ = integer k^2 (e.g. 26 for the sqrt(26) moat); R = search radius\n",
        argv[0], argv[0]);
    return 1;
}
