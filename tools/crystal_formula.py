#!/usr/bin/env python3
"""From 'the composites are a predictable crystal' to an exact prime formula --
and the honest accounting of what that formula costs.

Produces:
  - a verified printout of every claim in docs/the-crystal-formula.md
  - docs/fig-crystal-sieve.svg : the holes of a FIXED wheel are primes + impostors;
    only when the wheel grows to all primes <= sqrt(n) do the holes become exactly
    the primes.
Pure stdlib.
"""
from math import gcd, isqrt, log, factorial
import bisect, os

OUT = os.path.join(os.path.dirname(__file__), "..", "docs")

def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = b"\x00" * len(s[i*i::i])
    return s

N = 2_000_000
sv = sieve(N)
primes = [i for i in range(2, N + 1) if sv[i]]

# ---- exact prime predicate from the GROWING crystal: gcd(n, prod_{p<=sqrt n} p) ----
pp, pv, prod = [], [], 1
for p in primes:
    if p > isqrt(N): break
    prod *= p; pp.append(p); pv.append(prod)
def M_of(n):
    k = bisect.bisect_right(pp, isqrt(n))
    return pv[k-1] if k else 1
def is_prime_crystal(n):
    return n >= 2 and gcd(n, M_of(n)) == 1

# ---- Wilson indicator: exact, self-contained (no prior primes), but O(n) ----
def wilson_isprime(j):
    return 1 if j >= 2 and (factorial(j-1) % j) == j-1 else 0
def pi_wilson(n):
    return sum(wilson_isprime(j) for j in range(2, n+1))
def nth_prime_formula(n):
    # Willans (1964), exact closed form:  p_n = 1 + sum_{m=1}^{2^n} [ pi(m) < n ].
    # The summand uses ONLY the Wilson indicator, so no prior prime list is needed.
    # Equivalent to: smallest m with pi_wilson(m) = n. Computed that way here.
    m, hits = 1, 0
    while hits < n:
        m += 1
        hits += wilson_isprime(m)
    return m

def report():
    print("CLAIM 1  finite crystal mod 210 is exactly periodic:",
          all((gcd(n,210)>1)==(gcd(n+210,210)>1) for n in range(2,5000)),
          " | holes/period = phi(210) =", sum(gcd(r,210)==1 for r in range(210)))

    print("\nCLAIM 2  smallest impostor (composite hole) of each fixed wheel = (next prime)^2:")
    for wp in ([2],[2,3],[2,3,5],[2,3,5,7],[2,3,5,7,11]):
        P=1
        for p in wp: P*=p
        n=2
        while not (gcd(n,P)==1 and not sv[n]): n+=1
        print(f"   wheel {'*'.join(map(str,wp)):12s} (P={P:5d}):  smallest impostor {n:4d} = {primes[len(wp)]}^2")
    print("   210-wheel precision (holes truly prime) vs asymptotic 4.375/ln N:")
    for M in (10**3,10**4,10**5,10**6):
        holes=[n for n in range(2,M+1) if gcd(n,210)==1]
        pr=sum(sv[n] for n in holes)/len(holes)
        print(f"     N={M:>8d}: measured {pr:.4f}   asymptotic {4.375/log(M):.4f}")

    print("\nCLAIM 3  EXACT:  n prime  <=>  gcd(n, prod_{p<=sqrt n} p) = 1")
    print("         verified against the sieve for all n < 200000:",
          all(is_prime_crystal(n)==bool(sv[n]) for n in range(2,200000)))
    k=bisect.bisect_right(pp,isqrt(10**6))
    print(f"         but the modulus grows: M(10^6) has {len(str(pv[k-1]))} digits "
          f"(log M ~ sqrt n); building it already needs the primes <= sqrt n.")

    print("\nCLAIM 4  no FIXED period can isolate primes (primality is not periodic):")
    for T in (210,2310,30030):
        q=next(q for q in primes if q>T and not sv[q+T])
        print(f"   period {T:6d} broken: {q} prime, {q}+{T}={q+T} composite")

    print("\nCLAIM 5  exact self-contained formulas (no prior primes), but costly:")
    print("   Wilson (n-1)!=-1 mod n correct for n<300:",
          all(bool(wilson_isprime(n))==bool(sv[n]) for n in range(2,300)))
    got=[nth_prime_formula(n) for n in range(1,16)]
    print("   n-th prime via the Wilson/Willans counting formula, p_1..p_15:")
    print("     ", got, "==", primes[:15], "->", got==primes[:15])

# ----------------------------------------------------------------- figure ---- #
def fig():
    K = 330                                   # show n = 2..K
    wheels = [("wheel 6  (2·3)",[2,3]),
              ("wheel 30  (2·3·5)",[2,3,5]),
              ("wheel 210  (2·3·5·7)",[2,3,5,7]),
              ("wheel 2310  (·11)",[2,3,5,7,11]),
              ("EXACT: all primes ≤ √n",None)]
    cw, ch, lab = 4, 30, 150                   # cell w/h, label gutter
    W = lab + (K-1)*cw + 20
    H = 50 + len(wheels)*(ch+8) + 40
    GREEN="#2a9d4a"; RED="#d6283c"; GREY="#e0e3e7"; WHITE="#ffffff"
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'viewBox="0 0 {W} {H}" font-family="system-ui,Arial"><rect width="{W}" height="{H}" fill="white"/>']
    s.append(f'<text x="{W/2}" y="26" text-anchor="middle" font-size="16" font-weight="700" '
             f'fill="#222">Holes of a fixed crystal = primes + impostors; only the growing crystal gives primes exactly</text>')
    for bi,(name,wp) in enumerate(wheels):
        y=50+bi*(ch+8)
        if wp:
            P=1
            for p in wp: P*=p
        s.append(f'<text x="{lab-8}" y="{y+ch*0.68:.0f}" text-anchor="end" font-size="11.5" fill="#333">{name}</text>')
        for n in range(2,K+1):
            x=lab+(n-2)*cw
            if wp is None:                      # exact band: hole <=> prime
                col = GREEN if sv[n] else GREY
            else:
                hole = gcd(n,P)==1
                if not hole: col=GREY           # screened by the wheel
                elif sv[n]:   col=GREEN         # a real prime
                else:         col=RED           # IMPOSTOR (rough composite)
            s.append(f'<rect x="{x}" y="{y}" width="{cw-0.4:.1f}" height="{ch}" fill="{col}"/>')
    # legend
    ly=H-22
    for i,(c,t) in enumerate([(GREEN,"prime"),(RED,"impostor (composite in a hole)"),(GREY,"screened / not a hole")]):
        lx=lab+i*230
        s.append(f'<rect x="{lx}" y="{ly-11}" width="14" height="14" fill="{c}" stroke="#bbb"/>')
        s.append(f'<text x="{lx+20}" y="{ly}" font-size="12" fill="#333">{t}</text>')
    s.append('</svg>')
    p=os.path.join(OUT,"fig-crystal-sieve.svg"); open(p,"w").write("".join(s))
    print("\nwrote",os.path.normpath(p))

if __name__=="__main__":
    os.makedirs(OUT,exist_ok=True)
    report()
    fig()
