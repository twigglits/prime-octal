import numpy as np

def sieve(n):
    b=np.ones(n+1,bool); b[:2]=False
    for i in range(2,int(n**0.5)+1):
        if b[i]: b[i*i::i]=False
    return b

def prime_factors(m):
    f=set(); d=2
    while d*d<=m:
        while m%d==0: f.add(d); m//=d
        d+=1
    if m>1: f.add(m)
    return f

# Free digit-rule primes for base b:
#   last digit -> primes dividing b
#   digit sum  -> primes dividing b-1
#   alt sum    -> primes dividing b+1
#   ALSO 2 is always catchable if b even (last digit parity)
def free_primes(b):
    s=set()
    s|=prime_factors(b)        # last-digit residue mod b catches these
    s|=prime_factors(b-1)
    s|=prime_factors(b+1)
    return s

N=20_000_000
isp=sieve(N)
nums=np.arange(N+1)
primes_mask=isp.copy()
n_primes=int(primes_mask.sum())

for b in (8,16):
    P=sorted(free_primes(b))
    wheel=1
    for p in P: wheel*=p
    # candidate = coprime to all detected primes
    cand=np.ones(N+1,bool); cand[:2]=False
    for p in P:
        cand[::p]=False
        cand[p]=True  # keep the prime p itself as "not a candidate" -> actually rule rejects multiples incl p
    # rebuild cleanly: candidate iff gcd(n,wheel)==1
    cand=np.ones(N+1,bool); cand[0]=False
    for p in P: cand[::p]=False
    n_cand=int(cand.sum())
    # among candidates, how many prime; recall of primes (primes>maxP kept)
    cand_primes=int((cand & primes_mask).sum())
    keep_frac=n_cand/N
    concentration=(cand_primes/n_cand)/(n_primes/N)
    # primes missed = primes that are NOT candidates and not in P themselves
    missed=int((primes_mask & ~cand).sum())
    missed_small=sum(1 for p in P if p<=N)
    print(f"base {b:2d}: free primes {P}  wheel={wheel}")
    print(f"          keep {keep_frac*100:6.3f}% of integers | concentration {concentration:.4f}x | candidate prime-density {cand_primes/n_cand*100:.3f}%")
    print(f"          primes missed by wheel = {missed} (all are the rule-primes {P}; recall above max={missed-missed_small})")
print(f"\n[N={N}, pi(N)={n_primes}]")
def order(b,p):
    if b%p==0: return None  # detected by trailing digits, "order 0"
    x=1%p
    for k in range(1,p):
        x=(x*b)%p
        if x==1: return k
    return p-1

primes=[2,3,5,7,11,13,17,19,23,29,31]
print(f"{'p':>3} | {'ord_8(p)':>9} {'ord_16(p)':>9}   (smaller order = more local / cheaper rule)")
print("-"*50)
w8=w16=0
for p in primes:
    o8=order(8,p); o16=order(16,p)
    s8 = "trailing" if o8 is None else str(o8)
    s16= "trailing" if o16 is None else str(o16)
    win = "8" if (o8 or 0)<(o16 or 0) else ("16" if (o16 or 0)<(o8 or 0) else "tie")
    print(f"{p:>3} | {s8:>9} {s16:>9}   octal {'<' if win=='8' else ('>' if win=='16' else '=')} hex")
# Push further: for each base b, which small primes get a LOCAL rule?
#   - "trailing" (p | b): order 0, cheapest
#   - digit-sum (p | b-1): order 1
#   - alt-sum  (p | b+1): order 2
# Define locality cost = ord_p(b). A base "cheaply captures" p if ord<=L.
# Question: which base captures the most of the first-k primes at low cost?

def order(b,p):
    if b%p==0: return 0
    x=1%p
    for k in range(1,p):
        x=(x*b)%p
        if x==1: return k
    return p-1

SMALL=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]

def capture(b,L):
    return [p for p in SMALL if order(b,p)<=L]

print("=== Which base best captures small primes with LOCAL rules (ord<=2: trailing/digit-sum/alt-sum)? ===")
best=[]
for b in range(3,65):
    caps=capture(b,2)
    # value = product of captured primes' "wheel removal" ~ sum of 1/p weighting big-prime cost,
    # but the meaningful metric: how many of the smallest consecutive primes captured (2,3,5,7,...)
    consec=0
    for p in SMALL:
        if p in caps: consec+=1
        else: break
    best.append((consec, sum(1/p for p in caps), b, caps))
best.sort(reverse=True)
print(f"{'base':>4} | {'consec small primes (2,3,5,7..)':>32} | local-captured set")
for consec,score,b,caps in best[:12]:
    pretty="2"+"".join(f"·{p}" if p in caps else "" for p in SMALL[1:consec])
    print(f"{b:>4} | first {consec:>2} primes {str(SMALL[:consec]):>22} | {caps}")

print("\n=== Locality cost of capturing {2,3,5,7} (the mod-210 wheel) per base ===")
for b in range(3,33):
    if all(p in [2,3,5,7] or True for p in [2,3,5,7]):
        cost=[order(b,p) for p in (2,3,5,7)]
        maxc=max(cost)
        tag=" <-- octal" if b==8 else (" <-- hex" if b==16 else "")
        if maxc<=4:
            print(f"base {b:>2}: ord(2,3,5,7)={cost}  max-span={maxc}{tag}")
