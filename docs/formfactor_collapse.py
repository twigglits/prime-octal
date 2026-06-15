import re, math, numpy as np
txt=open('/tmp/claude-1000/-home-jeannaude-Documents-prime-octal/63550034-ca76-4301-9489-80fa38ef6280/tasks/bhdfuq5pk.output').read()
rows=[]; R=None; L=None
for line in txt.splitlines():
    m=re.match(r'=== R=(\d+) L=(\w+)',line)
    if m: R=int(m.group(1)); L=m.group(2); continue
    p=line.split()
    if len(p)==6 and p[0].isdigit():
        K,lam,psi_p,psi_c,delta,D=int(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4]),float(p[5])
        rows.append((R,L,K,lam,psi_p,psi_c,delta))

print("=== 1. UNIVERSALITY COLLAPSE: Psi_prime vs delta, across R and lattice ===")
# bin by delta, report spread (std) within bins -> small spread = collapse
import collections
byd=collections.defaultdict(list)
for R,L,K,lam,pp,pc,d in rows: byd[round(d,3)].append(pp)
print(f"{'delta':>7} {'n':>3} {'mean Psi':>9} {'std':>7}  (small std across R&lattice = collapse)")
for d in sorted(byd):
    v=byd[d]; print(f"{d:>7.3f} {len(v):>3} {np.mean(v):>9.3f} {np.std(v):>7.3f}")

print("\n=== 2. Is it CUE? Compare measured number variance V=Psi*lam to Poisson and CUE ===")
print(f"{'lam':>9} {'V_data':>10} {'V_Poisson':>10} {'V_CUE':>9} {'data/Pois':>9} {'data/CUE':>10}")
for R,L,K,lam,pp,pc,d in sorted(rows,key=lambda r:-r[3])[:6]:
    V=pp*lam; Vp=lam
    Vcue=(1/math.pi**2)*(math.log(2*math.pi*lam)+0.5772+1)  # CUE number variance (large-lam)
    print(f"{lam:>9.1f} {V:>10.1f} {Vp:>10.1f} {Vcue:>9.2f} {V/Vp:>9.3f} {V/Vcue:>10.1f}")

print("\n=== 3. Ramp fit Psi=min(s*delta+b,1): collapse master curve ===")
d=np.array([r[6] for r in rows]); y=np.array([r[4] for r in rows])
mask=y<0.96  # unsaturated part
A=np.vstack([d[mask],np.ones(mask.sum())]).T
s,b=np.linalg.lstsq(A,y[mask],rcond=None)[0]
print(f"  unsaturated slope s={s:.2f}, intercept b={b:.2f}; saturates to ~1 at delta~{(1-b)/s:.2f}")
resid=y[mask]-(s*d[mask]+b)
print(f"  residual RMS about the master ramp = {np.sqrt(np.mean(resid**2)):.3f}")
