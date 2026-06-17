#!/usr/bin/env python3
"""Render the prime wheels: octal vs hex vs primorial, plus the dyadic refinement.

Reads results/digits.csv (big-N corridor counts from cpu_survey) for the
quantitative figures, and sieves a small range itself for the prime-ray scatter.
Pure stdlib -- emits SVG (text) and an ASCII wheel. No numpy/matplotlib needed.
"""
import csv, math, os

OUT = os.path.join(os.path.dirname(__file__), "..", "docs")
RES = os.path.join(os.path.dirname(__file__), "..", "results")

# octal-corridor colour map: a prime is coloured by p mod 8 so the hex panel
# shows each octal corridor splitting into two adjacent hex rays (the "delta").
OCT_COLOR = {1: "#d1495b", 3: "#edae49", 5: "#00798c", 7: "#66a182"}
GREY = "#c9ccd1"


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = b"\x00" * len(sieve[i * i :: i])
    return [i for i in range(2, n + 1) if sieve[i]]


def coprime(d, b):
    return math.gcd(d, b) == 1


# --------------------------------------------------------------------------- #
#  Figure 1: prime-ray scatter -- THE geometric pattern of where primes land   #
# --------------------------------------------------------------------------- #
def ray_panel(cx, cy, R, base, primes, M, title, color_by_oct=True):
    """One wheel panel: spoke n mod base, radius grows with n. Only coprime
    spokes populate, so primes form rays with empty wedges between them."""
    s = []
    s.append(f'<text x="{cx}" y="{cy-R-26}" text-anchor="middle" '
             f'font-size="17" font-weight="600" fill="#222">{title}</text>')
    # faint guide spokes + rim labels
    for d in range(base):
        th = d * 2 * math.pi / base
        x = cx + (R + 4) * math.sin(th)
        y = cy - (R + 4) * math.cos(th)
        cop = coprime(d, base)
        s.append(f'<line x1="{cx}" y1="{cy}" x2="{cx+R*math.sin(th):.1f}" '
                 f'y2="{cy-R*math.cos(th):.1f}" stroke="{"#e8eaed" if cop else "#f4d7d7"}" '
                 f'stroke-width="1"/>')
        lab = "0123456789ABCDEF"[d] if base <= 16 else str(d)
        lx, ly = cx + (R + 16) * math.sin(th), cy - (R + 16) * math.cos(th)
        s.append(f'<text x="{lx:.1f}" y="{ly+4:.1f}" text-anchor="middle" '
                 f'font-size="10" fill="{"#888" if cop else "#d99"}">{lab}</text>')
    s.append(f'<circle cx="{cx}" cy="{cy}" r="{R}" fill="none" stroke="#ddd"/>')
    # the primes
    for p in primes:
        if p > M:
            break
        r = R * math.sqrt(p / M)          # area-uniform radius
        th = (p % base) * 2 * math.pi / base
        x = cx + r * math.sin(th)
        y = cy - r * math.cos(th)
        col = OCT_COLOR.get(p % 8, "#00798c") if color_by_oct else "#00798c"
        s.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.5" fill="{col}" opacity="0.85"/>')
    ncorr = sum(coprime(d, base) for d in range(base))
    s.append(f'<text x="{cx}" y="{cy+R+22}" text-anchor="middle" font-size="11.5" '
             f'fill="#555">{ncorr} of {base} spokes populated = '
             f'{100*ncorr/base:.1f}% admissible</text>')
    return "\n".join(s)


def fig_rays():
    M = 4096
    pr = primes_upto(M)
    R, pad, top = 150, 60, 122
    cell = 2 * R + pad
    W = 2 * cell + pad
    H = top + 2 * cell + 44
    # top row: the power-of-two family (evenly spaced rays, 50%)
    # bottom row: the primorial family (irregular rays, genuinely tighter)
    panels = [
        (8,  "OCTAL  (base 8 = 2³)"),
        (16, "HEX  (base 16 = 2⁴)"),
        (6,  "PRIMORIAL  (base 6 = 2·3)"),
        (30, "PRIMORIAL  (base 30 = 2·3·5)"),
    ]
    body = []
    for i, (b, t) in enumerate(panels):
        col, row = i % 2, i // 2
        cx = pad + R + col * cell
        cy = top + R + row * cell
        body.append(ray_panel(cx, cy, R, b, pr, M, t, color_by_oct=(b in (8, 16))))
    body.append(f'<text x="{W/2}" y="34" text-anchor="middle" font-size="18" '
                f'font-weight="700" fill="#222">Where primes land: the same primes on different wheels</text>')
    body.append(f'<text x="{W/2}" y="58" text-anchor="middle" font-size="12" fill="#555">'
                f'top row = powers of two (octal&#8594;hex is just 2&#215; finer); '
                f'bottom row = primorial bases (fewer, genuinely tighter corridors)</text>')
    cap = (f'<text x="{W/2}" y="{H-14}" text-anchor="middle" font-size="12" fill="#333">'
           f'Each dot is a prime &#8804; {M} on spoke (p mod base), radius growing with p. '
           f'Colour = octal corridor (p mod 8): note red = {{1,9}}, etc. split across two hex spokes.</text>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}" font-family="system-ui,Arial">'
           f'<rect width="{W}" height="{H}" fill="white"/>'
           + "\n".join(body) + cap + "</svg>")
    p = os.path.join(OUT, "fig-rays.svg")
    open(p, "w").write(svg)
    print("wrote", p)


# --------------------------------------------------------------------------- #
#  Figure 2: dyadic refinement tree -- 2 -> 4 -> 8 -> 16 -> 32 nested rings     #
# --------------------------------------------------------------------------- #
def fig_refinement():
    bases = [2, 4, 8, 16, 32]
    W, H = 560, 600
    cx, cy = W / 2, 300
    rings = []
    r0, dr = 55, 48
    for i, b in enumerate(bases):
        r = r0 + i * dr
        for d in range(b):
            th = d * 2 * math.pi / b
            x = cx + r * math.sin(th)
            y = cy - r * math.cos(th)
            cop = coprime(d, b)
            col = OCT_COLOR.get(d % 8, "#00798c") if (cop and b >= 8) else ("#00798c" if cop else GREY)
            rr = 5 if cop else 2.5
            rings.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rr}" '
                         f'fill="{col if cop else GREY}" stroke="#fff" stroke-width="0.7"/>')
        rings.append(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="#eee"/>')
        nc = sum(coprime(d, b) for d in range(b))
        rings.append(f'<text x="{cx+r*math.sin(0.30)+6:.1f}" y="{cy-r*math.cos(0.30):.1f}" '
                     f'font-size="11" fill="#666">b={b}: {nc}/{b} = 50%</text>')
    title = (f'<text x="{cx}" y="34" text-anchor="middle" font-size="16" font-weight="600" '
             f'fill="#222">Dyadic refinement: every 2^k wheel is the last one at 2&#215; resolution</text>')
    sub = (f'<text x="{cx}" y="{H-24}" text-anchor="middle" font-size="12" fill="#444">'
           f'Filled = prime corridor (odd spoke). Going out one ring doubles the spokes '
           f'but the corridor density stays pinned at 50%.</text>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}" font-family="system-ui,Arial">'
           f'<rect width="{W}" height="{H}" fill="white"/>'
           + title + "\n".join(rings) + sub + "</svg>")
    p = os.path.join(OUT, "fig-refinement.svg")
    open(p, "w").write(svg)
    print("wrote", p)


# --------------------------------------------------------------------------- #
#  Figure 3: Chebyshev deviation -- the one place octal and hex truly differ    #
# --------------------------------------------------------------------------- #
def load_counts():
    counts = {}
    path = os.path.join(RES, "digits.csv")
    if not os.path.exists(path):
        return counts
    with open(path) as fh:
        for row in csv.DictReader(fh):
            counts.setdefault(int(row["base"]), {})[int(row["residue"])] = int(row["count"])
    return counts


def dev_panel(cx, cy, w, h, base, counts, title):
    corr = [d for d in range(base) if coprime(d, base)]
    tot = sum(counts[d] for d in corr)
    even = tot / len(corr)                      # expected even share
    devs = [(d, counts[d] - even) for d in corr]
    mx = max(abs(v) for _, v in devs) or 1
    bw = w / len(corr)
    s = [f'<text x="{cx+w/2}" y="{cy-10}" text-anchor="middle" font-size="14" '
         f'font-weight="600" fill="#222">{title}</text>']
    mid = cy + h / 2
    s.append(f'<line x1="{cx}" y1="{mid}" x2="{cx+w}" y2="{mid}" stroke="#bbb"/>')
    s.append(f'<text x="{cx-6}" y="{mid+4}" text-anchor="end" font-size="9" '
             f'fill="#999">even share</text>')
    for i, (d, v) in enumerate(devs):
        bh = (v / mx) * (h / 2 - 6)
        x = cx + i * bw + bw * 0.18
        bwid = bw * 0.64
        qr = (base == 8 and d == 1) or (base == 16 and d in (1, 9))
        col = "#d1495b" if qr else "#00798c"
        y = mid - bh if bh >= 0 else mid
        s.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bwid:.1f}" height="{abs(bh):.1f}" '
                 f'fill="{col}"/>')
        lab = "0123456789ABCDEF"[d]
        s.append(f'<text x="{x+bwid/2:.1f}" y="{mid+h/2+14:.1f}" text-anchor="middle" '
                 f'font-size="10" fill="{"#d1495b" if qr else "#555"}">{lab}'
                 f'{"*" if qr else ""}</text>')
    return "\n".join(s)


def fig_chebyshev(counts):
    if 8 not in counts:
        print("no counts; skip chebyshev fig"); return
    W, H = 720, 400
    body = [dev_panel(60, 120, 250, 180, 8, counts[8], "Octal corridors"),
            dev_panel(410, 120, 250, 180, 16, counts[16], "Hex corridors")]
    title = (f'<text x="{W/2}" y="30" text-anchor="middle" font-size="16" font-weight="600" '
             f'fill="#222">Chebyshev bias: prime count minus the even share, per corridor</text>')
    sub = (f'<text x="{W/2}" y="{H-12}" text-anchor="middle" font-size="11.5" fill="#444">'
           f'Red* = quadratic-residue class (octal 1; hex 1 &amp; 9). It runs light. '
           f'Octal-1 = hex-1 + hex-9, so the deficit refines exactly.</text>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}" font-family="system-ui,Arial">'
           f'<rect width="{W}" height="{H}" fill="white"/>'
           + title + "\n".join(body) + sub + "</svg>")
    p = os.path.join(OUT, "fig-chebyshev.svg")
    open(p, "w").write(svg)
    print("wrote", p)


# --------------------------------------------------------------------------- #
#  ASCII wheels for the terminal                                               #
# --------------------------------------------------------------------------- #
def ascii_wheel(base, counts):
    print(f"\n  base-{base} wheel (last digit), prime corridors marked #, forbidden .")
    corr = [d for d in range(base) if coprime(d, base)]
    tot = sum(counts.get(d, 0) for d in corr) or 1
    for d in range(base):
        lab = "0123456789ABCDEF"[d] if base <= 16 else f"{d:2d}"
        if coprime(d, base):
            share = 100 * counts.get(d, 0) / tot
            bar = "#" * round(share / 2)
            print(f"   spoke {lab}: {bar:<26} {share:6.3f}%")
        else:
            print(f"   spoke {lab}: {'.'*2:<26} forbidden")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig_rays()
    fig_refinement()
    counts = load_counts()
    fig_chebyshev(counts)
    if counts:
        ascii_wheel(8, counts[8])
        ascii_wheel(16, counts[16])
