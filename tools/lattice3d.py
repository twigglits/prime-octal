#!/usr/bin/env python3
"""Lift the prime wheel into 3D: digits on the planar axes, magnitude up the z-axis.

Two views, each for octal (base 8) and hex (base 16):

  (A) digit x order-of-magnitude RASTER  -- the literal construction in the brief:
      x = last digit (0..b-1), row z = floor(n/b) (each row = "the next block",
      i.e. the next order of magnitude step). Cell coloured by what kills n:
      div by 2 / 3 / 5 / 7, a larger-factor composite, or PRIME (bright).
      Reveals the prime corridors as vertical stripes and the mod-3/5/7
      divisibility lattices as diagonals -- whose SLOPE differs octal vs hex.

  (B) 3D POINT CLOUD (bijective digit layout):
      X = n mod b           (units digit, 0..b-1)
      Y = (n/b) mod b       (next digit, 0..b-1)
      Z = n / b^2           (everything above = the order of magnitude)
      One voxel per prime. Exported as .ply (open in MeshLab/Blender) and
      rendered to an isometric .png preview.

Pure stdlib: hand-rolled PNG encoder, software isometric rasteriser. No deps.
Run:  python3 tools/lattice3d.py [NMAX]
"""
import math, os, struct, zlib, sys

OUT = os.path.join(os.path.dirname(__file__), "..", "docs")

# ----------------------------------------------------------------- primes ---- #
def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i*i::i] = b"\x00" * len(s[i*i::i])
    return s

# ----------------------------------------------------------------- PNG ------- #
def write_png(path, w, h, buf):
    def chunk(typ, data):
        return (struct.pack(">I", len(data)) + typ + data +
                struct.pack(">I", zlib.crc32(typ + data) & 0xffffffff))
    raw = bytearray()
    stride = w * 3
    for y in range(h):
        raw.append(0)                       # filter: none
        raw += buf[y*stride:(y+1)*stride]
    png = (b"\x89PNG\r\n\x1a\n"
           + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
           + chunk(b"IDAT", zlib.compress(bytes(raw), 9))
           + chunk(b"IEND", b""))
    open(path, "wb").write(png)
    print("wrote", os.path.normpath(path))

def mk(w, h, bg=(255, 255, 255)):
    b = bytearray(bg * (w * h))
    return b

def px(buf, w, h, x, y, c):
    if 0 <= x < w and 0 <= y < h:
        i = (y * w + x) * 3
        buf[i:i+3] = bytes(c)

def rect(buf, w, h, x0, y0, dw, dh, c):
    for y in range(y0, y0 + dh):
        for x in range(x0, x0 + dw):
            px(buf, w, h, x, y, c)

def disc(buf, w, h, cx, cy, r, c):
    for y in range(cy - r, cy + r + 1):
        for x in range(cx - r, cx + r + 1):
            if (x - cx)**2 + (y - cy)**2 <= r*r:
                px(buf, w, h, x, y, c)

def line(buf, w, h, x0, y0, x1, y1, c):
    x0, y0, x1, y1 = int(x0), int(y0), int(x1), int(y1)
    dx, dy = abs(x1-x0), -abs(y1-y0)
    sx, sy = (1 if x0 < x1 else -1), (1 if y0 < y1 else -1)
    err = dx + dy
    while True:
        px(buf, w, h, x0, y0, c)
        if x0 == x1 and y0 == y1:
            break
        e2 = 2*err
        if e2 >= dy: err += dy; x0 += sx
        if e2 <= dx: err += dx; y0 += sy

# ------------------------------------------------------ (A) raster view ----- #
COL = {                       # smallest-factor categories
    "p2":   (226, 232, 240),  # even          -> pale slate
    "p3":   (120, 170, 220),  # /3            -> blue
    "p5":   (130, 200, 150),  # /5            -> green
    "p7":   (235, 200, 120),  # /7            -> amber
    "big":  (170, 175, 182),  # composite, smallest factor > 7
    "prime":(214,  40,  60),  # PRIME         -> bright red
}
def category(n, sv):
    if n < 2:
        return "p2"
    for p in (2, 3, 5, 7):
        if n % p == 0:
            return "prime" if n == p else {2:"p2",3:"p3",5:"p5",7:"p7"}[p]
    return "prime" if sv[n] else "big"

def raster(base, rows, sv, cell=18, gap=1):
    W, H = base * cell, rows * cell
    buf = mk(W, H)
    for z in range(rows):
        for d in range(base):
            n = base * z + d
            if n >= len(sv):
                continue
            cat = category(n, sv)
            c = COL[cat]
            x0, y0 = d * cell, (rows - 1 - z) * cell      # z grows upward
            rect(buf, W, H, x0 + gap, y0 + gap, cell - gap, cell - gap, c)
            if cat == "prime":                            # make primes pop
                rect(buf, W, H, x0 + gap + 2, y0 + gap + 2,
                     cell - gap - 4, cell - gap - 4, (150, 10, 30))
    name = {8: "octal", 16: "hex", 30: "primorial30"}[base]
    write_png(os.path.join(OUT, f"lattice-{name}.png"), W, H, buf)

# ------------------------------------------------------ (B) point cloud ----- #
# octal corridor colours (consistent with the wheel figures)
PAL = [(120,125,130),(214,40,60),(120,125,130),(237,174,73),
       (120,125,130),(0,121,140),(120,125,130),(102,161,130)]  # by d mod 8

def cloud(base, nmax, sv):
    b2 = base * base
    pts = []                                   # (X,Y,Z, color)
    for n in range(2, min(nmax, len(sv)-1) + 1):
        if sv[n]:
            X, Y, Z = n % base, (n // base) % base, n // b2
            pts.append((X, Y, Z, PAL[n % 8]))
    maxz = max(p[2] for p in pts)
    name = {8: "octal", 16: "hex"}[base]

    # ---- PLY export (spin it yourself) ----
    ply = [f"ply\nformat ascii 1.0\nelement vertex {len(pts)}",
           "property float x\nproperty float y\nproperty float z",
           "property uchar red\nproperty uchar green\nproperty uchar blue",
           "end_header"]
    for X, Y, Z, c in pts:
        ply.append(f"{X} {Y} {Z*0.5:.3f} {c[0]} {c[1]} {c[2]}")
    open(os.path.join(OUT, f"cloud-{name}.ply"), "w").write("\n".join(ply) + "\n")
    print("wrote", os.path.normpath(os.path.join(OUT, f"cloud-{name}.ply")))

    # ---- isometric render ----
    ca = math.cos(math.radians(30)); sa = math.sin(math.radians(30))
    zk = 0.62                                   # vertical units per z-level
    def proj(X, Y, Z):
        return ((X - Y) * ca, (X + Y) * sa - Z * zk)
    # bounding-box corners to fit canvas
    corners = [proj(x, y, z) for x in (0, base-1) for y in (0, base-1) for z in (0, maxz)]
    xs = [c[0] for c in corners]; ys = [c[1] for c in corners]
    scale = 1.0
    span_x = max(xs) - min(xs); span_y = max(ys) - min(ys)
    targ_h = 1150; scale = targ_h / span_y
    W = int(span_x * scale) + 80
    H = int(span_y * scale) + 80
    def to_screen(X, Y, Z):
        sx, sy = proj(X, Y, Z)
        return (int((sx - min(xs)) * scale) + 40,
                int((max(ys) - sy) * scale) + 40)
    buf = mk(W, H, (250, 250, 251))
    # box wireframe for 3D legibility
    box = (180, 184, 190)
    def edge(a, b): line(buf, W, H, *to_screen(*a), *to_screen(*b), box)
    c000=(0,0,0); c700=(base-1,0,0); c070=(0,base-1,0); c770=(base-1,base-1,0)
    c00z=(0,0,maxz); c70z=(base-1,0,maxz); c07z=(0,base-1,maxz); c77z=(base-1,base-1,maxz)
    for a,b in [(c000,c700),(c000,c070),(c700,c770),(c070,c770),
                (c00z,c70z),(c00z,c07z),(c70z,c77z),(c07z,c77z),
                (c000,c00z),(c700,c70z),(c070,c07z),(c770,c77z)]:
        edge(a,b)
    # points, painter's order (far -> near = small x+y+z first), depth shade
    pts.sort(key=lambda p: p[0]+p[1]+p[2])
    for X, Y, Z, c in pts:
        f = 0.55 + 0.45 * (Z / maxz)           # nearer top = brighter
        col = tuple(int(v*f + 255*(1-f)*0.15) for v in c)
        sx, sy = to_screen(X, Y, Z)
        disc(buf, W, H, sx, sy, 2, col)
    # axis labels
    write_png(os.path.join(OUT, f"cloud-{name}.png"), W, H, buf)
    return len(pts), maxz

# --------------------------------------------------- structural report ------ #
def slopes(base):
    """For each small prime p, the column killed in row z is x ≡ s*z (mod p);
    report s = (-base) mod p -- the diagonal slope of the 'divisible by p' lattice."""
    out = {}
    for p in (3, 5, 7):
        out[p] = (-base) % p
    return out

def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 16384
    os.makedirs(OUT, exist_ok=True)
    sv = sieve(max(nmax, 30 * 50))
    # (A) rasters: 105 rows = lcm(3,5,7) of the vertical periods, structure is clear.
    raster(8, 105, sv)
    raster(16, 105, sv)
    # bonus control: a PRIMORIAL width tiles cleanly (period 210/30 = 7 rows),
    # because 30 is divisible by 3 and 5 so those lattices become vertical, not sheared.
    raster(30, 42, sv, cell=14)
    # (B) clouds -- octal slab is b^2=64 wide, hex 256, so quarter the octal
    # n-range to give both the same z-height (proportionate boxes for comparison).
    no, zo = cloud(8, max(64, nmax // 4), sv)
    nh, zh = cloud(16, nmax, sv)

    print("\n=== repeating-structure report ===")
    print(f"raster width = base; row z holds [base*z, base*z+base-1]; z = order of magnitude")
    for base in (8, 16, 30):
        forb = [d for d in range(base) if math.gcd(d, base) != 1]
        s = slopes(base)
        print(f"\nbase {base}:")
        print(f"  forbidden columns (always empty): {forb}  -> {len(forb)} vertical empty stripes")
        print(f"  corridors (can hold primes):       {[d for d in range(base) if math.gcd(d,base)==1]}")
        print(f"  divisibility diagonals  x ≡ s·z (mod p):  "
              f"p=3 slope {s[3]}, p=5 slope {s[5]}, p=7 slope {s[7]}")
    print(f"\ncloud points (primes ≤ {nmax}):  octal {no} (z up to {zo}),  hex {nh} (z up to {zh})")
    print("octal slopes (3,5,7) =", tuple(slopes(8)[p] for p in (3,5,7)),
          " vs hex =", tuple(slopes(16)[p] for p in (3,5,7)))
    print("the SAME ×2 delta: since 16 = 2·8, every diagonal slope doubles (mod p):")
    print("  2·(1,2,6) mod (3,5,7) = (2,4,5) = the hex slopes. (mod 3 this flips the tilt: 2 ≡ -1.)")

if __name__ == "__main__":
    main()
