"""MATLAB-compatible randn for rng(seed,'twister') (mt19937ar generator, default Ziggurat transform).

Reconstructed 2026-10-25 from MATLAB-generated noise archived in this repository
(validate_randn.py). Agrees with MATLAB R2024b to <= 1e-14 absolute on every checked sample.

Per attempt two 32-bit MT19937 words a, b are drawn:
  strip j = b >> 24 (0 = top strip next to x = 0, 255 = base strip with the tail);
  U = ((a >> 3) * 2^24 + (b & 0xFFFFFF)) / 2^53, u = 2U - 1, candidate x = u * x_j;
  fast path |u| < x_{j-1}/x_j (256-strip Marsaglia-Tsang table, r = 3.654152885361009).
Wedge test and tail use MATLAB rand uniforms (genrand_res53 of two words):
  wedge: accept if f_j + U (f_{j-1} - f_j) < exp(-x^2/2), else draw a new attempt;
  tail:  repeat x = -log(U1)/r, y = -log(U2) until 2y >= x^2; return sign(u) (r + x).
rng(0) maps to MT seed 5489 as in MATLAB.
"""
import math, sys
import numpy as np

N = 256


def _table():
    f = lambda x: math.exp(-.5 * x * x)
    dn, vn = 3.654152885361009, 0.004928673233974658
    x = [vn / f(dn), dn]
    for _ in range(2, N):
        x.append(math.sqrt(-2 * math.log(vn / x[-1] + f(x[-1]))))
    return x + [0.0]


X = _table()                       # X[0] base (virtual) width, X[1] = r, X[256] = 0; index i = 255 - j
F = [math.exp(-.5 * v * v) for v in X]
R01 = X[1] / X[0]
RATIO = [X[i + 1] / X[i] for i in range(N)]


def words(seed, n):
    seed = 5489 if seed == 0 else seed
    return np.random.RandomState(seed).randint(0, 2**32, size=n, dtype=np.uint32).astype(np.int64).tolist()


def res53(a, b):
    return ((a >> 5) * 67108864.0 + (b >> 6)) / 9007199254740992.0


def randn(seed, n, offset=0):
    """n normals of the stream rng(seed,'twister') after `offset` words; returns (values, new_offset)."""
    w = words(seed, offset + int(n * 2.1) + 256)
    k, m, out = offset, 0, np.zeros(n)
    while m < n:
        if k + 16 > len(w):
            w = words(seed, 2 * len(w))
        a, b = w[k], w[k + 1]
        k += 2
        i = 255 - (b >> 24)
        u = 2 * (((a >> 3) * 16777216.0 + (b & 0xFFFFFF)) / 9007199254740992.0) - 1
        x = u * X[i]
        if i == 0:
            if abs(u) < R01:
                out[m] = x
                m += 1
                continue
            while True:
                if k + 4 > len(w):
                    w = words(seed, 2 * len(w))
                xx = -math.log(res53(w[k], w[k + 1])) / X[1]
                yy = -math.log(res53(w[k + 2], w[k + 3]))
                k += 4
                if yy + yy >= xx * xx:
                    break
            out[m] = X[1] + xx if u > 0 else -(X[1] + xx)
            m += 1
            continue
        if abs(u) < RATIO[i]:
            out[m] = x
            m += 1
            continue
        U = res53(w[k], w[k + 1])
        k += 2
        if F[i] + U * (F[i + 1] - F[i]) < math.exp(-.5 * x * x):
            out[m] = x
            m += 1
    return out, k


if __name__ == "__main__":
    # CLI for the Octave shim: mlrandn.py seed offset n outfile -> raw float64 LE, prints new offset
    seed, offset, n, path = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    v, k = randn(seed, n, offset)
    v.astype("<f8").tofile(path)
    print(k)
