"""Download the sequence-79 half of the Volve ST0202 raw PZDETECTED file of shot line 5079 (field test I, raw data).
The file holds two shootings of line 5079: FFID 17... (cables 5113 / 5097) and FFID 18... (sequence 79, cables 5081 /
5065, shot points 2720 -> 1760). The FFID >= 181 000 000 part is contiguous at the end; it is copied with the 3600-byte
file header in front, by 8 MB range requests in 8 threads (a single 64 MB stream ran at 0.06 MB/s), resumable.
Token from the Databricks CLI at run time (volve_fetch.token), never stored.
usage: fetch_raw79.py <out.sgy>      (Git Bash: MSYS_NO_PATHCONV=1)"""
import sys, os, struct, time, threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0, "D:/ChatGPT/研究/地震衰减补偿研究/uq_rtm_20261117/scripts")
import volve_fetch as V  # noqa: E402

V.PATH = ("/Volumes/equinor_asa_volve_data_village/public/volve/Seismic/ST0202/Raw_data/ST0202_PZDETECTED_NAVMERGE_SHOT-SEQ-3-100/"
          "ST0202_PZDETECTED_NAVMERGE_SHOT-SEQ-3-100.00005079.sgy")
SIZE = 7144215312; TSZ = 7744; N = (SIZE - 3600) // TSZ; CH = 8 << 20; NTH = 8


def ffid(t):
    return struct.unpack(">i", V.rng(3600 + t * TSZ + 8, 3600 + t * TSZ + 11))[0]


def main(out):
    out = Path(out); meta = out.with_suffix(".first_trace")
    if meta.exists():
        t0 = int(meta.read_text())
    else:
        lo, hi = 0, N
        while lo < hi:
            mid = (lo + hi) // 2
            lo, hi = (mid + 1, hi) if ffid(mid) < 181_000_000 else (lo, mid)
        t0 = lo; meta.write_text(str(t0))
    a = 3600 + t0 * TSZ; total = 3600 + (N - t0) * TSZ
    print("first trace", t0, "ffid", ffid(t0), "previous", ffid(t0 - 1), "traces", N - t0, "bytes", total, flush=True)
    mode = "r+b" if out.exists() else "wb"
    if mode == "wb":
        with open(out, "wb") as f:
            f.write(V.rng(0, 3599)); f.truncate(total)
    done = out.with_suffix(".done8")              # start offsets (in the source file) of finished 8 MB pieces
    ok = set(int(x) for x in done.read_text().split()) if done.exists() else set()
    old = out.with_suffix(".chunks")              # first run: 64 MB pieces -> mark their 8 MB sub-pieces
    if old.exists():
        for k in (int(x) for x in old.read_text().split()):
            ok |= {a + k * (64 << 20) + i * CH for i in range(8)}
    todo = [s for s in range(a, SIZE, CH) if s not in ok]; lock = threading.Lock(); t = time.monotonic(); nb = [0]
    print("pieces to fetch", len(todo), "of", len(range(a, SIZE, CH)), flush=True)

    def job(s):
        e = min(s + CH, SIZE) - 1; d = V.rng(s, e)
        with lock:
            with open(out, "r+b") as f:
                f.seek(3600 + (s - a)); f.write(d); f.flush(); os.fsync(f.fileno())
            with open(done, "a") as g:
                g.write("%d\n" % s)
            nb[0] += len(d)
            print("piece at %d MB  %.2f MB/s" % ((s - a) >> 20, nb[0] / 2 ** 20 / (time.monotonic() - t)), flush=True)
    V.TOK = V.token()
    with ThreadPoolExecutor(NTH) as ex:
        list(ex.map(job, todo))
    print("done", out, os.path.getsize(out), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
