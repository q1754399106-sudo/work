"""Read-only access to the Volve ST0202 PreMig PZ shot gathers by HTTP range requests (stage 4).

The access token is read at run time from the Databricks CLI (`databricks auth token`), never stored.
usage:
  volve_fetch.py bounds <seq>                       first and last trace index of a shot sequence
  volve_fetch.py layout <first_trace> <n>           header summary of n consecutive traces
  volve_fetch.py fetch <first> <last> <dline> <out> SEG-Y with the traces of detector line <dline>
Run from Git Bash with MSYS_NO_PATHCONV=1 (the volume path starts with '/').
"""
import sys, os, json, struct, subprocess, time, glob, urllib.request, urllib.parse

HOST = "https://dbc-f8011bce-eef5.cloud.databricks.com"
PATH = "/Volumes/equinor_asa_volve_data_village/public/volve/Seismic/ST0202/PreMig_data/ST0202R08_PZ_Shot_gathers_in_PP_Time.segy"
NS, TSZ, NTR = 1876, 7744, 24963410
HDR = dict(fldr=(8, "i"), tracf=(12, "i"), offset=(36, "i"), sx=(72, "i"), sy=(76, "i"), gx=(80, "i"), gy=(84, "i"),
           seq=(180, "i"), sline=(184, "i"), dline=(188, "i"), spoint=(192, "i"), dpoint=(196, "i"))


def token():
    exe = glob.glob(os.path.expanduser("~/AppData/Local/Microsoft/WinGet/Packages/Databricks.DatabricksCLI_*/databricks.exe"))[0]
    return json.loads(subprocess.run([exe, "auth", "token"], capture_output=True, text=True).stdout)["access_token"]


TOK = None


def rng(a, b):
    global TOK
    TOK = TOK or token()
    url = HOST + "/api/2.0/fs/files" + urllib.parse.quote(PATH)
    for i in range(6):
        try:
            r = urllib.request.Request(url, headers={"Authorization": "Bearer " + TOK, "Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(r, timeout=300) as f:
                d = f.read()
            if len(d) == b - a + 1:
                return d
        except Exception as e:  # retry with backoff, refresh the token once
            if i == 2:
                TOK = token()
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"range {a}-{b} failed")


def hdr(t, raw=None):
    h = raw if raw is not None else rng(3600 + t * TSZ, 3600 + t * TSZ + 239)
    return {k: struct.unpack(">" + c, h[o:o + 4])[0] for k, (o, c) in HDR.items()}


def bounds(seq):
    def first_ge(s):  # first trace with sequence number >= s (sequences are stored in increasing order)
        lo, hi = 0, NTR
        while lo < hi:
            mid = (lo + hi) // 2
            if hdr(mid)["seq"] < s:
                lo = mid + 1
            else:
                hi = mid
        return lo
    a = first_ge(seq); b = first_ge(seq + 1) - 1
    print(json.dumps(dict(seq=seq, first=a, last=b, n=b - a + 1, first_hdr=hdr(a), last_hdr=hdr(b))))


def layout(t0, n):
    d = rng(3600 + t0 * TSZ, 3600 + (t0 + n) * TSZ - 1); prev = None; run = 0
    for i in range(n):
        h = hdr(0, d[i * TSZ:i * TSZ + 240]); key = (h["fldr"], h["dline"])
        if key != prev:
            if prev is not None:
                print(f"  fldr {prev[0]} dline {prev[1]}: {run} traces")
            prev, run = key, 0
        run += 1
    print(f"  fldr {prev[0]} dline {prev[1]}: {run} traces")


def fetch(t0, t1, dline, out):
    """Scan the sequence shot by shot: read one header per block to find the runs of <dline>, then read only those."""
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    head = rng(0, 3599); n_out = 0; t = t0; tmp = out + ".part"
    done = 0; last_len = {}
    if os.path.exists(tmp):  # resume
        done = (os.path.getsize(tmp) - 3600) // TSZ
    with open(tmp, "ab") as fo:
        if done == 0:
            fo.write(head)
        while t <= t1:
            h = hdr(t); key = (h["fldr"], h["dline"])
            same = lambda u: u <= t1 and (lambda g: (g["fldr"], g["dline"]) == key)(hdr(u))
            # predict the run length from the previous run of the same cable, check 2 headers, search nearby
            L = last_len.get(h["dline"], 80); end = None
            for d in [0, 1, -1, 2, -2, 3, -3, 4, -4, 6, -6, 8, -8, 12, -12]:
                e = t + L - 1 + d
                if e >= t and same(e) and not same(e + 1):
                    end = e; break
            if end is None:  # fall back to galloping + bisection
                step = 1
                while same(t + step):
                    step *= 2
                lo, hi = t + step // 2, min(t + step, t1 + 1)
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    lo, hi = (mid, hi) if same(mid) else (lo, mid)
                end = lo if step > 1 else t
            last_len[h["dline"]] = end - t + 1
            if h["dline"] == dline:
                if n_out >= done:
                    fo.write(rng(3600 + t * TSZ, 3600 + (end + 1) * TSZ - 1)); fo.flush()
                n_out += end - t + 1
                if n_out % 2400 < (end - t + 1):
                    print(f"  shot fldr {h['fldr']}: {n_out} traces of dline {dline}", flush=True)
            t = end + 1
    os.replace(tmp, out); print(f"done {out}: {n_out} traces")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "bounds":
        bounds(int(sys.argv[2]))
    elif cmd == "layout":
        layout(int(sys.argv[2]), int(sys.argv[3]))
    elif cmd == "fetch":
        fetch(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
