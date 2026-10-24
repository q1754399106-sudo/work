"""Download a verified uv wheel into this installation's private tools folder."""
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

root = Path(__file__).resolve().parent
tools = root / ".tools"
uv_dir = tools / "uv"
existing = list(uv_dir.rglob("uv.exe")) if uv_dir.exists() else []
if not existing:
    tools.mkdir(exist_ok=True)
    with urllib.request.urlopen("https://pypi.org/pypi/uv/0.13.5/json", timeout=60) as response:
        metadata = json.load(response)
    wheel = next(item for item in metadata["urls"] if item["filename"].endswith("win_amd64.whl"))
    with urllib.request.urlopen(wheel["url"], timeout=120) as response:
        payload = response.read()
    if hashlib.sha256(payload).hexdigest() != wheel["digests"]["sha256"]:
        raise RuntimeError("uv wheel SHA-256 verification failed")
    archive = tools / wheel["filename"]
    archive.write_bytes(payload)
    with zipfile.ZipFile(archive) as bundle:
        bundle.extractall(uv_dir)
    existing = list(uv_dir.rglob("uv.exe"))
print(existing[0])
