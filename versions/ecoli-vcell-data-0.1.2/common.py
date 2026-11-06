"""Shared helpers for the E. coli virtual-cell data pipeline."""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"          # tables produced by each step
CACHE = ROOT / "cache"        # cached API answers (safe to delete)
STRUCT = ROOT / "structures"  # downloaded mmCIF files
OUT = ROOT / "output"         # coarse-grained ingredients for packing / the viewer
for d in (DATA, CACHE, STRUCT, OUT):
    d.mkdir(exist_ok=True)

PROTEOME_ID = "UP000000625"   # E. coli K-12 MG1655 reference proteome
ORGANISM = "Escherichia coli K-12 MG1655"
COMPARTMENTS = ("OM", "PERI", "IM", "CYTO", "NUC", "EXT")


def session():
    """requests session that retries on rate limits and server errors."""
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    s = requests.Session()
    retry = Retry(total=6, backoff_factor=1.5,
                  status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=("GET",))
    s.mount("https://", HTTPAdapter(max_retries=retry, pool_maxsize=16))
    s.headers["User-Agent"] = "ecoli-virtual-cell-data/0.1 (academic research script)"
    return s


def safe_key(text):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def cached_json(sess, url, key):
    """GET url as JSON, caching the answer on disk. A 404 is cached as None."""
    f = CACHE / f"{safe_key(key)}.json"
    if f.exists():
        return json.loads(f.read_text())
    r = sess.get(url, timeout=60)
    if r.status_code == 404:
        data = None
    else:
        r.raise_for_status()
        data = r.json() if r.text.strip() else None
    f.write_text(json.dumps(data))
    return data


def download(sess, url, dest):
    """Stream url to dest. Returns True on success, False on 404."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return True
    r = sess.get(url, timeout=120, stream=True)
    if r.status_code == 404:
        return False
    r.raise_for_status()
    tmp = dest.with_suffix(dest.suffix + ".part")
    with open(tmp, "wb") as fh:
        for chunk in r.iter_content(1 << 16):
            fh.write(chunk)
    tmp.rename(dest)
    return True


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def fnum(x, default=0.0):
    try:
        return float(str(x).replace(",", ""))
    except (TypeError, ValueError):
        return default
