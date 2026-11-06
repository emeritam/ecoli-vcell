"""Step 4 - Download the chosen structures into structures/.

PDB entries: biological assembly 1 as mmCIF from RCSB (falls back to the
asymmetric unit if no assembly file exists). AlphaFold: the model URL is read
from the AlphaFold DB API, so new model versions are picked up automatically.
Files already on disk are skipped, so the step can be re-run after interruptions.

    python step4_download_structures.py --workers 6
"""
import argparse
from concurrent.futures import ThreadPoolExecutor

from common import DATA, STRUCT, cached_json, download, read_csv, session, write_csv


def fetch(sess, row):
    if row["source"] == "pdb" and row["pdb_id"]:
        pid = row["pdb_id"].lower()
        for name, url in ((f"{pid}-assembly1.cif.gz", f"https://files.rcsb.org/download/{pid}-assembly1.cif.gz"),
                          (f"{pid}.cif.gz", f"https://files.rcsb.org/download/{pid}.cif.gz")):
            if download(sess, url, STRUCT / name):
                return name, "assembly 1" if "assembly" in name else "asymmetric unit"
        return "", "PDB download failed"
    if row["source"] == "afdb":
        acc = row["members"]
        info = cached_json(sess, f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", f"afdb_{acc}")
        if not info:
            return "", "no AlphaFold entry"
        url = info[0].get("cifUrl") or info[0].get("pdbUrl")
        name = url.rsplit("/", 1)[-1]
        return (name, "AlphaFold model") if download(sess, url, STRUCT / name) else ("", "AlphaFold download failed")
    return "", "no structure"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()

    rows = read_csv(DATA / "structures.csv")
    sess = session()

    def job(row):
        try:
            return fetch(sess, row)
        except Exception as e:
            return "", f"error: {e}"

    with ThreadPoolExecutor(a.workers) as pool:
        for row, (fname, status) in zip(rows, pool.map(job, rows)):
            row["file"], row["file_status"] = fname, status
    write_csv(DATA / "structures.csv", rows)
    ok = sum(1 for r in rows if r["file"])
    print(f"{ok}/{len(rows)} structures available in structures/")
    for r in rows:
        if not r["file"] and r["source"] != "none":
            print(f"  missing: {r['ingredient_id']} ({r['file_status']})")


if __name__ == "__main__":
    main()
