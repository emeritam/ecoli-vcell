"""Step 4 - Download the chosen structures into structures/.

PDB entries: biological assembly 1 as mmCIF from RCSB (falls back to the
asymmetric unit if no assembly file exists). AlphaFold: the model URL is read
from the AlphaFold DB API, so new model versions are picked up automatically.
Files already on disk are skipped, so the step can be re-run after interruptions.
For PDB entries it also saves structures/<id>.entities.json, the chain-to-UniProt
map from the RCSB Data API (assembly files do not always include it).

    python step4_download_structures.py --workers 6
"""
import argparse
from concurrent.futures import ThreadPoolExecutor

import json

from common import DATA, STRUCT, cached_json, download, read_csv, session, write_csv


def entity_map(sess, pid):
    """entity_id -> UniProt accessions for a PDB entry (RCSB Data API), saved beside the structure."""
    dest = STRUCT / f"{pid}.entities.json"
    if dest.exists():
        return True
    entry = cached_json(sess, f"https://data.rcsb.org/rest/v1/core/entry/{pid.upper()}", f"rcsb_entry_{pid}")
    if not entry:
        return False
    eids = (entry.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids") or []
    out = {}
    for eid in eids:
        ent = cached_json(sess, f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid.upper()}/{eid}",
                          f"rcsb_entity_{pid}_{eid}") or {}
        ids = ent.get("rcsb_polymer_entity_container_identifiers") or {}
        accs = set(ids.get("uniprot_ids") or [])
        for ref in ids.get("reference_sequence_identifiers") or []:
            if str(ref.get("database_name", "")).lower() == "uniprot" and ref.get("database_accession"):
                accs.add(ref["database_accession"])
        poly = (ent.get("entity_poly") or {}).get("rcsb_entity_polymer_type", "")
        out[str(eid)] = {"uniprot": sorted(accs), "type": poly}
    dest.write_text(json.dumps(out))
    return True


def fetch(sess, row):
    if row["source"] == "pdb" and row["pdb_id"]:
        pid = row["pdb_id"].lower()
        for name, url in ((f"{pid}-assembly1.cif.gz", f"https://files.rcsb.org/download/{pid}-assembly1.cif.gz"),
                          (f"{pid}.cif.gz", f"https://files.rcsb.org/download/{pid}.cif.gz")):
            if download(sess, url, STRUCT / name):
                mapped = entity_map(sess, pid)
                kind = "assembly 1" if "assembly" in name else "asymmetric unit"
                return name, kind if mapped else kind + " (no UniProt map)"
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
