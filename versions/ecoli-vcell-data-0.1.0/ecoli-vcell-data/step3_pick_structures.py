"""Step 3 - Choose one structure for each "ingredient" of the cell.

Ingredients are either
  * curated assemblies from data/assemblies.csv (ribosome, ATP synthase, RNA
    polymerase, GroEL/ES ...), which use one PDB entry for several genes, or
  * single proteins, which use
      1. the best experimental PDB entry (PDBe "best structures" ranking, which
         weighs coverage of the UniProt sequence and resolution), if it covers
         at least --min-coverage of the sequence, otherwise
      2. the AlphaFold DB model, otherwise
      3. a low-coverage PDB entry, if that is all there is.

Output: data/structures.csv

    python step3_pick_structures.py                 # every protein with a copy number
    python step3_pick_structures.py --top 500       # the 500 proteins with the most mass in the cell
"""
import argparse
import re
from concurrent.futures import ThreadPoolExecutor

from common import DATA, cached_json, fnum, read_csv, session, write_csv

PDBE_URLS = ["https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{acc}",
             "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/{acc}"]


def pdbe_best(sess, acc):
    """List of PDB chains mapped to a UniProt accession, best first."""
    for i, url in enumerate(PDBE_URLS):
        try:
            d = cached_json(sess, url.format(acc=acc), f"pdbe_best_{i}_{acc}")
        except Exception as e:                  # try the other endpoint
            print(f"  PDBe lookup failed for {acc} ({e}); trying fallback")
            continue
        if isinstance(d, dict):
            return d.get(acc, [])
        if d is None:
            return []
    return []


def resolve_members(spec, by_gene):
    accs, missing = [], []
    for tok in [t.strip() for t in spec.split(";") if t.strip()]:
        if tok.startswith("regex:"):
            rx = re.compile(tok[6:])
            hits = [p for g, p in by_gene.items() if rx.search(g)]
            accs += [p["accession"] for p in hits]
        elif tok.lower() in by_gene:
            accs.append(by_gene[tok.lower()]["accession"])
        else:
            missing.append(tok)
    return list(dict.fromkeys(accs)), missing


def describe(e):
    return {"pdb_id": e.get("pdb_id", "").lower(), "chain_id": e.get("chain_id", ""),
            "coverage": f"{fnum(e.get('coverage')):.2f}",
            "resolution": e.get("resolution") or "",
            "method": e.get("experimental_method", "")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assemblies", default=str(DATA / "assemblies.csv"))
    ap.add_argument("--min-copies", type=float, default=1.0)
    ap.add_argument("--top", type=int, default=0, help="keep only the N proteins with the most mass (copies x mass)")
    ap.add_argument("--min-coverage", type=float, default=0.6)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()

    proteins = read_csv(DATA / "proteins.csv")
    by_acc = {p["accession"]: p for p in proteins}
    by_gene = {}
    for p in proteins:
        by_gene.setdefault(p["gene"].lower(), p)
    sess = session()
    rows, consumed = [], set()

    # 1. curated assemblies
    assemblies = read_csv(a.assemblies)
    anchors = {}
    for asm in assemblies:
        accs, missing = resolve_members(asm["members"], by_gene)
        if missing:
            print(f"  {asm['id']}: genes not in proteome: {', '.join(missing)}")
        if not accs:
            print(f"  {asm['id']}: skipped (no members found)")
            continue
        consumed.update(accs)
        row = {"ingredient_id": asm["id"], "kind": "assembly", "name": asm["name"],
               "compartment": asm["compartment"], "mode": asm.get("mode") or "complex",
               "members": ";".join(accs),
               "member_genes": ";".join(by_acc[x]["gene"] for x in accs),
               "source": "pdb", "pdb_id": asm["pdb_id"].lower(), "chain_id": "",
               "coverage": "", "resolution": "", "method": "", "note": asm.get("note", "")}
        if row["pdb_id"] == "auto":
            anchors[asm["id"]] = accs[0]
        rows.append(row)

    # 2. single proteins
    cands = [p for p in proteins if p["accession"] not in consumed and fnum(p["copies_per_cell"]) >= a.min_copies]
    cands.sort(key=lambda p: -fnum(p["copies_per_cell"]) * fnum(p["mass_da"]))
    if a.top:
        cands = cands[:a.top]
    lookups = list(dict.fromkeys([p["accession"] for p in cands] + list(anchors.values())))
    print(f"Looking up PDBe best structures for {len(lookups)} proteins...")
    with ThreadPoolExecutor(a.workers) as pool:
        best = dict(zip(lookups, pool.map(lambda x: pdbe_best(sess, x), lookups)))

    for row in rows:                          # fill in "auto" assemblies
        if row["ingredient_id"] in anchors:
            entries = best.get(anchors[row["ingredient_id"]]) or []
            if entries:
                row.update(describe(entries[0]))
                row["chain_id"] = ""          # keep every chain of the assembly
            else:
                row.update(source="none", pdb_id="", note=(row["note"] + "; no PDB entry for anchor").strip("; "))

    stats = {"pdb": 0, "afdb": 0, "pdb_low_coverage": 0, "none": 0}
    for p in cands:
        entries = best.get(p["accession"]) or []
        good = [e for e in entries if fnum(e.get("coverage")) >= a.min_coverage]
        row = {"ingredient_id": p["gene"], "kind": "protein", "name": p["protein_name"],
               "compartment": p["compartment"], "mode": "", "members": p["accession"],
               "member_genes": p["gene"], "source": "", "pdb_id": "", "chain_id": "",
               "coverage": "", "resolution": "", "method": "", "note": ""}
        if good:
            row.update(source="pdb", **describe(good[0]))
            stats["pdb"] += 1
        elif p["has_alphafold"] == "yes":
            row.update(source="afdb")
            if entries:
                row["note"] = f"best PDB {entries[0].get('pdb_id')} covers only {fnum(entries[0].get('coverage')):.0%}"
            stats["afdb"] += 1
        elif entries:
            row.update(source="pdb", note="low coverage", **describe(entries[0]))
            stats["pdb_low_coverage"] += 1
        else:
            row.update(source="none")
            stats["none"] += 1
        rows.append(row)

    write_csv(DATA / "structures.csv", rows)
    print(f"{len(rows)} ingredients -> data/structures.csv "
          f"({sum(1 for r in rows if r['kind'] == 'assembly')} curated assemblies)")
    print("Single proteins:", ", ".join(f"{k} {v}" for k, v in stats.items()))


if __name__ == "__main__":
    main()
