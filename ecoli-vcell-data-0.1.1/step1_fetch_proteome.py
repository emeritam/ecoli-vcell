"""Step 1 - Download the E. coli K-12 MG1655 reference proteome (UP000000625) from UniProt.

Output: data/proteome.csv (one row per protein), data/sequences.json and
data/proteome_raw.tsv. Transmembrane segments and topological domains are kept
so step 5 can orient membrane proteins.
Each protein gets a compartment (OM, IM, PERI, CYTO, NUC, EXT) from UniProt's
subcellular-location annotation, falling back to transmembrane segments.

    python step1_fetch_proteome.py               # all entries in the proteome
    python step1_fetch_proteome.py --reviewed    # Swiss-Prot entries only
"""
import argparse
import json
import re

from common import DATA, PROTEOME_ID, session, write_csv

FIELDS = ["accession", "gene_primary", "gene_names", "gene_oln", "protein_name",
          "length", "mass", "cc_subcellular_location", "ft_transmem",
          "xref_pdb", "xref_alphafolddb", "cc_subunit", "ft_topo_dom", "sequence"]

# Checked in this order; the first match wins. Membranes come before
# "periplasm" so that "inner membrane; periplasmic side" stays IM.
LOCATION_RULES = [
    ("cell outer membrane", "OM"),
    ("cell inner membrane", "IM"),
    ("cell membrane", "IM"),
    ("periplasm", "PERI"),
    ("nucleoid", "NUC"),
    ("fimbrium", "EXT"),
    ("flagellum", "EXT"),
    ("secreted", "EXT"),
    ("cytoplasm", "CYTO"),
]


def assign_compartment(location, tm_count):
    text = location.lower()
    for key, comp in LOCATION_RULES:
        if key in text:
            return comp, "uniprot"
    if tm_count > 0:
        return "IM", "transmembrane segments"
    return "CYTO", "default (no annotation)"


def membrane_mode(location, tm_count, comp):
    """How a protein in a membrane compartment is attached, from UniProt keywords.

    Returns (compartment, mode, side). mode: integral, lipid_anchor or peripheral.
    side: the compartment the protein's body faces (cytoplasm, periplasm, outside)."""
    if comp not in ("OM", "IM"):
        return comp, "", ""
    l = location.lower()
    if tm_count > 0 or "multi-pass" in l or "single-pass" in l:
        return comp, "integral", ""
    if "lipid-anchor" in l:
        if "periplasmic side" in l:            # explicit side wins over notes on minor surface-exposed forms
            side = "periplasm"
        elif "cytoplasmic side" in l:
            side = "cytoplasm"
        elif "extracellular side" in l or "cell surface" in l:
            side = "outside"
        else:
            side = "periplasm"
        return comp, "lipid_anchor", side
    if "peripheral membrane protein" in l:
        if "cytoplasm" in l.replace("cytoplasmic side", ""):
            return "CYTO", "", ""             # mainly soluble, also seen at the membrane
        if "periplasmic side" in l:
            side = "periplasm"
        elif "cytoplasmic side" in l:
            side = "cytoplasm"
        else:
            side = "cytoplasm" if comp == "IM" else "periplasm"
        return comp, "peripheral", side
    return comp, "integral", ""                 # e.g. outer-membrane beta-barrels without strand annotation


def parse_features(text, kind):
    """'TRANSMEM 10..30; /note="Helical"; ...' -> [(10, 30, 'Helical'), ...]"""
    out = []
    for m in re.finditer(kind + r' <?(\d+)\.\.>?(\d+)(;\s*/note="([^"]*)")?', text):
        out.append((int(m.group(1)), int(m.group(2)), m.group(4) or ""))
    return out


def parse_tsv(text):
    lines = [ln for ln in text.splitlines() if ln.strip()]
    rows = []
    for ln in lines[1:]:                      # first line is the header
        parts = ln.split("\t")
        parts += [""] * (len(FIELDS) - len(parts))
        r = dict(zip(FIELDS, parts))
        names = r["gene_names"].split()
        gene = r["gene_primary"].split(";")[0].strip() or (names[0] if names else r["accession"])
        bnum = re.search(r"\bb\d{4}\b", r["gene_oln"])
        tm_segs = parse_features(r["ft_transmem"], "TRANSMEM")
        topo = parse_features(r["ft_topo_dom"], "TOPO_DOM")
        tm = max(len(tm_segs), r["ft_transmem"].count("TRANSMEM"))
        comp, basis = assign_compartment(r["cc_subcellular_location"], tm)
        comp2, mmode, mside = membrane_mode(r["cc_subcellular_location"], tm, comp)
        if comp2 != comp:
            comp, basis = comp2, "uniprot (peripheral, mainly cytoplasmic)"
        pdbs = [p for p in r["xref_pdb"].split(";") if p.strip()]
        rows.append({
            "accession": r["accession"],
            "gene": gene,
            "b_number": bnum.group(0) if bnum else "",
            "synonyms": " ".join(n for n in names if n != gene),
            "protein_name": r["protein_name"],
            "length": r["length"],
            "mass_da": r["mass"].replace(",", ""),
            "compartment": comp,
            "compartment_basis": basis,
            "tm_count": tm,
            "n_pdb": len(pdbs),
            "has_alphafold": "yes" if r["xref_alphafolddb"].strip() else "no",
            "subcellular_location": r["cc_subcellular_location"].replace("SUBCELLULAR LOCATION: ", ""),
            "subunit": r["cc_subunit"].replace("SUBUNIT: ", ""),
            "membrane_mode": mmode,
            "membrane_side": mside,
            "tm_segments": ",".join(f"{a}-{b}" for a, b, _ in tm_segs),
            "tm_type": "beta" if any("beta" in n.lower() for *_, n in tm_segs) else ("helical" if tm_segs else ""),
            "topology": ",".join(f"{a}-{b}:{n}" for a, b, n in topo if n),
            "_sequence": r["sequence"].strip(),
        })
    return rows


def apply_overrides(rows):
    """data/membrane_overrides.csv: gene,compartment,membrane_mode,membrane_side (blank = keep)."""
    f = DATA / "membrane_overrides.csv"
    if not f.exists():
        return
    import csv
    by_gene = {r["gene"].lower(): r for r in rows}
    for o in csv.DictReader(open(f, encoding="utf-8")):
        r = by_gene.get(o.get("gene", "").strip().lower())
        if not r:
            print(f"  override: gene {o.get('gene')} not found")
            continue
        for col in ("compartment", "membrane_mode", "membrane_side"):
            if (o.get(col) or "").strip():
                r[col] = o[col].strip()
        r["compartment_basis"] = "data/membrane_overrides.csv"
        print(f"  override: {r['gene']} -> {r['compartment']} {r['membrane_mode']} {r['membrane_side']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reviewed", action="store_true", help="keep Swiss-Prot (reviewed) entries only")
    a = ap.parse_args()

    query = f"(proteome:{PROTEOME_ID})" + (" AND (reviewed:true)" if a.reviewed else "")
    print(f"Querying UniProt: {query}")
    r = session().get("https://rest.uniprot.org/uniprotkb/stream",
                      params={"query": query, "format": "tsv", "fields": ",".join(FIELDS)},
                      timeout=300)
    r.raise_for_status()
    (DATA / "proteome_raw.tsv").write_text(r.text, encoding="utf-8")

    rows = parse_tsv(r.text)
    apply_overrides(rows)
    seqs = {row["accession"]: row.pop("_sequence") for row in rows}
    (DATA / "sequences.json").write_text(json.dumps(seqs))
    write_csv(DATA / "proteome.csv", rows)
    counts = {}
    for row in rows:
        counts[row["compartment"]] = counts.get(row["compartment"], 0) + 1
    print(f"{len(rows)} proteins -> data/proteome.csv")
    print("By compartment:", ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    print(f"With PDB entries: {sum(1 for x in rows if x['n_pdb'])}; "
          f"with AlphaFold models: {sum(1 for x in rows if x['has_alphafold'] == 'yes')}")


if __name__ == "__main__":
    main()
