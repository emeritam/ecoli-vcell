"""Step 2 - Add copies per cell to every protein.

You supply the abundance file (these journals/databases don't offer a stable
download URL). Two formats are supported:

1. A table (.xlsx, .csv, .tsv) with a gene/ID column and a value column,
   e.g. Li et al. 2014 Cell Table S1 (synthesis rates in molecules per
   generation, which approximate copies per cell for stable proteins).
     python step2_add_abundance.py li2014_tableS1.xlsx --value-col "MOPS complete" --skiprows 1

2. A PaxDb abundance file for E. coli (species 511145), values in ppm:
     python step2_add_abundance.py 511145-integrated.txt --format paxdb --units ppm --total-proteins 3e6

IDs are matched by b-number (b0001...) when present, otherwise by gene name or
synonym. Unmatched IDs are written to data/abundance_unmatched.txt.
Run without --value-col first to see the available columns.
"""
import argparse
import re
import sys

from common import DATA, fnum, read_csv, write_csv


def build_index(proteins):
    by_b, by_name = {}, {}
    for p in proteins:
        if p["b_number"]:
            by_b[p["b_number"].lower()] = p
        by_name.setdefault(p["gene"].lower(), p)
    for p in proteins:                         # synonyms never override primary names
        for s in p["synonyms"].split():
            by_name.setdefault(s.lower(), p)
    return by_b, by_name


def match(ident, by_b, by_name):
    s = str(ident).strip()
    m = re.search(r"b\d{4}", s.lower())
    if m and m.group(0) in by_b:
        return by_b[m.group(0)]
    s = re.sub(r"^\d+\.", "", s)               # PaxDb-style "511145.xxx" prefix
    return by_name.get(s.lower())


def load_table(path, id_col, value_col, sheet, skiprows):
    import pandas as pd
    suffix = path.lower().rsplit(".", 1)[-1]
    if suffix in ("xlsx", "xls"):
        df = pd.read_excel(path, sheet_name=sheet if sheet is not None else 0, skiprows=skiprows)
    else:
        df = pd.read_csv(path, sep=None, engine="python", skiprows=skiprows)
    df.columns = [str(c).strip() for c in df.columns]
    print("Columns found:", " | ".join(df.columns))
    if id_col is None:
        cands = [c for c in df.columns if re.search(r"gene|locus|b.?number|\bid\b|name", c, re.I)]
        if not cands:
            sys.exit("Could not guess the gene/ID column; pass --id-col.")
        id_col = cands[0]
    if value_col is None:
        numeric = [c for c in df.columns if c != id_col and pd.api.types.is_numeric_dtype(df[c])]
        if len(numeric) != 1:
            sys.exit(f"Pick the abundance column with --value-col (numeric columns: {numeric}).")
        value_col = numeric[0]
    print(f"Using ID column '{id_col}' and value column '{value_col}'")
    return [(r[id_col], r[value_col]) for _, r in df.iterrows()], f"{path} [{value_col}]"


def load_paxdb(path):
    pairs = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.split()
            if len(parts) < 2:
                continue
            ident = next((t for t in parts[:-1] if re.search(r"b\d{4}", t)), parts[-2])
            pairs.append((ident, parts[-1]))
    return pairs, f"{path} [PaxDb]"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--format", choices=["table", "paxdb"], default="table")
    ap.add_argument("--id-col")
    ap.add_argument("--value-col")
    ap.add_argument("--sheet", help="Excel sheet name or index")
    ap.add_argument("--skiprows", type=int, default=0, help="header rows to skip before the column names")
    ap.add_argument("--units", choices=["copies", "ppm"], default="copies")
    ap.add_argument("--skip-bracketed", action="store_true",
                    help="ignore values written in brackets, e.g. [21] (lower-confidence values in Li et al. 2014)")
    ap.add_argument("--total-proteins", type=float, default=3.0e6,
                    help="protein molecules per cell, used to convert ppm (default 3e6)")
    a = ap.parse_args()

    proteins = read_csv(DATA / "proteome.csv")
    by_b, by_name = build_index(proteins)
    sheet = int(a.sheet) if a.sheet and a.sheet.isdigit() else a.sheet
    pairs, label = (load_paxdb(a.file) if a.format == "paxdb"
                    else load_table(a.file, a.id_col, a.value_col, sheet, a.skiprows))

    copies, unmatched, total_in, total_matched, n_brack = {}, [], 0.0, 0.0, 0
    for ident, val in pairs:
        sval = str(val).strip()
        bracketed = sval.startswith("[") and sval.endswith("]")
        if bracketed:
            if a.skip_bracketed:
                continue
            n_brack += 1
        v = fnum(sval.strip("[]"), None)
        if v is None or v != v or v < 0:      # skip blanks, NA and NaN
            continue
        if a.units == "ppm":
            v = v * 1e-6 * a.total_proteins
        # rows like "tufA+tufB": paralogs that could not be told apart; split evenly
        names = [n for n in re.split(r"\s*\+\s*", str(ident).strip()) if n]
        for name in names:
            share = v / len(names)
            total_in += share
            p = match(name, by_b, by_name)
            if p is None:
                unmatched.append(name)
                continue
            total_matched += share
            acc = p["accession"]
            copies[acc] = max(copies.get(acc, 0.0), share)   # duplicate rows: keep the larger value

    for p in proteins:
        c = copies.get(p["accession"])
        p["copies_per_cell"] = f"{c:.1f}" if c is not None else ""
        p["abundance_source"] = label if c is not None else ""
    write_csv(DATA / "proteins.csv", proteins)
    (DATA / "abundance_unmatched.txt").write_text("\n".join(unmatched), encoding="utf-8")

    have = [p for p in proteins if p["copies_per_cell"]]
    if n_brack:
        print(f"{n_brack} bracketed values were included (use --skip-bracketed to drop them)")
    print(f"Matched {len(have)} proteins ({len(unmatched)} IDs unmatched -> data/abundance_unmatched.txt)")
    if total_in:
        print(f"{100 * total_matched / total_in:.2f}% of the listed copies were matched")
    print(f"Total protein copies per cell: {sum(fnum(p['copies_per_cell']) for p in have):,.0f}")
    print("Most abundant:")
    for p in sorted(have, key=lambda x: -fnum(x["copies_per_cell"]))[:15]:
        print(f"  {p['gene']:<8} {fnum(p['copies_per_cell']):>10,.0f}  {p['protein_name'][:60]}")


if __name__ == "__main__":
    main()
