"""Step 5 - Turn each structure into a packing ingredient.

For every ingredient this step
  * keeps the right chains (all chains for curated assemblies; only the chains
    of the target protein for single proteins, so partners in a co-crystal are
    not counted twice),
  * trims AlphaFold residues below --min-plddt (floppy tails would otherwise
    inflate the shape),
  * clusters heavy atoms with k-means into a set of spheres (about one sphere
    per --res-per-sphere residues), the multi-sphere representation autoPACK
    and cellPACK use for collision tests,
  * reads the stoichiometry of each member protein from the assembly file and
    turns protein copy numbers into copies of the ingredient per cell.

Output (lengths in angstroms):
  output/ingredients.json        everything, for packing and the viewer
  output/ingredients_summary.csv one line per ingredient
  output/ingredients/<id>.json   one file per ingredient

    python step5_coarse_grain.py --res-per-sphere 12 --workers 4
"""
import argparse
import datetime
import difflib
import json
import statistics
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from common import DATA, ORGANISM, OUT, PROTEOME_ID, STRUCT, fnum, read_csv, write_csv

VDW = 1.7  # angstrom, added to cluster radii


# ---------- structure reading ----------
def entity_accessions(block):
    """entity_id -> set of UniProt accessions, from _struct_ref."""
    out = {}
    tab = block.find("_struct_ref.", ["entity_id", "?db_name", "?pdbx_db_accession"])
    for r in tab:
        db = r.str(1) if r.has(1) else ""
        acc = r.str(2) if r.has(2) else ""
        if acc and (not db or db.upper().startswith("UNP")):
            out.setdefault(r.str(0), set()).add(acc)
    return out


def sidecar_map(path):
    """entity_id -> set of accessions from structures/<pdb>.entities.json (made by step 4)."""
    pid = path.name.split("-assembly")[0].split(".")[0].lower()
    f = path.parent / f"{pid}.entities.json"
    if not f.exists():
        return {}, {}
    d = json.loads(f.read_text())
    return ({k: set(v.get("uniprot", [])) for k, v in d.items()},
            {k: v.get("type", "") for k, v in d.items()})


def load(path):
    import gemmi
    doc = gemmi.cif.read(str(path))
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    acc_map = entity_accessions(block)
    side, types = sidecar_map(path)
    for k, v in side.items():
        if v:
            acc_map.setdefault(k, set()).update(v)
    for ent in st.entities:                    # fallback: SIFTS annotation
        sifts = getattr(ent, "sifts_unp_acc", None)
        if sifts and ent.name not in acc_map:
            acc_map[ent.name] = set(sifts)
    chains = []
    for ch in st[0]:
        span = ch.get_polymer()
        if not span or len(span) == 0:
            continue
        ent = st.get_entity_of(span)
        ptype = types.get(ent.name, "") if ent else ""
        if not ptype and ent is not None:
            ptype = "RNA" if "Rna" in str(ent.polymer_type) else ("DNA" if "Dna" in str(ent.polymer_type) else "Protein")
        chains.append({"name": ch.name, "entity": ent.name if ent else "?", "type": ptype,
                       "accs": acc_map.get(ent.name, set()) if ent else set(), "span": span})
    return chains


def one_letter(name):
    import gemmi
    t = gemmi.find_tabulated_residue(name)
    code = t.one_letter_code.upper() if t else "X"
    return code if code.isalpha() and code != " " else "X"


def parse_ranges(text, with_label=False):
    out = []
    for part in [p for p in str(text or "").split(",") if p]:
        rng_, _, label = part.partition(":")
        a, _, b = rng_.partition("-")
        try:
            out.append((int(a), int(b), label) if with_label else (int(a), int(b)))
        except ValueError:
            pass
    return out


def membrane_frame(keep, members, comp, kind, center):
    """Membrane normal and centre from UniProt TM segments mapped onto the structure.

    The normal points out of the cell: towards the periplasm for inner-membrane
    proteins, towards the outside for outer-membrane proteins."""
    by_acc = {m["accession"]: m for m in members}
    seg_vecs, tm_pts, inside, outside = [], [], [], []
    for c in keep:
        acc = next((a for a in c["accs"] if a in by_acc), None)
        if acc is None and (kind == "protein" or not c["accs"]):
            acc = members[0]["accession"] if kind == "protein" else None
        info = by_acc.get(acc)
        if not info or not info.get("seq") or not info.get("tm"):
            continue
        seq, cas = [], []
        for res in c["span"]:
            seq.append(one_letter(res.name))
            ca = res.find_atom("CA", "*")
            cas.append((ca.pos.x, ca.pos.y, ca.pos.z) if ca else None)
        uni_pos = {}
        sm = difflib.SequenceMatcher(None, "".join(seq), info["seq"], autojunk=False)
        for i, j, n in sm.get_matching_blocks():
            if n >= 5:
                for k in range(n):
                    uni_pos[j + k + 1] = i + k          # UniProt position (1-based) -> residue index
        def pts(a, b):
            return [cas[uni_pos[p]] for p in range(a, b + 1) if p in uni_pos and cas[uni_pos[p]] is not None]
        for a, b in info["tm"]:
            P = np.array(pts(a, b))
            if len(P) < 5:
                continue
            Pc = P - P.mean(0)
            v = np.linalg.svd(Pc, full_matrices=False)[2][0]
            if np.dot(v, P[-1] - P[0]) < 0:
                v = -v
            seg_vecs.append(v)
            tm_pts.extend(P.tolist())
        for a, b, label in info.get("topo", []):
            lab = label.lower()
            side = None
            if "cytoplasm" in lab:
                side = "in"
            elif "periplasm" in lab:
                side = "out" if comp == "IM" else "in"
            elif "extracellular" in lab or "exterior" in lab:
                side = "out"
            if side:
                (inside if side == "in" else outside).extend(pts(a, b))
    if not seg_vecs:
        return None
    ref = seg_vecs[0]
    normal = np.sum([v if np.dot(v, ref) >= 0 else -v for v in seg_vecs], axis=0)
    normal /= np.linalg.norm(normal)
    mc = np.mean(tm_pts, axis=0)
    known = True
    if inside and outside:
        d = np.dot(np.mean(outside, 0) - np.mean(inside, 0), normal)
    elif inside:
        d = -np.dot(np.mean(inside, 0) - mc, normal)
    elif outside:
        d = np.dot(np.mean(outside, 0) - mc, normal)
    else:
        d, known = 1.0, False
    if d < 0:
        normal = -normal
    return {"normal": [round(float(v), 4) for v in normal],
            "center_A": [round(float(v), 2) for v in (mc - center)],
            "outward_known": known, "n_segments": len(seg_vecs),
            "method": "UniProt transmembrane segments" + (" + topology" if known else "")}


def coords_of(chains, min_plddt=None):
    xyz, nres = [], 0
    for c in chains:
        for res in c["span"]:
            if min_plddt is not None and len(res) and res[0].b_iso < min_plddt:
                continue
            nres += 1
            xyz.extend((a.pos.x, a.pos.y, a.pos.z) for a in res)
    return np.asarray(xyz, dtype=np.float64).reshape(-1, 3), nres


# ---------- shape ----------
def _assign(X, C):
    lab = np.empty(len(X), dtype=np.int64)
    for i in range(0, len(X), 20000):
        d = ((X[i:i + 20000, None, :] - C[None, :, :]) ** 2).sum(-1)
        lab[i:i + 20000] = d.argmin(1)
    return lab


def kmeans(X, k, seed=0, iters=30):
    rng = np.random.default_rng(seed)
    S = X if len(X) <= 30000 else X[rng.choice(len(X), 30000, replace=False)]
    k = max(1, min(k, len(S)))
    C = [S[rng.integers(len(S))]]
    d2 = ((S - C[0]) ** 2).sum(1)
    for _ in range(1, k):                      # k-means++ seeding
        tot = d2.sum()
        i = rng.choice(len(S), p=d2 / tot) if tot > 0 else rng.integers(len(S))
        C.append(S[i])
        d2 = np.minimum(d2, ((S - S[i]) ** 2).sum(1))
    C = np.array(C)
    for _ in range(iters):
        lab = _assign(S, C)
        new = np.array([S[lab == j].mean(0) if np.any(lab == j) else C[j] for j in range(k)])
        if np.allclose(new, C, atol=1e-2):
            break
        C = new
    return C


def spheres_from(X, n_res, res_per_sphere, max_spheres):
    center = X.mean(0)
    Xc = X - center
    k = int(np.clip(round(n_res / res_per_sphere), 1, max_spheres))
    C = kmeans(Xc, k)
    lab = _assign(Xc, C)
    sph = []
    for j in range(len(C)):
        pts = Xc[lab == j]
        if len(pts) == 0:
            continue
        r = float(np.percentile(np.linalg.norm(pts - C[j], axis=1), 95)) + VDW
        sph.append([round(float(v), 2) for v in C[j]] + [round(r, 2)])
    dist = np.linalg.norm(Xc, axis=1)
    evals, evecs = np.linalg.eigh(np.cov(Xc.T)) if len(Xc) > 3 else (None, np.eye(3))
    return {
        "spheres_A": sph,
        "encapsulating_radius_A": round(float(dist.max()) + VDW, 2),
        "radius_of_gyration_A": round(float(np.sqrt((dist ** 2).mean())), 2),
        "principal_vector": [round(float(v), 4) for v in evecs[:, -1]],
    }


# ---------- one ingredient ----------
def process(job):
    row, members, opts = job
    out = {"id": row["ingredient_id"], "name": row["name"], "kind": row["kind"],
           "compartment": row["compartment"], "source": row["source"], "pdb_id": row["pdb_id"],
           "file": row.get("file", ""), "flags": [], "members": []}
    if not row.get("file"):
        out["flags"].append("no structure file")
        return out
    try:
        chains = load(STRUCT / row["file"])
    except Exception as e:
        out["flags"].append(f"could not read structure: {e}")
        return out
    accs = [m["accession"] for m in members]
    is_af = row["source"] == "afdb"

    if row["kind"] == "assembly":
        mapped = [c for c in chains if c["accs"] & set(accs)]
        if not mapped:
            keep = chains
        else:
            keep_rna = str(row.get("keep_nucleic_acid", "")).lower() in ("yes", "true", "1")
            keep_other = str(row.get("keep_other_chains", "")).lower() in ("yes", "true", "1")
            keep = [c for c in chains if (c["accs"] & set(accs))
                    or (keep_rna and ("RNA" in c["type"] or "DNA" in c["type"]))
                    or (keep_other and c["accs"])]
            dropped = len(chains) - len(keep)
            if dropped:
                out["flags"].append(f"dropped {dropped} chains not in the member list")
    elif is_af:
        keep = chains
    else:
        target = accs[0]
        keep = [c for c in chains if target in c["accs"]]
        if not keep and row["chain_id"]:
            ent = next((c["entity"] for c in chains
                        if c["name"] == row["chain_id"] or c["name"].startswith(row["chain_id"] + "-")), None)
            keep = [c for c in chains if c["entity"] == ent] if ent else []
        if not keep:
            keep = chains
            out["flags"].append("could not identify target chains; kept all chains")
    if not keep:
        out["flags"].append("no polymer chains")
        return out
    if row["kind"] == "protein" and len(keep) > opts["max_oligomer"]:
        out["flags"].append(f"structure has {len(keep)} copies (large particle?); kept one chain")
        keep = keep[:1]
    X, n_res = coords_of(keep, opts["min_plddt"] if is_af else None)
    if is_af:
        X_all, n_all = coords_of(keep)
        if n_res < 0.3 * n_all:               # mostly disordered: keep everything rather than a fragment
            X, n_res = X_all, n_all
            out["flags"].append(f"mostly low pLDDT (<{opts['min_plddt']:g}); kept all residues")
    if len(X) < 4:
        out["flags"].append("too few atoms after filtering")
        return out
    out.update(spheres_from(X, n_res, opts["res_per_sphere"], opts["max_spheres"]))
    if row["compartment"] in ("OM", "IM"):
        try:
            frame = membrane_frame(keep, members, row["compartment"], row["kind"], X.mean(0))
        except Exception as e:                  # never let orientation break the ingredient
            frame = None
            out["flags"].append(f"membrane orientation failed: {e}")
        if frame:
            out["membrane"] = frame
    out["chains_kept"] = [c["name"] for c in keep]
    out["n_residues"] = n_res

    # stoichiometry and copies per cell
    override = {}
    for part in [p for p in str(row.get("stoichiometry", "")).split(";") if ":" in p]:
        g, _, n = part.partition(":")
        override[g.strip().lower()] = int(n)
    have_refs = any(c["accs"] for c in keep)
    for m in members:
        if is_af or row["kind"] == "protein":
            stoich = len(keep)
        else:
            stoich = sum(1 for c in keep if m["accession"] in c["accs"])
        m["stoichiometry"] = stoich if (have_refs or is_af or row["kind"] == "protein") else 1
        if m["gene"].lower() in override:
            m["stoichiometry"] = override[m["gene"].lower()]
    if override:
        out["flags"].append("stoichiometry set in assemblies.csv")
    elif row["kind"] == "assembly" and not have_refs:
        out["flags"].append("no UniProt mapping in file (homolog?); stoichiometry assumed 1")
    with_data = [m for m in members if m["copies_per_cell"] > 0 and m["stoichiometry"] > 0]
    if row.get("mode") == "sum":
        # paralogs that make the same protein: pool all copies, divide by chains in the structure
        total = sum(m["copies_per_cell"] for m in members)
        chains_per = sum(m["stoichiometry"] for m in members) or len(keep)
        copies = total / max(1, chains_per)
        if not total:
            out["flags"].append("no copy number")
    elif not with_data:
        copies = 0.0
        out["flags"].append("no copy number")
    else:
        copies = statistics.median(m["copies_per_cell"] / m["stoichiometry"] for m in with_data)
    out["copies_per_cell"] = int(round(copies))
    mass = sum(m["mass_da"] * m["stoichiometry"] for m in members if m["stoichiometry"])
    na_res = sum(len(c["span"]) for c in keep if "RNA" in c["type"] or "DNA" in c["type"])
    if na_res:
        mass += na_res * 330                   # ~330 Da per nucleotide
        out["nucleic_acid_residues"] = na_res
    out["mass_da"] = int(mass) if mass else int(n_res * 110)
    absent = [m["gene"] for m in members if m["stoichiometry"] == 0]
    if absent and row["kind"] == "assembly" and row.get("mode") != "sum":
        out["flags"].append(f"members not in structure: {', '.join(absent)}")
    if row["compartment"] in ("OM", "IM") and "membrane" not in out:
        out["flags"].append("no TM segments mapped: orientation from the longest axis (approximate)")
    elif out.get("membrane") and not out["membrane"]["outward_known"]:
        out["flags"].append("membrane normal from TM segments, but which side faces out is unknown")
    if is_af:
        out["flags"].append("AlphaFold monomer: oligomeric state not modelled")
    out["members"] = [{k: v for k, v in m.items() if k not in ("seq", "tm", "topo")} for m in members]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--res-per-sphere", type=float, default=12)
    ap.add_argument("--max-spheres", type=int, default=80)
    ap.add_argument("--min-plddt", type=float, default=50)
    ap.add_argument("--max-oligomer", type=int, default=24,
                    help="single proteins with more copies than this in their structure keep one chain")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    opts = {"res_per_sphere": a.res_per_sphere, "max_spheres": a.max_spheres, "min_plddt": a.min_plddt,
            "max_oligomer": a.max_oligomer}

    prot = {p["accession"]: p for p in read_csv(DATA / "proteins.csv")}
    seq_file = DATA / "sequences.json"
    seqs = json.loads(seq_file.read_text()) if seq_file.exists() else {}
    if not seqs:
        print("Note: data/sequences.json not found; re-run steps 1 and 2 to orient membrane proteins.")
    rows = read_csv(DATA / "structures.csv")
    asm_opts = {a["id"]: a for a in read_csv(DATA / "assemblies.csv")}
    for r in rows:                             # chain-keeping options come from assemblies.csv
        a_ = asm_opts.get(r["ingredient_id"]) if r["kind"] == "assembly" else None
        if a_:
            r["keep_nucleic_acid"] = a_.get("keep_nucleic_acid", "")
            r["keep_other_chains"] = a_.get("keep_other_chains", "")
            r["stoichiometry"] = a_.get("stoichiometry", "")
    jobs = []
    for r in rows:
        members = [{"accession": acc, "gene": prot[acc]["gene"],
                    "copies_per_cell": fnum(prot[acc]["copies_per_cell"]),
                    "mass_da": fnum(prot[acc]["mass_da"])} for acc in r["members"].split(";") if acc in prot]
        if r["compartment"] in ("OM", "IM"):
            for m in members:
                p = prot[m["accession"]]
                m["seq"] = seqs.get(m["accession"], "")
                m["tm"] = parse_ranges(p.get("tm_segments", ""))
                m["topo"] = parse_ranges(p.get("topology", ""), with_label=True)
        jobs.append((r, members, opts))

    (OUT / "ingredients").mkdir(exist_ok=True)
    results = []
    with ProcessPoolExecutor(a.workers) as pool:
        for i, res in enumerate(pool.map(process, jobs, chunksize=4), 1):
            results.append(res)
            (OUT / "ingredients" / f"{res['id'].replace('/', '_')}.json").write_text(json.dumps(res))
            if i % 100 == 0:
                print(f"  {i}/{len(jobs)} done")

    usable = [r for r in results if r.get("spheres_A") and r.get("copies_per_cell", 0) > 0]
    bundle = {"organism": ORGANISM, "proteome": PROTEOME_ID, "units": "angstrom",
              "generated": datetime.date.today().isoformat(),
              "settings": opts, "ingredients": usable}
    (OUT / "ingredients.json").write_text(json.dumps(bundle))
    write_csv(OUT / "ingredients_summary.csv", [{
        "id": r["id"], "kind": r["kind"], "compartment": r["compartment"],
        "copies_per_cell": r.get("copies_per_cell", ""), "mass_kda": round(r.get("mass_da", 0) / 1000, 1),
        "n_spheres": len(r.get("spheres_A", [])), "encapsulating_radius_A": r.get("encapsulating_radius_A", ""),
        "source": r["source"], "pdb_id": r["pdb_id"], "file": r["file"], "flags": " | ".join(r["flags"]),
    } for r in results])

    print(f"{len(usable)} usable ingredients -> output/ingredients.json "
          f"({len(results) - len(usable)} skipped; see output/ingredients_summary.csv)")
    # crowding check against ~30% w/v (180.7 Da per nm^3) in a ~0.69 um^3 cytoplasm
    cyto = sum(r["copies_per_cell"] * r["mass_da"] for r in usable if r["compartment"] in ("CYTO", "NUC"))
    vol_nm3 = 0.69e9
    print(f"Cytoplasm + nucleoid macromolecular mass: {cyto / vol_nm3 / 180.7 * 30:.1f}% w/v "
          f"(target ~30%, assuming a 0.69 um^3 cytoplasm; rRNA included via the ribosome)")


if __name__ == "__main__":
    main()
