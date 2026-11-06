# E. coli virtual cell – data pipeline (v0.1.5)

Builds packing ingredients for a whole-cell model of *E. coli* K-12 MG1655
(UniProt proteome UP000000625): which molecules, how many per cell, where they
live, and a coarse 3D shape for each. The output feeds autoPACK/cellPACK and
the browser viewer.

## Setup

```bash
python -m venv venv && source venv/bin/activate     # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Steps

| Step | Command | Makes | Time (approx.) |
|---|---|---|---|
| 1 | `python step1_fetch_proteome.py` | `data/proteome.csv` – all ~4,400 proteins with gene, b-number, mass, compartment | 1 min |
| 2 | `python step2_add_abundance.py FILE ...` | `data/proteins.csv` – adds copies per cell | seconds |
| 3 | `python step3_pick_structures.py` | `data/structures.csv` – one structure per ingredient | 10–20 min first run (cached afterwards) |
| 4 | `python step4_download_structures.py` | `structures/*.cif(.gz)` | 20–60 min, roughly 1–3 GB for everything |
| 5 | `python step5_coarse_grain.py` | `output/ingredients.json` + summary CSV | 5–20 min |

Or run everything: `./run_all.sh li2014_tableS1.xlsx --value-col "MOPS complete" --skiprows 1`
(set `TOP=500` in front to limit step 3 to the 500 proteins carrying the most mass).

For a first test, use `python step3_pick_structures.py --top 300`. The few hundred
most abundant proteins make up most of the cell's protein mass, so this gives a
realistic-looking cell quickly.

### Step 2: abundance data

You supply the file, because neither source offers a stable direct link.

- **Li et al. 2014, Cell 157:624** (the dataset Goodsell used): download the
  supplementary table with absolute synthesis rates. Values are molecules per
  generation, which approximate copies per cell for stable proteins. Run the
  script once without `--value-col` to list the columns, then choose the MOPS
  complete column. Use `--skiprows` if the column names are not on the first row.
- **PaxDb** (pax-db.org): download an *E. coli* (511145) abundance file, then
  `--format paxdb --units ppm --total-proteins 3e6`. The total protein count per
  cell depends on growth rate (roughly 2–4 million), so treat it as a parameter.

Check `data/abundance_unmatched.txt` afterwards; a match rate below ~95% usually
means the wrong ID column.

### Step 3: how structures are chosen

1. **Curated assemblies** in `data/assemblies.csv` (ribosome, ATP synthase,
   complex I, RNA polymerase, GroEL/ES and others from Goodsell's Table 1) use a
   single PDB entry for several genes. `pdb_id = auto` takes the best entry for
   the first listed gene and keeps all its chains. `members` accepts
   `regex:` patterns (the ribosome uses `regex:^rp[lsm][A-Z]$`). Edit this
   file freely; it is the main place for expert curation.
2. **Single proteins** use the top-ranked entry from PDBe's best-structures
   service if it covers at least 60% of the sequence (`--min-coverage`),
   otherwise the AlphaFold DB model.

### Step 5: shapes and counts

- Heavy atoms are clustered into spheres (about one per 12 residues, at most
  80), the multi-sphere representation used for collision tests in packing.
- Single proteins keep only their own chains from the PDB file, so the
  stoichiometry of homo-oligomers (e.g. the OmpF trimer) is read from the
  biological assembly and copies per cell are divided accordingly.
- Assemblies get a count from the median of (copies of member ÷ chains of that
  member), so a 14:7 GroEL:GroES complex is counted correctly.
- AlphaFold residues with pLDDT < 50 are dropped unless that would remove most
  of the protein.

## Known gaps (flagged in `output/ingredients_summary.csv`)

- **Membrane orientation**: membrane proteins are not yet aligned to the
  bilayer. The OPM database (opm.phar.umich.edu) provides membrane-oriented
  coordinates; adding that is the natural next step.
- **AlphaFold models are monomers**: proteins without an experimental
  structure lose their oligomeric state. AlphaFold 3 predictions of the
  oligomer, using UniProt's "subunit" annotation in `data/proteome.csv`, would fix the most abundant cases.
- **Homolog structures**: assemblies built from another species' PDB entry have
  no E. coli UniProt mapping, so their stoichiometry is assumed to be 1.
- **Not yet included**: mRNA, tRNA, DNA, lipids, LPS, peptidoglycan and
  metabolites. The viewer handles DNA, membranes and the envelope separately.
- PDB biological assembly 1 is the authors' or software's best guess and is
  occasionally not the physiological state. Spot-check large complexes.

## Output format (`output/ingredients.json`)

Lengths in Å, coordinates centred on each ingredient's centroid:

```json
{"id": "groELS", "kind": "assembly", "compartment": "CYTO",
 "copies_per_cell": 3929, "mass_da": 858000, "source": "pdb", "pdb_id": "1aon",
 "members": [{"gene": "groL", "stoichiometry": 14, "copies_per_cell": 30000}, ...],
 "spheres_A": [[x, y, z, r], ...], "encapsulating_radius_A": 92.1,
 "radius_of_gyration_A": 58.3, "principal_vector": [0, 0, 1], "flags": []}
```

Compartments: OM outer membrane, PERI periplasm, IM inner membrane, CYTO
cytoplasm, NUC nucleoid, EXT outside the cell (fimbriae, flagella, secreted).
