# TamI validation checklist v0.1

Validation date: 2026-06-12

System: TamI / tirandamycin biosynthetic gene cluster / Streptomyces sp. 307-9

## Primary local files

- NCBI FASTA: `data/raw/ncbi/ADC79647.1_TamI_Streptomyces_307-9.fasta`
- UniProt FASTA: `data/raw/uniprot/D3Y1J3_TamI_UniProt.fasta`
- PDB structures: `data/raw/pdb/6XA2_TamI_tirandamycinC_bound.pdb`, `data/raw/pdb/6XA3_TamI_ligand_free.pdb`
- mmCIF structures: `data/raw/pdb/6XA2.cif`, `data/raw/pdb/6XA3.cif`
- Core metadata: `metadata/core_bgc_p450_dataset_v0.1_TamI.csv`
- Literature note: `literature_notes/TamI_literature_note.md`

## Checklist

- [x] NCBI and UniProt FASTA sequences were compared after removing FASTA headers and line wrapping.
  - Result: identical protein sequences.
  - Length: 413 aa in both files.
  - Sequence-only SHA256: `fbc563a7aeebba572bac52607531eb227351b690037bf8108f176a3f8701027a`.
  - Note: file-level hashes differ because headers and line wrapping differ.
- [x] Leu101, Leu244, and Leu295 are locatable in 6XA2 and 6XA3 structure coordinates.
- [x] TamI key residue mutation/evidence table added to `literature_notes/TamI_literature_note.md`.
- [x] README updated with v0.1 benchmark readiness statement.
- [x] Archived legacy reference FASTA copies remain traceable, but the official raw NCBI and UniProt downloads are the primary FASTA sources for v0.1.

## Structure residue localization

Residue numbering below uses PDB author residue numbering.

| Residue | 6XA2 coordinate evidence | 6XA3 coordinate evidence | Notes |
| --- | --- | --- | --- |
| Leu101 | `LEU A 101`, ATOM records present in chain A; also listed in SITE records near HEM/V0A sites. | `LEU A 101`, ATOM records present in chain A; also listed in the HEM SITE record. | Bicyclic ketal binding-site residue. |
| Leu244 | `LEU A 244`, ATOM records present in chain A; listed in the 6XA2 HEM SITE record. | `LEU A 244`, ATOM records present in chain A. | Bicyclic ketal binding-site residue. |
| Leu295 | `LEU A 295`, ATOM records present in chain A; also listed in SITE records near HEM/V0A sites. | `LEU A 295`, ATOM records present in chain A. | Bicyclic ketal binding-site residue. |

## Validation outcome

TamI v0.1 passes the sequence and structure checks needed for the first benchmark system.
