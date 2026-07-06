# TxtC literature note v0.1

## System

- System ID: TxtC
- BGC/product: thaxtomin A biosynthetic gene cluster / thaxtomin A
- BGC organism: Streptomyces turgidiscabies Car8
- MIBiG: BGC0000444.5
- BGC locus: AY707081.1
- Protein IDs: NCBI GenPept AAW49320.1; UniProt Q8VS75
- Structure source: Streptomyces acidiscabies 84-104
- Structure: PDB 6F0B and 6F0C

## Why TxtC is the second benchmark

TxtC is a bacterial P450 with BGC information, a defined thaxtomin pathway context, substrate-bound structures, and a published mechanism centered on substrate conformational switching. It is useful as the second standardized system after TamI because it tests whether future residue-level analyses can detect functionally relevant regions that are not limited to the heme-proximal conserved motifs.

## Source split and sequence caveat

The BGC/NCBI source and structure source are not the same strain:

- MIBiG and NCBI FASTA: Streptomyces turgidiscabies Car8, AAW49320.1.
- UniProt and PDB DBREF: Streptomyces acidiscabies 84-104, Q8VS75.

Both sequences are 395 aa, but they are not identical. The local comparison found 40 residue substitutions. Ser280 is Ser in Q8VS75/PDB numbering but Ala at position 280 in AAW49320.1. Thr385 is Thr in both sequences. Use `metadata/TxtC_AAW49320.1_vs_Q8VS75_sequence_mismatches_v0.1.csv` before residue-level mapping or benchmark interpretation.

## Known or candidate functional residues

Primary evidence: Alkhalaf et al., 2019, JACS, "Binding of Distinct Substrate Conformations Enables Hydroxylation of Remote Sites in Thaxtomin D by Cytochrome P450 TxtC" (PMID: 30516965; DOI: 10.1021/JACS.8B08864).

| Candidate residue | Numbering context | Structural evidence | Benchmark interpretation |
| --- | --- | --- | --- |
| Ser280 | Q8VS75/PDB numbering; Ala280 in AAW49320.1 | ATOM coordinates present in 6F0B and 6F0C. Listed in 6F0B SITE records near substrate/heme-associated sites. | Candidate mechanism residue for conformational switching; do not treat as a direct BGC-sequence Ser without mapping. |
| Thr385 | Q8VS75/PDB and AAW49320.1 both Thr385 | ATOM coordinates present in 6F0B and 6F0C. Listed in SITE records in both PDB files. | Candidate mechanism residue for remote-site hydroxylation / substrate conformation logic. |

## Known reaction logic

TxtC is associated with sequential aliphatic and aromatic hydroxylation in thaxtomin biosynthesis. The key mechanistic theme for this dataset is substrate conformational switching: distinct substrate conformations enable hydroxylation at remote sites.

## Data to compare against future predictions

- Do EVcouplings/MIp top-ranked residues recover Ser280, Thr385, or their local structural network in the Q8VS75/PDB coordinate system?
- Are predicted residues mapped to the substrate-binding pocket, access channel, or dynamic conformational-control regions?
- Does the signal distinguish remote conformational switching residues from only heme-proximal conserved P450 motifs?
- If using AAW49320.1 as the BGC-sequence anchor, are PDB/UniProt residues mapped through the local mismatch table first?

## Local raw data package

- Collection date: 2026-07-06
- MIBiG JSON: `data/raw/mibig/BGC0000444.5_annotations.json`
- NCBI FASTA: `data/raw/ncbi/AAW49320.1_TxtC_Streptomyces_turgidiscabies_Car8.fasta`
- UniProt FASTA: `data/raw/uniprot/Q8VS75_TxtC_UniProt.fasta`
- PDB/mmCIF structures: `data/raw/pdb/6F0B_TxtC.pdb`, `data/raw/pdb/6F0C_TxtC.pdb`, `data/raw/pdb/6F0B.cif`, `data/raw/pdb/6F0C.cif`
- Data collection log: `metadata/data_collection_log_v0.1_TxtC.csv`
- Sequence mismatch table: `metadata/TxtC_AAW49320.1_vs_Q8VS75_sequence_mismatches_v0.1.csv`

The MIBiG JSON parsed successfully and records BGC0000444, thaxtomin A, Streptomyces turgidiscabies Car8, AY707081.1, quality `questionable`, status `active`, and completeness `unknown`. The PDB headers distinguish 6F0B and 6F0C as TxtC structures for substrate conformational switching and record Q8VS75 as the UniProt sequence source.
