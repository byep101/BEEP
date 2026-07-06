# TxtC validation checklist v0.1

Validation date: 2026-07-06

System: TxtC / thaxtomin biosynthetic gene cluster / Streptomyces turgidiscabies Car8 BGC with Streptomyces acidiscabies 84-104 structure reference

## Primary local files

- MIBiG JSON: `data/raw/mibig/BGC0000444.5_annotations.json`
- NCBI FASTA: `data/raw/ncbi/AAW49320.1_TxtC_Streptomyces_turgidiscabies_Car8.fasta`
- UniProt FASTA: `data/raw/uniprot/Q8VS75_TxtC_UniProt.fasta`
- PDB structures: `data/raw/pdb/6F0B_TxtC.pdb`, `data/raw/pdb/6F0C_TxtC.pdb`
- mmCIF structures: `data/raw/pdb/6F0B.cif`, `data/raw/pdb/6F0C.cif`
- Core metadata: `metadata/core_bgc_p450_dataset_v0.1_TxtC.csv`
- Data collection log: `metadata/data_collection_log_v0.1_TxtC.csv`
- Sequence mismatch table: `metadata/TxtC_AAW49320.1_vs_Q8VS75_sequence_mismatches_v0.1.csv`
- Literature note: `literature_notes/TxtC_literature_note.md`

## Checklist

- [x] MIBiG JSON downloaded and parsed.
  - Result: BGC0000444.5, thaxtomin A, Streptomyces turgidiscabies Car8, AY707081.1.
  - MIBiG quality: `questionable`; status: `active`; completeness: `unknown`.
- [x] NCBI and UniProt FASTA files downloaded.
  - NCBI AAW49320.1 length: 395 aa.
  - UniProt Q8VS75 length: 395 aa.
  - NCBI sequence-only SHA256: `b7a3e9b9629cf6e8ff7882bfc4b151ec8e54a795e96fe673eafc8699fab4b9db`.
  - UniProt sequence-only SHA256: `ce85d2d1dcaf0e85a048be8d41e482fd547122e644f8fe218e7db367c29c4aa0`.
- [x] NCBI and UniProt sequence relationship documented.
  - Result: same length but not identical.
  - Difference count: 40 positions.
  - Important caveat: position 280 is Ala in AAW49320.1 but Ser in Q8VS75/PDB. Position 385 is Thr in both sequences.
  - Mapping support: `metadata/TxtC_AAW49320.1_vs_Q8VS75_sequence_mismatches_v0.1.csv`.
- [x] PDB structure files downloaded and checked.
  - 6F0B: X-ray diffraction, 2.80 A resolution, DBREF Q8VS75 residues 1-395.
  - 6F0C: X-ray diffraction, 1.70 A resolution, DBREF Q8VS75 residues 1-395.
- [x] Ser280 and Thr385 candidate residues are locatable in structure coordinates.
- [x] Core metadata CSV, data collection log, literature note, README, and operation manual are present.
- [x] Candidate residue status is conservative.
  - Ser280 and Thr385 are recorded as candidate mechanism residues, not as direct point-mutant validated benchmark residues.

## Candidate residue localization

Residue numbering below uses PDB/UniProt Q8VS75 numbering.

| Residue | 6F0B coordinate evidence | 6F0C coordinate evidence | Notes |
| --- | --- | --- | --- |
| Ser280 | `SER A 280`, ATOM records present in chain A; listed in 6F0B SITE records near substrate/heme-associated sites. | `SER A 280`, ATOM records present in chain A. | Q8VS75/PDB residue; the same position is Ala280 in AAW49320.1. |
| Thr385 | `THR A 385`, ATOM records present in chain A; listed in 6F0B SITE records. | `THR A 385`, ATOM records present in chain A; listed in 6F0C SITE records. | Conserved as Thr385 in both checked FASTA files. |

## Validation outcome

TxtC v0.1 passes the raw data collection and structure-localization checks needed for the second standardized benchmark package. It is suitable for structure-aware testing and schema validation.

Residue-level benchmark use requires care: the BGC/NCBI sequence and UniProt/PDB sequence are equal length but non-identical, so PDB candidate residues should be mapped explicitly before being interpreted against the BGC/NCBI sequence.
