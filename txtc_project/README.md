# txtc_project

This folder contains the initial standardized dataset skeleton for TxtC-thaxtomin BGC.

## Contents

- metadata/core_bgc_p450_dataset_v0.1_TxtC.csv: first TxtC row for the core dataset.
- metadata/data_collection_log_v0.1_TxtC.csv: source URLs, collection date, file sizes, SHA256 hashes, and validation notes for the local raw data package.
- metadata/TxtC_AAW49320.1_vs_Q8VS75_sequence_mismatches_v0.1.csv: residue-level mismatch table for the BGC/NCBI sequence versus the UniProt/PDB structure sequence.
- metadata/validation_checklist_TxtC_v0.1.md: final v0.1 validation checklist for source split, FASTA consistency, structure residue localization, and benchmark readiness.
- literature_notes/TxtC_literature_note.md: summary of TxtC function, candidate mechanism residues, sequence caveats, and local source files.
- data/raw/mibig/BGC0000444.5_annotations.json: official MIBiG annotation JSON.
- data/raw/ncbi/AAW49320.1_TxtC_Streptomyces_turgidiscabies_Car8.fasta: official NCBI Protein FASTA.
- data/raw/uniprot/Q8VS75_TxtC_UniProt.fasta: official UniProt FASTA.
- data/raw/pdb/: official RCSB PDB and mmCIF files for 6F0B and 6F0C.
- manuals/: operation manual.
- reports/: reserved for future reports.
- coevolution/round1_pilot/: isolated TxtC round-one input snapshot, alignment/QC, logs, and starter toolkit. See its README before using the results.
- coevolution/round1_pilot/results/round1_20261110/: expanded first-round MIp/DCA and structure-contact analysis with input/QC manifest.
- manuals/BEEP_TxtC_共进化分析操作手册_20261028.docx: round-one coevolution operation manual.

## Local status

TxtC raw data collection for v0.1 is complete as of 2026-07-06. The current package includes MIBiG JSON, NCBI and UniProt FASTA files, RCSB PDB/mmCIF structure files, metadata CSV, data collection log, literature note, validation checklist, sequence mismatch table, and the operation manual.

Important caveat: the BGC/NCBI sequence source is Streptomyces turgidiscabies Car8 (AAW49320.1), while the UniProt/PDB structure source is Streptomyces acidiscabies 84-104 (Q8VS75). Both sequences are 395 aa, but they are not identical. Ser280 is Ser in Q8VS75/PDB numbering and Ala at the same position in AAW49320.1. Use the sequence mismatch table before residue-level mapping.

TxtC v0.1 is ready as the second standardized benchmark system for data-organization and structure-aware testing, with the sequence-source caveat documented.

The original round-one pilot snapshot remains low depth (8 sequences, Neff 1.0). An expanded analysis is now recorded under `coevolution/round1_pilot/results/round1_20261110/`; it reaches Neff/L=6.534 and includes MIp/APC, plmc DCA, rank comparisons, and 6F0B/6F0C contact checks. Treat it as exploratory broad-P450 output only: the NCBI remote BLAST request failed at the server CPU limit, TxtC/TxtE and BGC-context review is incomplete, and plmc did not certify convergence (`MAXIMUMITERATION`). See that run's README for detailed limits and reproducibility hashes. All round-one inputs remain separate from the canonical `data/` package; existing canonical structures were not overwritten.
