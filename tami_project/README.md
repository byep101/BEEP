# tami_project

This folder contains the initial standardized dataset skeleton for TamI-tirandamycin BGC.

## Contents
- metadata/core_bgc_p450_dataset_v0.1_TamI.csv: first TamI row for the core dataset.
- metadata/data_collection_log_v0.1_TamI.csv: source URLs, collection date, file sizes, SHA256 hashes, and validation notes for the local raw data package.
- literature_notes/TamI_literature_note.md: summary of known TamI function, benchmark residues, and local source files.
- validation_checklist_TamI_v0.1.md: final v0.1 validation checklist for FASTA consistency, structure residue localization, and benchmark readiness.
- data/raw/mibig/BGC0001052.5_annotations.json: official MIBiG annotation JSON.
- data/raw/ncbi/ADC79647.1_TamI_Streptomyces_307-9.fasta: official NCBI Protein FASTA.
- data/raw/uniprot/D3Y1J3_TamI_UniProt.fasta: official UniProt FASTA.
- data/raw/pdb/: official RCSB PDB and mmCIF files for 6XA2 and 6XA3.
- data/processed/legacy_reference_fastas/: archived pre-download reference FASTA copies that were previously placed under raw directories.
- reports/: group meeting report.
- manuals/: operation manual.
- coevolution/: round 2 MSA expansion, Neff, MIp/DCA inputs, results and scripts (round1_pilot/ and round2/; see coevolution/README.md). Final results: coevolution/round2/TamI_round2_supplemental_search_v1/README_FINAL.md.

## Local status
TamI raw data collection for v0.1 is complete as of 2026-06-09. The current package includes MIBiG JSON, NCBI and UniProt FASTA files, RCSB PDB/mmCIF structure files, metadata CSV, literature note, data collection log, operation manual, and the 6.10 progress report.

TamI v0.1 dataset is ready for use as the first benchmark system.
