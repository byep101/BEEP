# core_bgc_p450_dataset_v0.1 master validation checklist

Date: 2026-07-09

## Scope

This checklist validates the integration-layer dataset after adding TamI, TxtC, RosC, and MpaDE_MpaD to a shared master CSV.

## Files

- `metadata/core_bgc_p450_dataset_v0.1_master.csv`
- `metadata/per_system/core_bgc_p450_dataset_v0.1_TamI.csv`
- `metadata/per_system/core_bgc_p450_dataset_v0.1_TxtC.csv`
- `metadata/per_system/core_bgc_p450_dataset_v0.1_RosC.csv`
- `metadata/per_system/core_bgc_p450_dataset_v0.1_MpaDE_MpaD.csv`
- `metadata/data_collection_log_master_v0.1.csv`

## Schema checks

- [x] Master CSV uses the recommended v0.1 integration schema with 41 columns.
- [x] Master CSV contains exactly four systems: TamI, TxtC, RosC, and MpaDE_MpaD.
- [x] Each per-system CSV uses the same header as the master CSV.
- [x] Required identity fields are filled for every row: System_ID, Product, BGC_ID, Organism, P450_name, Ref_sequence_source, Data_status.
- [x] Dataset_role values are non-duplicated across the four systems.
- [x] Reaction_category, Evidence_level, Mapping_status, and Data_status use the controlled vocabulary from the operation manual.

## System-specific checks

- [x] TamI is retained as the template_system row and keeps Leu101, Leu244, and Leu295 as known-residue benchmark targets.
- [x] TxtC reference sequence convention is explicit: Q8VS75/UniProt/PDB is the structure-mapping reference and AAW49320.1 remains the BGC/NCBI anchor.
- [x] TxtC residue-numbering caveat is explicit: Ser280 is Ser in Q8VS75/PDB but Ala280 in AAW49320.1; Thr385 is conserved as Thr in both sequences.
- [x] RosC includes loop/mutant benchmark evidence: Pro107, Leu176 in FG loop, and the P107S/L176Q mutant are recorded with 2020/2025 literature support.
- [x] RosC includes structure accessions 9VGM, 9VGN, 9VGQ, 9VGP, and 9VGO, with local chain/residue validation left in Issues_to_check.
- [x] MpaDE_MpaD is recorded as one fusion/domain comparison system, not as unrelated MpaD and MpaE rows.
- [x] MpaDE_MpaD reference sequence convention is explicit: BK008023.1 / DAA35181.1 is the fusion-level reference; HQ731031.1 ADY00131.1 and ADY00132.1 are retained for domain-specific cross-reference.

## Open issues intentionally retained

- RosC: resolve the source split between MIBiG taxonomy Micromonospora carbonacea and RosC/RosD literature source Micromonospora rosaria before treating BGC and RosC sequence as identical biological evidence.
- RosC: verify chain IDs and residue coordinates in local PDB/mmCIF files before scoring or structural coevolution mapping.
- MpaDE_MpaD: define exact domain boundaries and interface-residue numbering before MSA or fusion-interface scoring.
- MpaDE_MpaD: decide whether future analysis uses a full-fusion MSA, P450-domain-only MSA, or both.

## Freeze decision

The master CSV satisfies the v0.1 completion rule from the operation manual: all four systems are registered under a shared schema, reference sequence conventions are explicit, and unresolved points are documented in Issues_to_check rather than left as blank cells.
