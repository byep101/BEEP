# TxtC coevolution round-one pilot

This directory is an isolated snapshot of the TxtC round-one inputs and preparation outputs. The original working directory, `work/txtc_round1`, is retained unchanged.

## Status

- Candidate FASTA curation and MAFFT/focus-alignment preparation outputs are present.
- Recorded focus QC: 8 input sequences, Neff 1.0, and `LOW_DEPTH_REVIEW`.
- This original eight-sequence snapshot remains a low-depth pipeline pilot; its scores must not be interpreted biologically.
- A separate expanded run was executed on 2026-11-10. Its authoritative MIp/DCA, QC, and structure outputs are documented in [`results/round1_20261110/`](./results/round1_20261110/).
- The expanded run is exploratory broad-P450 analysis, not a curated TxtC ortholog set: the NCBI remote BLAST attempt failed at the server CPU limit, and BGC-context/phylogenetic review remains outstanding.
- The expanded run passed the project Neff/L threshold, but plmc reached `MAXIMUMITERATION`; its optimizer is not certified converged.
- `raw/structure/6F0B.cif` and `raw/structure/6F0C.cif` are retained as the exact pilot inputs. They differ byte-for-byte from the canonical copies under `../../data/raw/pdb/`; neither set has been overwritten.

## Layout

- `raw/`: reference, search results/candidates, and structure inputs.
- `align/`: raw and focus MSAs, position map, site QC, and sequence weights.
- `qc/`: candidate manifest and SHA256 records.
- `logs/`: preparation, environment, and smoke-test records.
- `TxtC_pipeline_starter/`: documented scripts, dependencies, and synthetic-data smoke test.
- `reports/`: reserved for additional reports.
- `results/round1_20261110/`: expanded first-round analysis; use its README to distinguish authoritative outputs from superseded preflight runs.

For prerequisites and commands, see `TxtC_pipeline_starter/README.md`. Run those commands from this directory so the documented relative input/output paths resolve correctly.
