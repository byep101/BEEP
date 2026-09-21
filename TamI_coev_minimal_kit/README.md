# TamI coevolution pilot (reconstructed prepare/MIp stage)

This directory is an independent pilot workspace created from
`TamI_MIp_DCA_pilot_Manual_v1.0.docx`. The starter-kit archive named by the
manual was not available, so the standard-library `prepare` and `mi` stages
were reconstructed from the documented rules and covered by focused tests.

The original BEEP MSA and reference FASTA are copied into `inputs/`; they are
not modified. Results are written under `results/TamI_v1/`.

DCA is intentionally not replaced by another method. It must be run with the
manual's pinned EVcouplings/plmc Linux environment before the cross-method
comparison stage can be completed.
