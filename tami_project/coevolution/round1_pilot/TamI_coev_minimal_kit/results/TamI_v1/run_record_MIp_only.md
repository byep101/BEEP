# TamI coevolution pilot run record: TamI_v1

Run date: 2026-09-21  
System ID: TamI  
Overall status: **partial - prepare and MIp/APC complete; DCA not run**

## Scope and provenance

This run reused the completed TamI full MSA and reference FASTA without
modifying the source files. The starter archive named in
`TamI_MIp_DCA_pilot_Manual_v1.0.docx` was not present locally and was not in
the current official BEEP repository, so the standard-library prepare and
MIp/APC stages were reconstructed from the manual's documented rules. Five
focused unit tests for weighting, overlap handling, MI, and pairwise missing
data passed before and after the real-data run.

No substitute DCA algorithm was used.

## Inputs

| Input | SHA256 |
|---|---|
| `inputs/TamI_full_msa.faa` | `af2cb39fe1415cb426c21f6147040b2eb9893a09a288aa05e85e870b8355cced` |
| `inputs/TamI_reference.faa` | `2b3fb203b546637ae6892666697ddd305bfde639254c84cf2d42f673506eb69c` |
| Ungapped TamI reference sequence | `fbc563a7aeebba572bac52607531eb227351b690037bf8108f176a3f8701027a` |

The copied input hashes exactly match the source files under
`tami_project/data/processed/msa/TamI/`.

## Fixed pilot settings

- Reference ID: `TamI_REF_ADC79647.1`
- Sequence reference coverage: at least 70%
- Column occupancy: at least 70%
- Pairwise identity threshold: 0.8
- Minimum comparable positions for neighbors: 70% of model length
- Weight: reciprocal neighbor count, including self
- MIp: weighted empirical MI, natural logarithm, no pseudocount, pairwise
  exclusion of gaps/non-standard residues, APC over all model-column pairs
- Display eligibility: both sites variable and reference-position separation
  at least 6 residues

## Prepare result

| Metric | Result |
|---|---:|
| Input sequences | 69 |
| Input MSA columns | 654 |
| Reference length | 413 aa |
| Sequences passing coverage | 69 |
| Exact focus duplicates removed | 0 |
| Model sequences | 69 |
| Model columns | 384 |
| Neff | 57.333333 |
| Neff / L_model | 0.149306 |

This is a low-depth method pilot. The value of Neff/L_model is not a
reliability guarantee.

Leu101, Leu244, and Leu295 all match the reference and remain in the model:

| Known site | Full-MSA column | Model column | Status |
|---|---:|---:|---|
| Leu101 | 239 | 76 | retained, variable |
| Leu244 | 418 | 219 | retained, variable |
| Leu295 | 473 | 270 | retained, variable |

## MIp/APC result

| Technical check | Result |
|---|---|
| Expected/scored pairs | 73,536 / 73,536 |
| Eligible nonlocal variable pairs | 69,074 |
| Variable/constant model sites | 377 / 7 |
| Finite MI, MIp, and observable-weight values | pass |
| Full pair-table SHA256 | `4e7a4e6245d4ad65e78f014c682106e50c0c843e31062540c5379d218b487047` |

The highest eligible MIp pair is V148-L410 (MIp 0.414156). The known-site
results below are descriptive only. Residue rank uses the maximum eligible
nonlocal MIp for each site, with equal scores assigned the same minimum rank.

| Known site | Best partner | Maximum MIp | Pair rank | Residue rank | Residue top 20 |
|---|---|---:|---:|---:|---|
| Leu101 | Asn255 | 0.225836 | 126 | 156 | no |
| Leu244 | Asn255 | 0.283596 | 19 | 35 | no |
| Leu295 | Thr252 | 0.276925 | 23 | 41 | no |

These ranks do not prove or disprove functional importance. In particular,
the three known sites are not sufficient to estimate predictive accuracy.

## DCA status

Status: `not_run_environment_missing`

The workstation has no WSL distribution, Docker/Podman, gcc, make, plmc, or
EVcouplings installation. The manual requires the pinned EVcouplings `run_plmc`
interface and plmc with external-weight support in an approved Linux/WSL
environment. Therefore DCA was not approximated, and the cross-method compare
stage was not run.

To finish this run without changing its input definition, copy the complete
`TamI_coev_minimal_kit` directory to an approved Linux/WSL environment, install
the manual's pinned EVcouplings/plmc versions, run DCA against
`results/TamI_v1/focus_cov70.faa` and the existing `weights.txt`, verify the
optimizer log and pair completeness, and only then run the comparison stage.

## Completion statement

The input/numbering checks, one Neff calculation, and the complete MIp/APC
baseline are finished and reproducible. The full first-round completion
statement from the manual must not yet be used because DCA technical checks are
pending.
