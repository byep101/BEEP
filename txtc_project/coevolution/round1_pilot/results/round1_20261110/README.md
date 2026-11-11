# TxtC round-one coevolution analysis

**Run date:** 2026-11-10 (Asia/Tokyo)  
**Status:** Scores generated; optimization review required. This is an exploratory broad-P450 analysis, not a validated TxtC-specific or functional result.

The authoritative, operation-manual-compatible outputs are under [`10_official_method/`](./10_official_method/). Earlier numbered directories in this run are preflight, smoke, superseded implementation checks, or an interrupted run; do not use their scores as the final result.

## Results at a glance

| Metric | Result |
|---|---:|
| Filtered homolog candidates before motif review | 4,007 |
| Unique exclusions (P450 signature and reference-anchor review) | 148 |
| Final focus alignment | 3,860 sequences × 361 positions |
| Sequence-weighted Neff | 2,358.644 |
| Neff/L | 6.534 |
| MIp and DCA all-pair scores | 64,980 each |
| Eligible nonlocal variable pairs | 61,899 |
| MIp–DCA Spearman correlation | 0.3876 |
| Top-20 / Top-50 / Top-100 overlap | 3 / 20, 12 / 50, 23 / 100 |
| plmc optimizer status | `MAXIMUMITERATION`; **not certified converged** |

Neff/L exceeds the project's 0.3 execution threshold. The large value does not demonstrate TxtC-family purity: this is explicitly a broad bacterial P450 set.

## Input curation and residue numbering

The expansion combines the existing UniProtKB bacterial cytochrome-P450 candidate collection with the available UniRef90 Q8VS75 cluster snapshot. The filtered collection contains 4,007 exact-unique homolog sequences. Numeric intake thresholds were 25–95% identity, query coverage ≥70%, E-value ≤1e-20, bits/reference-residue ≥0.5, length ratio 0.75–1.25, and ambiguous residues ≤1%. The Q8VS75 sequence was added as the structural focus reference, while AAW49320.1 remains in the alignment as the BGC/NCBI reference.

The final model excludes 138 sequences failing the recorded PROSITE PS00086 P450 signature screen and 13 sequences with at least one E272/R275/C344 reference-anchor mismatch; three records were in both groups, for 148 unique exclusions. See [`exclusion_audit.tsv`](./10_official_method/input/exclusion_audit.tsv) for per-sequence reasons.

Only one of the expanded homolog records has an explicit `GN=txtC` annotation. The prior remote NCBI BLAST request failed at the server's CPU-use limit, and no NCBI BLAST hits were included in this run. A phylogenetic tree and BGC-neighborhood review of the expanded set were not completed. Therefore, the depth threshold is met for the broad-P450 pilot, not established for a curated TxtC ortholog set.

AAW49320.1 and Q8VS75 are both 395 aa and align without an insertion/deletion, but differ at 40 positions. Most notably, AAW49320.1 position 280 is Ala, whereas Q8VS75/PDB position 280 is Ser. The final focus MSA and structure mapping use Q8VS75 numbering; the mapping retains the AAW49320.1 reference row. Q8VS75 Thr237 maps to Thr237 in both PDB chain-A structures. The nearest mapped cysteine SG to heme Fe is Cys344 (2.518 Å in 6F0B; 2.352 Å in 6F0C). This distance is a structural candidate, not independent proof of catalytic assignment.

## Scoring and interpretation

The final run uses the repository's documented `02_focus_weights.py`, 21-state gap-inclusive `03_mip.py`, and `04_compare.py`, with the same external weights for MIp and plmc. The plmc initial setting follows the TxtC operation manual: lambda_h=0.01, lambda_e=0.01×21×361=75.81, four threads, and a 100-iteration cap. All 64,980 pairs passed file/coordinate checks. The optimizer reached `MAXIMUMITERATION`; do not call the model converged.

The highest-scoring eligible pairs were:

| Rank | MIp pair | MIp | DCA-CN pair | DCA-CN |
|---:|---|---:|---|---:|
| 1 | 126–365 | 0.48009 | 126–365 | 0.89022 |
| 2 | 16–23 | 0.37548 | 115–360 | 0.74848 |
| 3 | 101–207 | 0.37452 | 258–328 | 0.73192 |
| 4 | 16–24 | 0.36982 | 274–336 | 0.68533 |
| 5 | 15–24 | 0.35941 | 121–142 | 0.61271 |

Three pairs are shared between the two Top-20 lists: 126–365, 101–207, and 115–360. The MIp/DCA Top-20 contact check found 5/20 and 17/20 pairs with minimum heavy-atom distance ≤5 Å, respectively, in both structures; every Top-20 pair was evaluable. The eligible-pair background contact counts were 690/59,795 in 6F0B and 683/59,795 in 6F0C. These are descriptive structural checks, not a significance test or functional validation.

Thr237 is variable and maps correctly, but its best eligible partner ranks are 2,192 by MIp and 320 by DCA. Q8VS75 Ser280 (AAW49320.1 Ala280) ranks 12 by MIp and 221 by DCA. Thr385 ranks 1,048 and 247, respectively. Cys344 is invariant in this alignment and consequently has no eligible covariation rank; a low or absent rank must not be interpreted as lack of functional importance.

The manual's lambda_e sensitivity settings were also run with the focus MSA and weights frozen:

| lambda_e | Optimizer status | Score Spearman vs 1× | Top-20 / 50 / 100 DCA-pair overlap vs 1× |
|---:|---|---:|---:|
| 0.5× (37.905) | `ROUNDING_ERROR` | 0.9789 | 19 / 20, 49 / 50, 95 / 100 |
| 1× (75.81) | `MAXIMUMITERATION` | — | reference |
| 2× (151.62) | `MAXIMUMITERATION` | 0.9626 | 18 / 20, 48 / 50, 96 / 100 |

Top-pair scores are relatively stable across these two regularization perturbations, but the optimizer warnings remain and the underlying homolog-family scope is unresolved. See each `dca_sensitivity/` subdirectory for complete logs and couplings.

## Reproducibility and files

- [`10_official_method/run_manifest.json`](./10_official_method/run_manifest.json): thresholds, environment, hashes, run settings, and QC summary.
- [`10_official_method/commands.json`](./10_official_method/commands.json) and [`10_official_method/environment.txt`](./10_official_method/environment.txt): exact command lines and software versions.
- [`10_official_method/frozen_sha256.tsv`](./10_official_method/frozen_sha256.tsv): SHA-256 manifest for primary and sensitivity outputs.
- [`10_official_method/results/validation.json`](./10_official_method/results/validation.json): independent full-pair, score-finiteness, weight, and coordinate validation.
- [`10_official_method/align/focus/`](./10_official_method/align/focus/): focus MSA, external sequence weights, Q8VS75 position map, site QC, and Neff report.
- [`10_official_method/results/mip/`](./10_official_method/results/mip/): complete MIp/APC pair table.
- [`10_official_method/results/dca/`](./10_official_method/results/dca/): primary plmc couplings, parameter file, and full log.
- [`10_official_method/results/compare/`](./10_official_method/results/compare/): complete pair validation, eligible rankings, rank correlation, and top-set intersections.
- [`10_official_method/results/structure/`](./10_official_method/results/structure/): 6F0B/6F0C residue mappings, pair contacts, ligand distances, and summaries.

The support archive's six synthetic tests and the repository TxtC pipeline smoke test both passed. These tests validate software behavior only, not biological interpretation.

## Remaining work before a TxtC-specific conclusion

Obtain a successful NCBI Protein BLAST expansion, manually verify candidate BGC context and TxtC/TxtE paralog separation, and inspect a phylogeny before treating this broad-P450 model as a TxtC homolog set. Then repeat the frozen analysis on the curated set. Review the `MAXIMUMITERATION`/`ROUNDING_ERROR` logs and convergence stability before interpreting pair ranks; experimental validation remains separate.
