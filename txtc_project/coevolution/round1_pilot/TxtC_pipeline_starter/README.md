# BEEP / TxtC reproducible pipeline starter (2026-10-28)

**Status:** This round-one snapshot includes curated candidates and preparation outputs. The recorded focus alignment has only 8 sequences (Neff 1.0; `LOW_DEPTH_REVIEW`). The toolkit smoke test passed on synthetic data; MIp, DCA, and structure-contact analyses are not completed biological results.

Biological reference: UniProt Q8VS75 / Streptomyces acidiscabies 84-104. Experimental structures: PDB 6F0B (thaxtomin D, ligand C8H) and PDB 6F0C (monohydroxylated intermediate, ligand C8B). TxtC is not TxtE.

## Preflight

Run within WSL2/Ubuntu, from this round-one working directory (`$HOME/BEEP/txtc_project/coevolution/round1_pilot`). Requires MAFFT, an available `plmc` executable, and Python with `biopython`, `numpy`, `pandas`, `scipy`. Reuse the existing TamI environment if it contains these packages.

```
command -v mafft
command -v plmc
python3 --version
python3 -m pip install -r TxtC_pipeline_starter/requirements.txt
```

If `plmc` is not on PATH, use the full executable path or prepend its bin directory to PATH. `plmc --help` confirms options. Do not reinstall/rebuild already working EVcouplings/plmc without reason.

## Verify reference + structures, then curate candidate hits manually

```
cd "$HOME/BEEP/txtc_project/coevolution/round1_pilot"
test -s raw/reference/TxtC_Q8VS75.fasta
test -s raw/structure/6F0B.cif
test -s raw/structure/6F0C.cif
```

The round-one reference and structures are already present in this snapshot. If refreshing them, save new versions separately and record their provenance and SHA256 rather than overwriting these inputs. Search UniRef90, NCBI BLASTP nr/ClusteredNR separately with the reference protein. Export a **verified protein FASTA** containing manually curated candidates to `raw/search/txtc_curated_candidates.fasta` (no nucleotide sequences). Keep original search tables and accession mapping. Filter initially using E-value <= 1e-10, explicit `qcovs` >= 70, identity >= 25, length ratio 0.70..1.30, but **these are provisional screening cutoffs**: require biological/BGC context to distinguish TxtC from other P450s and TxtE. BLAST `ppos` must never be read as query coverage. Keep explicit column headers; don't identify fields by ordinal position.

## Start from curated FASTA

```
python3 TxtC_pipeline_starter/scripts/01_prepare_fasta.py \
  --reference raw/reference/TxtC_Q8VS75.fasta \
  --candidates raw/search/txtc_curated_candidates.fasta \
  --output raw/search/txtc_clean.fasta --manifest qc/candidate_manifest.csv

mafft --auto --thread 4 raw/search/txtc_clean.fasta > align/mafft_raw.fasta 2> logs/mafft.log

python3 TxtC_pipeline_starter/scripts/02_focus_weights.py \
  --msa align/mafft_raw.fasta --outdir align/focus \
  --min-row-coverage 0.70 --max-column-gap 0.50 --identity-threshold 0.80
```

The `--anchors 237:T` option may be used ONLY after the residue numbering and amino-acid at that position are independently checked for this reference and PDB. **Never copy TamI E287/R290/C359 indices to TxtC.** The focus script writes a `position_map.csv` that tracks the new model positions to Q8VS75 indices, avoiding positional substitution errors. It also writes `align/focus/weights.txt` in the exact FASTA row order. Verify `qc.json` and the alignment before proceeding.

## MIp / APC

```
python3 TxtC_pipeline_starter/scripts/03_mip.py \
  --msa align/focus/txtc_focus.fasta --weights align/focus/weights.txt \
  --position-map align/focus/position_map.csv \
  --output results/mip/mip_allpairs.csv
```

Computes weighted 21-state mutual information (including gap as a state), then applies standard off-diagonal APC. This is the documented toolkit definition; it is not certified equivalent to the user's previous TamI MIp implementation. Freeze this choice and, if desired, run a gap-handling sensitivity analysis.

## plmc DCA / CN

EVcouplings default-style initial coupling regularization (`lambda_J=0.01`, `lambda_J_times_Lq=True`, `q=21`) gives a literal CLI coupling penalty of `0.01*21*L`. This is an **initial pilot setting**, not a claim that the user's TamI run used the same penalty.

```
L=$(python3 -c "import json;print(json.load(open('align/focus/qc.json'))['L'])")
LE=$(python3 -c "print(0.01*21*int('$L'))")
set -o pipefail
plmc -c results/dca/plmc_couplings.txt \
  -o results/dca/plmc.params \
  -f TXT_C_REF_Q8VS75 -w align/focus/weights.txt -g \
  -m 100 -n 4 -lh 0.01 -le "$LE" \
  align/focus/txtc_focus.fasta 2>&1 | tee logs/plmc.log
```

Check `logs/plmc.log`: if the optimizer stopped with `ROUNDING_ERROR`, `MINIMUMSTEP`, etc., record this and label all scores exploratory until stability checks pass; writing a coupling file does not imply formal convergence. plmc outputs six whitespace-delimited fields per pair, with last field the APC-corrected Frobenius norm/CN score.

## Compare / benchmark

```
python3 TxtC_pipeline_starter/scripts/04_compare.py \
  --mip results/mip/mip_allpairs.csv \
  --plmc results/dca/plmc_couplings.txt \
  --position-map align/focus/position_map.csv \
  --site-qc align/focus/site_qc.csv \
  --outdir results/compare --min-separation 6
```

Strictly requires exactly `L*(L-1)/2` unique residue pairs from each method; the report's `Top20/Top50` uses variable sites and reference-sequence separation >=6. The `residue_best_pair_rank.csv` metric is *defined here* as best ranked eligible pair per residue and MUST NOT be conflated with the historic TamI residue-rank definition unless recomputed identically.

## Structure mapping / contacts

```
python3 TxtC_pipeline_starter/scripts/05_structure_contacts.py \
  --reference raw/reference/TxtC_Q8VS75.fasta \
  --position-map align/focus/position_map.csv \
  --pairs results/compare/dca_ranked.csv \
  --structure raw/structure/6F0B.cif --chain A \
  --out-prefix results/structure/6F0B_dca

python3 TxtC_pipeline_starter/scripts/05_structure_contacts.py \
  --reference raw/reference/TxtC_Q8VS75.fasta \
  --position-map align/focus/position_map.csv \
  --pairs results/compare/dca_ranked.csv \
  --structure raw/structure/6F0C.cif --chain A \
  --out-prefix results/structure/6F0C_dca
```

Output: reference-to-PDB mapping, pairwise heavy atom and C-beta distances, residue distances to HEM/substrate, top-K contact enrichment. Unknown/missing coordinates are **non-evaluable**, not zero-distance or non-contact. The script suggests the nearest cysteine-SG to heme Fe on structural evidence; verify visually before naming it an axial ligand.

## Smoke test (mock data only)

```
python3 TxtC_pipeline_starter/tests/smoke_test.py
```

## Citation and reference links

- UniProt Q8VS75 https://www.uniprot.org/uniprotkb/Q8VS75/entry
- RCSB PDB 6F0B https://www.rcsb.org/structure/6F0B
- RCSB PDB 6F0C https://www.rcsb.org/structure/6F0C
- Alkhalaf LM et al. JACS (2019) DOI https://doi.org/10.1021/jacs.8b08864
- Healy et al. (2002) DOI https://doi.org/10.1128/JB.184.7.2019-2029.2002
- plmc https://github.com/debbiemarkslab/plmc
- NCBI BLAST field definitions https://www.ncbi.nlm.nih.gov/books/NBK279684/

No repository or GitHub data was edited by this toolkit.
