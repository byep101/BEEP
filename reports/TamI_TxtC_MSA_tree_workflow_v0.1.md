# TamI / TxtC MSA・系統樹構築 実行記録 v0.1

実行日: 2026-08-04

## 結論

TamI と TxtC を分離し、reference 確認、homolog 検索・filter、MAFFT full-length MSA、残基/motif QC、smart-gap tree MSA、IQ-TREE ModelFinder + SH-aLRT 1000 + UFBoot 1000 まで完了した。

本成果物は pilot である。TamI は推奨配列数を満たすが、TxtC は 44 配列で手順書の推奨 50 配列未満のため、最終共進化解析の前により広いデータベースで homolog set を拡張する。

## Homolog 検索

- 採用データベース: NCBI preformatted Swiss-Prot、更新日 2026-07-31
- 対象分類群: Bacteria (taxid 2)
- BLAST: local `blastp-fast`, max_target_seqs 500, E-value <= 1e-20
- filter: identity 25–95%、query coverage >= 70%、長さ比 0.75–1.25、X <= 1%、完全重複除外
- remote nr/refseq_protein の Bacteria/Actinomycetota 検索は NCBI CPU limit または長時間 queue/stall により採用せず、失敗ログを各 system の `logs/` に保存した。

## 結果サマリー

| System | raw hits | MSA sequences | full MSA | tree MSA | IQ-TREE model | MSA QC | tree QC |
| --- | ---: | ---: | ---: | ---: | --- | --- | --- |
| TamI | 68 | 69 | 654 columns | 521 columns | Q.pfam+F+R6 | qc_passed | qc_passed |
| TxtC | 43 | 44 | 518 columns | 518 columns | Q.pfam+F+R5 | qc_passed | qc_passed |

## Reference・残基 QC

- TamI reference: `ADC79647.1` / `TamI_REF_ADC79647.1`, 413 aa。Leu101 -> MSA 239、Leu244 -> 418、Leu295 -> 473。
- TxtC homolog search reference: `AAW49320.1` / `TxtC_REF_AAW49320.1`, 395 aa。
- TxtC structure reference `Q8VS75` は AAW49320.1 と同長で indel なし、identity 89.8734%（40 substitutions）。Q8VS75/PDB Ser280 は AAW49320.1 Ala280 に対応し、Thr385 は双方 Thr385。
- P450 motif anchor は TamI E305/C359、TxtC E272/C344 が reference 上で non-gap、MSA mapping はすべて pass。

## MSA QC

- MAFFT L-INS-i 相当 (`--localpair --maxiterate 1000`) を使用。
- TamI: mean column gap fraction 0.373266、gap >= 0.5 の列 250。smart-gap threshold 0.971、654 -> 521 columns。
- TxtC: mean column gap fraction 0.22135、gap >= 0.5 の列 117。smart-gap threshold 0.9773、518 -> 518 columns。閾値判定の結果、TxtC では削除列なし。
- ClipKIT 2.13.2 は Windows application-control policy により Biopython compiled extension をロードできなかったため、同版 `smart_gap_helper.py` の gap-distribution/slope rule を標準ライブラリのみで互換実装した。full MSA は保持した。

## 系統樹 QC

- IQ-TREE 2.4.0、amino-acid、ModelFinder Plus、SH-aLRT 1000、ultrafast bootstrap 1000。
- 両 treefile で reference の存在、leaf count、branch length を機械検査した。詳細は各 `metadata/msa/<System>/<System>_tree_qc_summary.json` を参照。
- IQ-TREE composition chi-square は候補配列の一部を warning とした。この検定は除外 filter としては使わず、pilot tree の解釈上の注意として `.iqtree` と `.log` に保持する。

## 再現性

- parameter: `config/p450_phylo_pilot_v0.1.json`
- processing utilities: `scripts/p450_phylo.py`, `scripts/concat_fasta.py`
- tool versions: `metadata/msa/tool_versions.txt`
- per-system search log: `metadata/msa/<System>/homolog_search_log.csv`
- per-system summary: `metadata/msa/<System>/msa_tree_summary.csv`
- IQ-TREE の完全 command line、seed、model、support 設定は各 `.iqtree` / `.log` に保存される。
