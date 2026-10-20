# TamI 共進化解析 (round 2) ディレクトリ索引

TamI の MSA/MIp/DCA 解析 round 2 の全成果物。各スクリプトは兄弟ディレクトリを相対参照するため、このフォルダ内の相対配置を保つこと。ディレクトリ名は変更していない。

| ディレクトリ | 役割 | 状態 |
| --- | --- | --- |
| `../round1_pilot/TamI_coev_minimal_kit/` | v1 最小キット（pilot MIp、元 MSA）。round 1 に属する | 参照用 |
| `TamI_round2_manual_input/` | 手動取得した RefSeq 入力と baseline curation | 入力 |
| `TamI_round2_results/` | round 2 first pass の結果 | 完了（履歴） |
| `TamI_round2_expanded_input_v1/` | 拡張 250 件の RefSeq / ClusteredNR 入力 | 入力 |
| `TamI_round2_expanded_results_v1/` | 拡張セットの intake・MSA・QC | 完了（履歴） |
| `TamI_round2_supplemental_search_v1/` | UniRef90 補充検索、ClusteredNR qcovs 補正、最終 MSA、MIp/DCA | **最終結果 (2026-10-20)** |
| `TamI_round2_execution/` | round 2 全工程の実行スクリプト | 実行コード |

## 最終結果の入口

- 要約: `TamI_round2_supplemental_search_v1/README_FINAL.md`
- 最終 MSA / prepare: `TamI_round2_supplemental_search_v1/06_final_msa/`
- MIp / DCA: `TamI_round2_supplemental_search_v1/07_coevolution/`
- 整合性: `FINAL_SHA256SUMS.txt`, `SUPPLEMENTAL_SHA256SUMS.txt`

## 最終数値 (2026-10-20)

- 最終 full-length MSA 入力 462 配列 → focus 重複 5 を除去 → 457 配列 × 401 位置
- Neff = 190.20 (θ=0.8)、Neff/L = 0.474（目標 0.3 を達成）
- MIp/APC と plmc DCA はともに 80,200 pair を出力。DCA の最適化状態は `MINIMUMSTEP`（収束とは扱わない）
- Leu101/Leu244/Leu295 は residue top20 に入らない（MIp rank 253/238/105、DCA rank 363/157/73）。記述的 benchmark であり、機能重要性の証明・否定には使わない
- `TamI.params`（各 256 MB）は WSL に保管し、本フォルダには含めない

## 注意

`TamI_round2_execution/` の一部スクリプトは WSL の `~/projects/tami_dca/` を参照する（実行履歴として保持）。
