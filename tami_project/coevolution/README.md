# TamI 共進化解析 索引

| 段階 | 場所 | 内容 | 状態 |
| --- | --- | --- | --- |
| Round 1 (pilot) | [round1_pilot/](./round1_pilot/) | 最小キットによる MIp pilot（元 MSA 69 配列ベース） | 完了・参照用 |
| Round 1 (pilot MSA/系統樹) | `../metadata/msa/TamI/`, `../data/processed/msa/`, `../data/processed/tree/`, `../data/raw/homologs/TamI/` | Swiss-Prot/Bacteria 検索、MAFFT MSA、IQ-TREE。parameter は `../../config/p450_phylo_pilot_v0.1.json`、手順は `../../reports/TamI_TxtC_MSA_tree_workflow_v0.1.md` | 完了。他体系と共通 script を使うため元の場所に保持 |
| Round 2 | [round2/](./round2/) | homolog set 拡張、Neff、MIp/APC、DCA | **最終結果 2026-10-20** |

Round 2 の入口: [round2/README.md](./round2/README.md) と `round2/TamI_round2_supplemental_search_v1/README_FINAL.md`。Round 2 Neff レポートは `../reports/TamI_round2_MSA_Neff_report_v1.md`。

Round 1 → Round 2 の変化: 配列数 69 → 457 (Neff/L 0.474)、MIp のみ → MIp/APC + DCA。
