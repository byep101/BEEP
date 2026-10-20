# TamI 第二轮扩充（RefSeq 250 + clusteredNR 250）

## 最终结果

已完成序列验证、固定阈值过滤、跨库合并、完整序列 MSA、模型位点过滤、θ=0.8 权重和 Neff 计算。最终分析使用 `combined_curated_expanded_v2.faa`；初次 MSA 仅作为 motif 审查记录保留。

| 指标 | 最终值 |
|---|---:|
| 原始来源记录 | 500（RefSeq 250 + clusteredNR 250） |
| 跨库唯一 accession | 439 |
| 固定数值 QC 合格的新序列 | 169 |
| MSA 后 conserved-motif 审查剔除 | 2 |
| 最终 full-length MSA 输入 | 191（旧序列 24 + 新序列 167） |
| exact focus 重复剔除 | 4 |
| 最终模型序列数 | 187（旧序列 24 + 新序列 163） |
| full MSA 列数 | 596 |
| model L | 411 |
| Neff | 33.076000287977 |
| Neff/L | 0.080476886345 |
| 三个 P450 motif anchors | 187/187 |

## Intake 与人工 QC

两组各 250 条记录中有 61 个 accession 重叠；这些重叠记录的蛋白序列、比对指标和 QC 判定全部一致。按数据库优先级保留 RefSeq 记录后得到 439 个唯一 accession。

固定 QC 未改变：identity 25–95%（含端点）、query-span coverage ≥0.70、reported qcovs ≥70%、长度比 0.75–1.25、非标准氨基酸比例 ≤1%，随后按完整蛋白质序列精确去重。439 个 accession 中 169 个通过；270 个被固定数值规则剔除。

初次 MSA 后剔除两条不完整/结构域截短记录：

- `WP_453160342.1`：FASTA 注释明确标记为 partial，参考投影中缺失 C359。
- `MDT7786182.1`：参考投影中缺失 C359 及整个 heme motif window。

最终 prepare 又识别出 4 条 exact focus duplicate：

- `WP_364525319.1` → `WP_382696941.1`
- `WP_382670764.1` → `WP_427736678.1`
- `WP_388324696.1` → `WP_342668605.1`
- `WP_427759588.1` → `WP_405444257.1`

## 权重与 Neff

权重定义保持第一轮项目规则：共同观测比例至少 0.7，pairwise identity ≥0.8 为邻居，包含自身，`w_i = 1 / n_neighbors_i`。

| 来源 | 模型 N | 权重和 |
|---|---:|---:|
| reference | 1 | 0.026315789474 |
| 第一轮合格旧序列 | 23 | 17.000000000000 |
| 新增 RefSeq | 157 | 11.535795609615 |
| 新增 clusteredNR | 6 | 4.513888888889 |

新增 163 条模型序列的权重和为 16.049684498504。相对 24 条 curated-legacy control（Neff=18.000000000000），最终 Neff 净增加 15.076000287977。

相对上一版 60 条新增模型，模型 N 增加 103，Neff 增加 14.503719524033，Neff/L 相对提高 78.09%。

## QC 结论与停止条件

reference hash 匹配；L101、L244、L295 以及 E287/R290/E305/C359 均保留；最终 187 条模型序列全部通过三 anchor 检查。

Neff/L=0.080477 高于上一版 0.045188，但仍低于冻结 V1 的 0.149306，也未达到规划目标 0.3。依照手册停止条件，本次止于 MSA、weights 和 Neff，未运行 MIp/DCA。由于冻结 V1 含 45 条明确异功能 P450，两者范围不同，不能把 V1 更高的 Neff/L 解读为更好的 TamI 功能数据集。

新增 FASTA 注释仍主要是通用 P450 或 hypothetical protein；这些记录是暂定 homolog，未由现有文件证明为 TamI ortholog，也没有 BGC 邻域证据。

## 主要文件

- 最终合并 FASTA：`02_intake/combined_curated_expanded_v2.faa`
- 最终 full MSA：`03_msa/TamI_round2_expanded250_curated_v2_full_MSA.faa`
- 最终 model MSA：`04_prepare/final_curated_v2/focus_cov70.faa`
- 最终权重：`04_prepare/final_curated_v2/weights.txt`、`sequence_weights.tsv`
- 位点映射：`04_prepare/final_curated_v2/column_map.tsv`
- 汇总与审计：`05_qc/expanded_final_summary.json`、`expanded_intake_flow.tsv`、`final_prepare_exclusions.tsv`、`final_sequence_weights_with_origin.tsv`
