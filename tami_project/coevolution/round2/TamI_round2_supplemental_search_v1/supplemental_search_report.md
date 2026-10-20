# TamI 第二轮代表同源序列补充搜索：阶段审计

## 当前结论

已从可信同源群选择 6 条不同距离分组的代表序列，完成 UniRef90 搜索。6 个 query 各返回 50 条，共 300 个命中行；按 UniRef90 cluster accession 合并为 147 个唯一簇，其中 73 个被至少两个代表 query 命中。

147 个簇的 representative 完整序列均已重新与 TamI reference ADC79647.1 比对。冻结 QC 后，119 条为新的、完整序列唯一的暂存同源序列；28 条未进入暂存集合，其中 26 条与当前第二轮输入完全同序列，1 条 identity 高于 95%，1 条长度比低于 0.75。未发现命中既定异功能关键词的 UniRef 描述。

这 119 条仍不是已确认的 TamI ortholog。它们尚需与手动 ClusteredNR 补充结果合并后进行最终 full-length MSA、reference-coordinate motif 审查、focus 去重、θ=0.8 权重和 Neff 计算。

## 代表 query

| 分支 | accession | TamI identity | 分组大小 |
|---|---|---:|---:|
| TamI_branch_01 | ADC79647.1 | 100.00% | 129 |
| TamI_branch_02 | WP_364843215.1 | 58.29% | 24 |
| TamI_branch_03 | WP_433263681.1 | 53.22% | 5 |
| TamI_branch_04 | GDY32794.1 | 56.05% | 3 |
| TamI_branch_05 | WP_133850432.1 | 55.97% | 2 |
| TamI_branch_06 | HXP54575.1 | 56.71% | 1 |

分组使用项目既有模型位点距离和 UPGMA，只用于分散检索起点，不作为系统发育或功能证据。

## UniRef90 的 TamI-reference QC

纳入规则未改变：E-value ≤1e-20；25% ≤ identity ≤95%；query-span coverage ≥0.70；qcovs ≥70%；长度比 0.75–1.25；非标准氨基酸比例 ≤1%；排除既定异功能序列；按标准化完整蛋白质序列精确去重。辅助 query 的 identity/E-value 仅作发现来源记录。

暂存序列的 TamI-reference identity 范围为 48.020%–80.049%，中位数 52.475%；qcovs 范围为 93.0%–99.0%，中位数 97.0%。

| TamI-reference identity | 新增暂存数 |
|---|---:|
| 25–<40% | 0 |
| 40–<50% | 17 |
| 50–<60% | 100 |
| 60–<70% | 0 |
| 70–<80% | 1 |
| 80–<90% | 1 |
| 90–≤95% | 0 |

| TamI-reference qcovs | 新增暂存数 |
|---|---:|
| 70–<80% | 0 |
| 80–<90% | 0 |
| 90–<95% | 2 |
| 95–≤100% | 117 |

逐条判定见 `04_uniref_tami_qc/uniref90_TamI_reference_QC_audit.tsv`；119 条暂存 FASTA 见 `04_uniref_tami_qc/uniref90_retained_after_TamI_QC.faa`。`combined_round2_plus_uniref_provisional.faa` 为 191 条当前输入加 119 条 UniRef 暂存序列的 310 条临时合集，不是最终 MSA 输入。

## ClusteredNR 待办

NCBI 远程 BLAST 对 ClusteredNR 的自动调用不可用或长期无返回。本目录已生成可直接上传的 multi-FASTA 与手动下载说明：

- `01_clusterednr_manual/TamI_representative_queries_for_ClusteredNR.faa`
- `01_clusterednr_manual/README_手动下载.md`

收到 ClusteredNR hit table 和 complete-sequence FASTA 后，按同一 TamI-reference QC 合并并继续最终 MSA、motif 审查和 Neff/权重计算。
