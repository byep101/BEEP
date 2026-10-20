# TamI 第二轮 homolog MSA 与 Neff 分析报告 v1

执行日期：2026-10-20

## 结论

第二轮已完成 homolog 数据核验、数值筛选、全长 MSA、reference 映射、序列权重和 Neff 计算。按照操作手册，本轮没有运行 MIp、DCA 或结构分析。

主数据集得到 84 条模型序列、411 个模型位点，Neff 为 18.5723，Neff/L 为 0.04519。该值低于冻结的第一轮 0.14931，也低于只清理旧数据的对照 0.04712，因此触发操作手册规定的停止条件。新增近缘 RefSeq 序列在 identity 阈值 θ=0.8 下高度冗余，没有提高有效样本深度。

## 输入数据与完整性检查

用户手动提供了 100 条 RefSeq 完整蛋白序列和 100 条 BLAST 命中记录。FASTA 与命中表的 accession 一一对应，没有缺失 accession、重复 accession、非标准氨基酸或完全重复的原始序列。

命中表没有表头，经逐列核验为 13 列：`qseqid, saccver, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcovs`。query 均为 ADC79647.1。

输入文件 SHA256：

- FASTA：`b7617e9d8b7715f4115d83986b609ee9b8556c7fa5103fec22167c98a954a12a`
- 命中表：`b9bd339248540a785ef29a0e6c670050f6466deff50490dc29ca118150ae7d1e`

## 筛选和数据集定义

按手册固定条件检查 identity 25–95%、query coverage ≥70%、长度比 0.75–1.25、非标准氨基酸比例 ≤1% 和精确序列重复：

- 100 条下载记录中，37 条因 identity >95% 排除。
- 63 条通过数值筛选，identity 为 77.262–94.891%，长度为 397–460 aa。
- 新增 63 条均注释为 Streptomyces 的泛化 `cytochrome P450 family protein`。
- 新增序列没有与此前排除的 45 条明确异功能旧序列完全相同。
- 所给文件没有 BGC 基因邻域或实验功能字段，因此这些记录只能称为高相似度 provisional homolog，不能称为已确认的 TamI ortholog。

旧数据按已批准的 TamI 功能聚焦方案处理：第一轮 69 条序列保持冻结；45 条带明确异功能证据的序列从第二轮候选中排除；保留 24 条，包括 reference 和 23 条功能仍未解决的旧候选。

建立两套数据集：

1. 主数据集：24 条清理后的旧候选加 63 条新 RefSeq 候选，共 87 条。
2. 近缘 RefSeq 敏感性集：reference 加 63 条新 RefSeq 候选，共 64 条，用于隔离旧候选对 Neff 的影响。

## MSA 和 reference QC

两套数据均使用 MAFFT 7.490，参数为 `--localpair --maxiterate 1000 --thread 4`，重新构建全长 MSA。主数据集的 full MSA 为 583 列。

主数据集进入 prepare 后有 3 条在 reference focus 投影上成为精确重复并被排除：

- WP_364525319.1
- WP_382670764.1
- WP_388324696.1

最终模型包含 84 条序列。reference 位点 1 和 2 的 occupancy 均为 0.2976，低于 0.7，因此被移除；其余 411 个 reference 位点进入模型。

以下检查全部通过：

- reference ID 为 `TamI_REF_ADC79647.1`。
- reference 长度为 413 aa，序列 SHA256 与第一轮一致。
- L101、L244、L295 均保留，分别映射到模型列 99、242、293。
- P450 ExxR 锚点 E287/R290、手册检查点 E305 和 heme-Cys C359 均保留。
- E287/R290/C359 三锚点在两套最终 focus alignment 中的匹配比例均为 1.0。
- model MSA、weights 和 column map 的 SHA256 校验全部通过。

## Neff 结果

| 数据集 | input N | model N | model L | Neff | Neff/L |
| --- | ---: | ---: | ---: | ---: | ---: |
| V1 冻结基线 | 69 | 69 | 384 | 57.3333 | 0.14931 |
| TamI 功能聚焦清理对照 | 24 | 24 | 382 | 18.0000 | 0.04712 |
| 第二轮主数据集 | 87 | 84 | 411 | 18.5723 | 0.04519 |
| 近缘 RefSeq 敏感性集 | 64 | 61 | 411 | 1.5723 | 0.00383 |

权重规则与第一轮保持一致：identity θ=0.8；共同观测位置不少于 0.7L；包含自身；权重为 `1 / neighbor count`。

在主数据集中：

- reference 的 neighbor count 为 37，权重为 0.02703。
- 23 条旧候选的权重总和为 17.0000。
- 60 条进入模型的新 RefSeq 候选的权重总和仅为 1.5453；单条权重中位数为 0.02。
- 新序列使 Neff 相对清理对照只增加 0.5723。
- Neff/L 相对清理对照下降 4.10%，相对冻结 V1 下降 69.73%。

冻结 V1 含大量明确异功能 P450，其较高 Neff 主要反映更广的 P450 功能和序列多样性，不能解释为更可靠的 TamI 特异共进化深度。第二轮主数据集与 V1 的组成发生了功能清理，因此跨轮差异不能只归因于新增序列。

## 停止条件与下一步数据需求

第二轮的 Neff/L 没有改善，且远低于规划目标 0.3，故在本轮停止。没有生成或解释第二轮 MIp、DCA、残基重要性或功能机制结果。

若继续扩充，优先手动收集以下数据：

1. 对 63 条通过筛选的 RefSeq accession，导出包含基因邻域的 GenBank 文件，建议覆盖目标 P450 上下游各 10 个基因或约 ±20 kb。
2. 标记邻域是否包含与 tirandamycin/TamI 类 BGC 相符的 PKS、调控、转运和修饰酶基因；记录排除理由。
3. 搜索更多序列多样性时，应使用 BGC/基因邻域约束，而不是通过放宽 filter 混入远缘异功能 P450。
4. 新一轮 FASTA 应保留 accession version，并附带 search database、日期、identity、coverage、完整注释和 BGC 证据表。

在取得这类功能和邻域证据前，本结果可用于说明“新增序列数量增加但有效深度未改善”，不应作为 TamI 残基共进化结论。
