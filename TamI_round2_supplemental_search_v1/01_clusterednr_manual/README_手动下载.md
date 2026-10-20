# ClusteredNR 代表序列补充搜索：手动下载说明

自动远程 BLAST 已分别测试 NCBI 数据库名 `nr_clustered` 和 `nr_cluster_seq`。前者被服务器拒绝为不存在，后者接受任务但长时间无结果，因此本分支按用户要求转为手动收集。

## 上传文件

`TamI_representative_queries_for_ClusteredNR.faa`

该文件包含 6 条 query，来自 TamI 可信同源群的不同距离分组。SHA-256：

`18b25eb77bd6d7fc079a46998b25ef18ce458df317a21dd76aedc0064e42f8b8`

## NCBI Protein BLAST 设置

1. 上传上述 multi-FASTA。
2. 数据库选择 **ClusteredNR（nr_clustered / experimental）**。
3. Organism 设为 **Bacteria (taxid:2)**。
4. Program 选择 `blastp`，matrix 为 `BLOSUM62`。
5. Expect threshold 设为 `1e-20`。
6. 关闭 low-complexity/SEG filter，与本轮其他搜索保持一致。
7. 每条 query 保留最多 50 个命中。不要预先按 TamI identity 删除结果。

## 需要下载的两个文件

完成搜索后，请下载并放到本目录：

- `TamI_R2_representative_ClusteredNR_hits.csv`：全部 6 条 query 的 Hit Table。NCBI 网页实际导出为无表头的 13 列 CSV，这是正常格式，请直接保留原文件，不要手工添加列名或改列。本次文件的第 13 列经核对为 `ppos`，不是 `qcovs`；辅助 query 的这些字段只用于记录搜索来源，正式纳入仍以 TamI reference 的本地重比对结果为准。
- `TamI_R2_representative_ClusteredNR_complete.txt`：上述命中的 **FASTA (complete sequence)**。NCBI 会把 FASTA 内容保存为 `.txt`；只要文件中的记录以 `>` 开头，就是可直接解析的标准 FASTA，无需改扩展名。不要下载 aligned regions。

如果 NCBI 页面只能按 query 分别导出，也可保留 6 组原始 CSV 和 FASTA 文本文件，文件名中写入 `TamI_branch_01` 至 `TamI_branch_06`；后续脚本会先合并 accession 和完整序列。

这些命中只是发现候选。下载后会统一与 TamI reference `ADC79647.1` 重新比对，并应用冻结规则：TamI-reference E-value ≤ `1e-20`、identity 25–95%（含端点）、query-span coverage ≥0.70、qcovs ≥70%、长度比 0.75–1.25、非标准氨基酸比例 ≤1%、已知异功能序列排除及完整蛋白序列精确去重。辅助 query 的 BLAST 指标不会直接作为纳入依据。
