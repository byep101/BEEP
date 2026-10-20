# TamI 第一轮—第二轮 intake 审计

## 结论

第二轮合并正确。第一轮69条中，按既定TamI功能聚焦范围判定合格的24条（参考序列1条、暂定同源序列23条）全部按accession进入第二轮合并输入，标准化后的完整蛋白质序列逐条完全一致，并全部进入主模型。没有合格旧序列漏合并，也没有旧序列在合并时发生序列替换。

第二轮“当前60条模型序列”是新增RefSeq子集；完整主模型为84条 = 第一轮合格24条 + 新增60条。第一轮其余45条因明确属于其他功能P450而被剔除，且未按accession或标准化完整序列被重新引入。每条旧序列的判定和具体证据见审计CSV。

## 核对方法与固定规则

- accession核对使用FASTA首字段；参考内部ID `TamI_REF_ADC79647.1` 标准化为accession `ADC79647.1`。
- 蛋白质序列标准化仅用于核对：转为大写并移除空白、`-` 和 `.`；未改写任何源数据。CSV记录标准化长度和SHA-256。
- 新命中保持既定数值QC：25% ≤ identity ≤ 95%，reported qcovs与query-span coverage均≥70%。未放宽阈值。
- prepare阶段保持既定精确focus去重：在TamI参考非gap位点投影后完全相同的序列只保留一个代表。

## Intake与合并结果

| 检查项 | 数量 | 结果 |
|---|---:|---|
| 第一轮冻结基线 | 69 | 已逐条核对 |
| 第一轮合格并应合并 | 24 | 24/24 accession匹配；24/24标准化完整序列完全一致；24/24进入主模型 |
| 第一轮范围剔除 | 45 | 40条由原注释明确异功能；5条由UniProt补充功能审查确认异功能 |
| 第二轮下载的新候选 | 100 | 完整序列均为下载集内唯一序列 |
| 新候选通过数值QC | 63 | 37条因identity >95%剔除；coverage无额外淘汰 |
| 新候选进入主模型 | 60 | 3条在focus投影后精确重复而剔除 |
| 第二轮完整主模型 | 84 | 24条第一轮合格序列 + 60条新增序列 |

第一轮45条范围剔除的证据标签计数：

- Steroid: 7
- UniProt补充功能审查: 5
- DoxA: 3
- Vancomycin: 3
- Biflaviolin: 2
- Mycinamicin: 2
- Pentalenolactone: 2
- Steroid;Methyl-branched lipid: 2
- Vitamin D: 2
- Azinomycin: 1
- Biotin biosynthesis: 1
- Cyclomarin: 1
- Epothilone: 1
- Erythromycin: 1
- Linalool: 1
- Mevastatin: 1
- Neocarzinostatin: 1
- Nocardicin: 1
- Pentalenic: 1
- PikC: 1
- PksS: 1
- Pulcherriminic: 1
- TylH1: 1
- deoxyerythronolide;Erythromycin: 1
- dihydromenaquinone: 1
- methoxybenzoate: 1

## 新增60条模型序列的identity分布

范围 77.262%–94.891%，中位数 80.121%，均值 81.964%。

| identity区间 | 序列数 | Neff权重和 |
|---|---:|---:|
| 75–<80% | 24 | 0.488171533832 |
| 80–<85% | 28 | 0.723812001065 |
| 85–<90% | 0 | 0.000000000000 |
| 90–95% | 8 | 0.333270202020 |

## 新增60条模型序列的覆盖率分布

reported qcovs范围 87.29%–97.32%，中位数 89.06%；query-span coverage范围 99.031%–99.758%，中位数 99.274%。

| reported qcovs区间 | 序列数 |
|---|---:|
| 70–<80% | 0 |
| 80–<90% | 39 |
| 90–<95% | 15 |
| 95–100% | 6 |

## 去重与Neff贡献

63条数值QC合格序列在完整蛋白质层面均不重复；其中3条只在投影至TamI参考focus位点后成为精确重复：

- `WP_364525319.1` → focus完全相同的保留代表：WP_382696941.1
- `WP_382670764.1` → focus完全相同的保留代表：WP_427736678.1
- `WP_388324696.1` → focus完全相同的保留代表：WP_342668605.1

进入模型的60条新增序列，其逐序列Neff权重合计为 1.545253736917549；单条权重范围 0.016666666666667–0.250000000000000，中位数 0.019230769230769。逐条权重和邻居数见CSV。

第一轮合格序列单独建模的Neff为 18.000000000000000，合并新增序列后的主模型Neff为 18.572280763944590，净增量为 0.572280763944590。新增权重和大于净增量，是因为重新计算邻居图后，参考/旧序列权重合计变化 -0.972972972972958；Neff是全局权重之和，不能把净增量机械分摊成彼此独立的单序列效应。

## 审计文件使用

CSV共169行：69条第一轮记录和100条第二轮新候选。筛选 `record_scope=round1_baseline` 且 `round2_primary_model_present=FALSE` 可得到45条旧序列剔除清单；筛选 `round2_current_60_new_model=TRUE` 可得到当前60条新增模型序列及其identity、coverage、去重状态和Neff权重。
