# TamI 第二轮补充共进化分析

## 最终数据集

从最终第二轮输入 191 条开始，代表 query 的 UniRef90 搜索新增 119 条完整序列唯一候选。ClusteredNR 代表 query 的 FASTA 与上一批 250 条字节相同，但字段审计发现旧流程把 percent positives 误当成 qcovs。重新对完整序列做 TamI-reference BLAST 后，恢复 152 条满足冻结阈值的序列。

最终 full-length MSA 输入 462 条。prepare 去除 5 条 exact focus duplicates，得到 457 条 × 401 位点模型；457/457 条模型序列均通过 E287/R290/C359 三锚点检查。L101、L244、L295、E287、R290、E305 和 C359 均保留。

θ=0.8 权重得到 Neff=190.197493，Neff/L=0.474308，超过 0.3 规划目标。

## Neff 来源贡献

| 来源 | 输入 N | 模型 N | Neff 贡献 | 占总 Neff |
|---|---:|---:|---:|---:|
| curated_legacy | 24 | 24 | 17.030303 | 9.0% |
| round2_existing_new | 167 | 163 | 13.653415 | 7.2% |
| UniRef90_supplement | 119 | 119 | 63.318536 | 33.3% |
| ClusteredNR_qcov_correction | 152 | 151 | 96.195238 | 50.6% |

贡献按最终模型中的序列权重求和，四组之和等于总 Neff。ClusteredNR 组指因 qcovs 字段纠正而首次进入模型的旧下载序列。

## 字段纠正

TamI 分支的 250 行逐行复算显示：下载 CSV 第 13 列与 BLAST `ppos` 250/250 一致，与 `qcovs` 0/250 一致。原始文件和 QC 阈值均未修改；最终判定统一使用本地 TamI-reference BLAST 的真实 qcovs。详情见 `FIELD_SEMANTICS_ERRATUM.md`。

## MIp 与 DCA

MIp/APC 和 plmc DCA 各生成 80,200 个完整位点对。DCA 使用同一份 457×401 focus alignment 和 external weights；plmc effective samples=190.2。DCA 文件结构、pair 完整性、有限值、坐标和哈希检查通过。

100-iteration DCA 的优化器状态为 `MINIMUMSTEP`，未标记为收敛。该状态不是 `MAXIMUMITERATION`，因此按手册未触发 200 次补跑；结果状态为 `scores_generated_optimization_review`。

MIp 与 DCA 的 top20 重叠 4/20，top50 重叠 10/50。

| 位点 | MIp rank | DCA rank |
|---|---:|---:|
| L101 | 253 | 363 |
| L244 | 238 | 157 |
| L295 | 105 | 73 |

三个已知 Leu 位点均保留并可评分，但均未进入任一方法的 residue top20。排名只作描述性 benchmark，不能单独证明或否定功能重要性。

### 各方法前 10 位点对

| 方法 | rank | pair | score |
|---|---:|---|---:|
| MIp | 1 | T75-R86 | 0.315412 |
| MIp | 2 | R35-P49 | 0.272232 |
| MIp | 3 | P40-P84 | 0.25964 |
| MIp | 4 | D76-R86 | 0.252886 |
| MIp | 5 | D76-K88 | 0.252519 |
| MIp | 6 | P84-G90 | 0.244504 |
| MIp | 7 | F16-Y24 | 0.232654 |
| MIp | 8 | R278-D378 | 0.231507 |
| MIp | 9 | D76-E83 | 0.222498 |
| MIp | 10 | P40-F79 | 0.219206 |
| DCA_CN | 1 | R278-D378 | 0.101613 |
| DCA_CN | 2 | F16-Y24 | 0.099658 |
| DCA_CN | 3 | L281-G374 | 0.098202 |
| DCA_CN | 4 | V59-L322 | 0.09778 |
| DCA_CN | 5 | R35-P49 | 0.096031 |
| DCA_CN | 6 | L151-G258 | 0.089215 |
| DCA_CN | 7 | E209-E235 | 0.088609 |
| DCA_CN | 8 | A33-D311 | 0.084269 |
| DCA_CN | 9 | A62-I309 | 0.080635 |
| DCA_CN | 10 | L78-F319 | 0.077355 |

## 解释范围

Neff/L 已达到规划目标，使完整 MIp/DCA 技术流程可以运行。新增序列仍主要依据通用 P450 注释、TamI-reference 数值 QC 和 conserved motif；TamI orthology、BGC 邻域与具体底物功能尚未建立。DCA 优化器为 MINIMUMSTEP，因此 pair 排名应作为候选生成结果，并结合结构映射、BGC 证据和稳定性分析复核。

两个 256 MB `TamI.params` 文件未复制进 Windows 包。正式 DCA 参数文件保留在 WSL，SHA-256 为 `0b9eb9f5a6f255474024060527019e0a7ae546f0767bbca9761dbb6c2eb578a8`；完整 pair 表、原始 couplings、迭代记录和日志均已复制。
