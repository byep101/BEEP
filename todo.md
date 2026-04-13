# BEEP 项目任务清单 (To-Do List)
## 课题：生合成遺伝子クラスタにおける共進化解析に基づくヘム酵素の機能決定残基の推定

---

## 🟢 阶段一：数据准备与预处理 (Data Foundation)
- [ ] **BGC 种子序列提取**
    - [ ] 在 MIBiG 数据库中检索目标ヘム酵素相关的 BGC 编号。
    - [ ] 从 NCBI (GenBank) 下载对应的完全基因组序列。
    - [ ] 提取目标ヘム酵素（P450, ジオキシゲナーゼ等）的原始氨基酸序列。
- [ ] **同源序列扩增 (Homolog Expansion)**
    - [ ] 使用 Jackhmmer 或 HHblits 在 UniRef90 数据库中进行多轮迭代搜索。
    - [ ] 设置合理的 E-value 阈值，筛选具有代表性的同源序列（目标 1,000+ 条）。
- [ ] **多序列比对 (MSA) 构建**
    - [ ] 使用 MAFFT 或 Clustal Omega 生成 MSA 文件。
    - [ ] **质量检查**：去除序列覆盖度低于 70% 的碎片序列。

---

## 🔵 阶段二：并行核心分析 (Core Analysis)
- [ ] **共进化分析 (Co-evolution Analysis)**
    - [ ] 配置并运行 EVcouplings 流程（使用 PLMC 统计模型）。
    - [ ] 提取 DI (Direct Information) 或耦合强度前 50 的残基对。
    - [ ] 生成共进化接触矩阵图（Contact Map）。
- [ ] **三维结构预测 (Structural Modeling)**
    - [ ] 运行 AlphaFold2 或调用 ColabFold 进行单体结构预测。
    - [ ] **评估指标**：检查 pLDDT 分数，确保活性中心附近的建模置信度 > 85。
    - [ ] 提取 PDB 格式模型文件。

---

## 🟡 阶段三：集成解析与验证 (Integration & Validation)
- [ ] **结构-进化映射 (Mapping)**
    - [ ] 编写 PyMOL 脚本，将共进化强度分数映射至 AlphaFold 模型的 B-factor 列。
    - [ ] 可视化高分残基在空间上的分布（是否聚集在基质口袋周边）。
- [ ] **FDR (功能决定残基) 推定**
    - [ ] 结合 BGC 的基质信息，识别在特定支系中保守但在全局不保守的“开关残基”。
    - [ ] 定义功能差异化的分子基础假设。
- [ ] **计算设计验证 (ProteinMPNN)**
    - [ ] 使用 ProteinMPNN 对推定残基进行“反向设计”模拟。
    - [ ] 验证在固定 FDR 残基的情况下，序列生成的收敛性。

---

## 🔴 阶段四：产出与实验对接 (Output & Handoff)
- [ ] **自动化报告生成**
    - [ ] 生成包含 3D 映射图、共进化网络图和残基列表的 PDF 报告。
- [ ] **突变实验设计**
    - [ ] 挑选 3-5 个关键突变位点（Point Mutations），拟定 wet-lab 验证方案。
- [ ] **项目归档**
    - [ ] 将解析流程（Scripts/Notebooks）上传至内部 Git 仓库。
    - [ ] 将集成数据集按 FAIR 原则进行本地备份。

---

## 📅 里程碑 (Milestones)
- [ ] **M1 (第4周)**：完成所有目标 BGC 的 MSA 构建与质量评估。
- [ ] **M2 (第8周)**：完成首批 10 组ヘム酵素的共进化与结构映射分析。
- [ ] **M3 (第12周)**：提交 FDR 推定报告及后续生化验证建议。

---
*Status: 🚀 In Progress | Priority: High | Owner: BEEP Research Team*
