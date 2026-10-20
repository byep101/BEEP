"""Write a truthful checkpoint without presenting controls as round-2 results."""
from pathlib import Path
import csv, datetime, hashlib, json, shutil, subprocess, sys

r=Path.home()/'projects/tami_dca/TamI_round2'
stage=Path('/mnt/c/Users/31441/Desktop/work/BEEP/TamI_round2_execution')
for p in stage.glob('*.py'): shutil.copyfile(p,r/'scripts'/p.name)
with (r/'logs/intake_tests.log').open('w') as log:
    subprocess.run([sys.executable,str(r/'scripts/test_intake.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
base=json.loads((r/'05_qc/baseline_audit_summary.json').read_text())
control=json.loads((r/'05_qc/baseline_curated_prepare_control/input_qc.json').read_text())
search=json.loads((r/'01_search/search_record.json').read_text())
search_status='running_no_completed_result' if 'returncode' not in search else ('completed' if search['returncode']==0 and search['output_bytes'] else 'failed_or_empty_result')
time=datetime.datetime.now(datetime.timezone.utc).isoformat()
status={'recorded_utc':time,'stage':'awaiting_search_and_required_download_email','search_status':search_status,'user_email_received':False,'scope_decision':'TamI-focused approved; 40 explicit other-function baseline sequences excluded','baseline_frozen_and_reproduced':True,'baseline_sources_unmodified':True,'baseline_curation_control':{'N':control['model_rows'],'L':control['model_length'],'Neff':control['neff'],'Neff_over_L':control['neff_over_model_length']},'round2_expanded_dataset_created':False,'round2_full_MSA_created':False,'round2_Neff_computed':False,'MIp_DCA_structure_run':False,'pending':['completed RefSeq search result','user valid email for NCBI download','new-hit annotation/taxonomy/motif review','full MSA, prepare, Neff and final QC']}
(r/'05_qc/run_status.json').write_text(json.dumps(status,indent=2,ensure_ascii=False)+'\n')
rows=[['dataset','N','L','Neff','Neff_over_L','status'],['V1_frozen',69,384,57.333333333333336,0.14930555555555555,'historical broad P450 set'],['V1_curation_only_control',29,387,23.0,23/387,'same inherited alignment; 40 exclusions; NOT round2 expansion'],['Round2_expanded','','','','','not computed; no completed new sequence intake']]
with (r/'05_qc/current_comparison.tsv').open('w',newline='') as f: csv.writer(f,delimiter='\t').writerows(rows)
readme=f'''TamI 第二轮执行记录（阶段性检查点）

记录时间 UTC：{time}
主目录：/home/byep101/projects/tami_dca/TamI_round2

状态
本轮尚未完成。已完成环境部署、第一轮冻结与重现、功能范围审计、用户批准的旧集清理和全长序列恢复。
NCBI RefSeq 检索状态：{search_status}。以 01_search/search_record.json 为最新检索状态依据。
未收到 NCBI 批量下载所需的用户有效邮箱，因此未执行命中蛋白下载。
尚未生成扩充后的第二轮 MSA、weights 或 Neff；没有运行 MIp、DCA 或结构分析。
空的 BLAST TSV 在任务未完成时不是零命中结果。

已完成的验证
1. 原始 WSL TamI_v1 文件未修改；focus、weights、input_qc、列映射、MIp/DCA pairs 已复制并记录 SHA256。
2. 用 WSL 已有 starter-kit 的 tami_minimal.py prepare 重现第一轮：N=69，L=384，Neff=57.333333333333336，Neff/L=0.14930555555555555。
   重现的模型 FASTA 与 weights 和第一轮文件逐字节相同。验证见 05_qc/baseline_reproduction_checks.json。
3. 既有注释表显示 40 条有明确不同功能注释；28 条只有泛 P450 等注释、仍不能确认 TamI 功能；另有 1 条 TamI reference。
   用户已选择 TamI 功能聚焦：旧基线冻结，第二轮剔除明确异功能序列。完整证据见 05_qc/baseline_function_scope_audit.tsv。
   40 条排除规则基于已有注释，不声称所有功能注释均为新的实验证实；未因缺少注释就推断其为 TamI ortholog。
4. 全部 68 条非 reference 旧序列在 focus 投影中都丢失过插入片段。直接去 gap 不能还原全长。
   已逐条核对旧 full MSA、原 unaligned FASTA 和 focus 投影，恢复同一组 69 条完整序列。
   恢复文件为 00_baseline/TamI_v1_verified_full_length.faa，逐条审计见 05_qc/baseline_full_length_audit.tsv。
5. 清理后 29 条的旧比对对照：N=29，L=387，Neff=23，Neff/L=0.059431524547803614。
   这只是剔除 40 条序列并重算 occupancy mask 的对照，不是新检索扩充后的第二轮结果。
   29 条中 28 条非 reference 暂为待验证 homolog；BGC 邻域和同功能/直系同源关系未证实。

主要文件
00_baseline/source_manifest.json                     第一轮原始文件与副本哈希
00_baseline/TamI_v1_verified_full_length.faa          核实后的原始全长序列
01_search/search_record.json                        BLAST 参数、开始/结束时间、返回码
logs/remote_blast.stderr.log                         BLAST 错误日志
02_intake/baseline_curated_v1/curation_audit.tsv       用户批准后的清理决策
02_intake/baseline_curated_v1/baseline_curated_unaligned.faa  29 条全长候选旧序列
05_qc/baseline_prepare_reproduction/                 首轮重现（不是 V2）
05_qc/baseline_curated_prepare_control/               清理对照 model、weights、column map
05_qc/current_comparison.tsv                         当前对照；未算的 V2 留空
05_qc/run_status.json                                当前阶段状态
scripts/                                            可复现脚本

环境与参数
Ubuntu 22.04，Python 3.10，MAFFT 7.490，BLAST+ 2.12.0。
工具通过 Ubuntu 软件源下载 .deb 后解压在 tools/prefix，未改动系统软件或第一轮环境。
软件包 SHA256 收录于 CHECKPOINT_SHA256SUMS.txt。版本原始输出见 05_qc/tool_versions.json。
检索：refseq_protein，Bacteria[Organism]，max_target_seqs=900，E-value<=1e-20；12 列与手册相同。
仅运行一个 remote 请求；输出空白时不重复提交。该任务设 1200 秒超时，结束状态写入 search_record.json。
数据库为检索当日远程库；TSV 不提供数据库 release，因此不得声称已冻结数据库版本。
下载过滤：25–95% identity；best-scoring HSP 的 query span/413 >=0.7；长度比0.75–1.25；非标准 AA <=1%；精确序列去重。
多 HSP 不直接相加，防止重叠区段导致 coverage 虚高；四项边界与 XML 解析测试通过，见 logs/intake_tests.log。
新命中以完整 GenBank XML 下载并保留 accession version、annotation、taxonomy、taxid 和 coded_by；需要有效邮箱。
新条目的功能/低复杂度初筛不能替代人工复核，motif 列表只作审查线索。
序列数<200使用 MAFFT --localpair --maxiterate 1000 --thread 4，否则 --auto --thread 4，记录策略变化。
prepare 保持 θ=0.8、共同观测 overlap>=0.7L、coverage/occupancy>=0.7、包含自身、w=1/neighbors；上限1000条，不截断。
L101/L244/L295、E305/R308/C359、reference hash、列映射、weights行数及Neff和必须核查。
Neff/L>=0.3 仅为规划目标；若相较冻结 V1 无改善则标记停止，不推进耦合推断。

与手册实施差异
手册提到的 TamI_round2_dataset_kit ZIP 未找到。本次复用本地已有 prepare，并以首轮输出逐字节复现验证。
collect_and_filter.py 与外层 QC 脚本为依照手册重建，不能声称来自缺失 ZIP。
程序使用 Python 标准库，因此没有安装手册列出的 Biopython/NumPy 依赖；权重计算仍使用原有实现。
基于用户已批准的功能聚焦决策，40 条异功能旧序列仅从第二轮候选中剔除，旧基线不变。
第二轮与 V1 的差异同时包含新序列、功能清理和可能的比对策略/列 mask 变化，不能将 Neff 变化全归因于扩充。

接续步骤
1. 确认 search_record.json returncode=0 且 TSV 有有效记录。失败时保留日志，按手册改用网页 BLAST 或明确版本的本地数据库。
2. 得到真实邮箱后，从 WSL 执行（MAIL 仅作 shell 变量示意，必须用用户提供的实际地址赋值）：
   R2=/home/byep101/projects/tami_dca/TamI_round2
   python3 "$R2/scripts/collect_and_filter.py" --blast-tsv "$R2/01_search/TamI_refseq_bacteria_blast.tsv" --baseline "$R2/02_intake/baseline_curated_v1/baseline_curated_unaligned.faa" --reference "$R2/00_baseline/TamI_reference.fasta" --out "$R2/02_intake/intake_v1" --email "$MAIL"
3. 审查 intake_audit.tsv 的所有新增注释、分类和 motif，再另存审定 FASTA 与审查说明。orthology/BGC 未证实须明确标注。
4. python3 "$R2/scripts/align_prepare.py" --input /真实路径/审定全长.faa --curation-note /真实路径/审查说明.txt
5. 检查 05_qc/validation_summary.json、v1_v2_comparison.tsv 和 round2_input_sha256.txt；本轮仍不运行 MIp/DCA。

技术资料
NCBI BLAST URL API：https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html
NCBI E-utilities 邮箱参数：https://www.ncbi.nlm.nih.gov/books/NBK25499/
'''
(r/'README.txt').write_text(readme,encoding='utf-8')
# Immutable-at-checkpoint manifest: mutable search/status files may change later.
manifest=[]
for p in sorted(r.rglob('*')):
    if not p.is_file() or p.name=='CHECKPOINT_SHA256SUMS.txt' or '__pycache__' in p.parts or 'prefix' in p.parts: continue
    manifest.append(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(r)))
(r/'CHECKPOINT_SHA256SUMS.txt').write_text('\n'.join(manifest)+'\n')
print(json.dumps(status,ensure_ascii=False,indent=2))
