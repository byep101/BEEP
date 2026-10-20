from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "TamI_round2_supplemental_search_v1"
QC = RUN / "04_uniref_tami_qc"
AUDIT = QC / "uniref90_TamI_reference_QC_audit.tsv"
COLLECTION = RUN / "03_uniref_candidates" / "uniref90_collection_summary.json"
SELECTION = RUN / "00_queries" / "representative_selection.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bin_counts(values: list[float], edges: list[float]) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for low, high in zip(edges, edges[1:]):
        inclusive = high == edges[-1]
        count = sum(low <= value <= high if inclusive else low <= value < high for value in values)
        result.append({"range": f"{low:g}–{'≤' if inclusive else '<'}{high:g}", "count": count})
    return result


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    collection = json.loads(COLLECTION.read_text(encoding="utf-8"))
    rows = read_tsv(AUDIT)
    retained = [row for row in rows if row["final_decision"] == "retain_provisional_homolog"]
    excluded = [row for row in rows if row["final_decision"] == "exclude"]
    identities = [float(row["TamI_reference_identity_percent"]) for row in retained]
    qcovs = [float(row["TamI_reference_qcovs_percent"]) for row in retained]
    query_support = Counter(int(row["query_support_count"]) for row in rows)
    exact_existing = sum("exact_full_length_duplicate_of:" in row["exclusion_reasons"] for row in rows)
    exact_candidate_groups = Counter(row["sequence_sha256"] for row in rows)
    explicit_other = [row for row in rows if row["explicit_other_function_terms"]]
    numeric_only = [row for row in excluded if "exact_full_length_duplicate_of:" not in row["exclusion_reasons"]]
    identity_bins = bin_counts(identities, [25, 40, 50, 60, 70, 80, 90, 95])
    qcov_bins = bin_counts(qcovs, [70, 80, 90, 95, 100])
    representatives = read_tsv(SELECTION)

    status = {
        "stage": "UniRef90 automated branch complete; ClusteredNR representative-query branch awaiting manual download",
        "representative_queries": len(representatives),
        "uniref90": {
            **collection,
            "excluded_candidates": len(excluded),
            "retained_provisional_after_TamI_reference_QC": len(retained),
            "exact_duplicates_of_existing_round2_removed": exact_existing,
            "numeric_or_length_exclusions_without_existing_duplicate": len(numeric_only),
            "explicit_other_function_annotations_removed": len(explicit_other),
            "identity_distribution": identity_bins,
            "qcovs_distribution": qcov_bins,
            "query_support_distribution": dict(sorted(query_support.items())),
            "exact_sequence_duplicate_groups_within_uniref_candidates": sum(count > 1 for count in exact_candidate_groups.values()),
        },
        "clusterednr": {
            "status": "manual_download_required",
            "query_file": "01_clusterednr_manual/TamI_representative_queries_for_ClusteredNR.faa",
            "instructions": "01_clusterednr_manual/README_手动下载.md",
            "required_hit_table": "01_clusterednr_manual/TamI_R2_representative_ClusteredNR_hits.csv",
            "required_complete_fasta": "01_clusterednr_manual/TamI_R2_representative_ClusteredNR_complete.fasta",
        },
        "final_msa_weights_status": "deferred until ClusteredNR manual branch is received, merged, re-aligned to TamI reference and curated",
        "qc_thresholds_changed": False,
    }
    (RUN / "supplemental_search_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    representative_table = [
        "| 分支 | accession | TamI identity | 分组大小 |",
        "|---|---|---:|---:|",
    ]
    for row in representatives:
        representative_table.append(
            f"| {row['query_label']} | {row['accession']} | {float(row['project_identity_to_TamI_percent']):.2f}% | {row['cluster_size']} |"
        )
    identity_table = ["| TamI-reference identity | 新增暂存数 |", "|---|---:|"]
    identity_table += [f"| {row['range']}% | {row['count']} |" for row in identity_bins]
    qcov_table = ["| TamI-reference qcovs | 新增暂存数 |", "|---|---:|"]
    qcov_table += [f"| {row['range']}% | {row['count']} |" for row in qcov_bins]

    report = [
        "# TamI 第二轮代表同源序列补充搜索：阶段审计",
        "",
        "## 当前结论",
        "",
        f"已从可信同源群选择 {len(representatives)} 条不同距离分组的代表序列，完成 UniRef90 搜索。6 个 query 各返回 50 条，共 {collection['raw_hit_rows']} 个命中行；按 UniRef90 cluster accession 合并为 {collection['unique_uniref90_clusters']} 个唯一簇，其中 {collection['clusters_supported_by_multiple_queries']} 个被至少两个代表 query 命中。",
        "",
        f"147 个簇的 representative 完整序列均已重新与 TamI reference ADC79647.1 比对。冻结 QC 后，{len(retained)} 条为新的、完整序列唯一的暂存同源序列；{len(excluded)} 条未进入暂存集合，其中 {exact_existing} 条与当前第二轮输入完全同序列，1 条 identity 高于 95%，1 条长度比低于 0.75。未发现命中既定异功能关键词的 UniRef 描述。",
        "",
        "这 119 条仍不是已确认的 TamI ortholog。它们尚需与手动 ClusteredNR 补充结果合并后进行最终 full-length MSA、reference-coordinate motif 审查、focus 去重、θ=0.8 权重和 Neff 计算。",
        "",
        "## 代表 query",
        "",
        *representative_table,
        "",
        "分组使用项目既有模型位点距离和 UPGMA，只用于分散检索起点，不作为系统发育或功能证据。",
        "",
        "## UniRef90 的 TamI-reference QC",
        "",
        "纳入规则未改变：E-value ≤1e-20；25% ≤ identity ≤95%；query-span coverage ≥0.70；qcovs ≥70%；长度比 0.75–1.25；非标准氨基酸比例 ≤1%；排除既定异功能序列；按标准化完整蛋白质序列精确去重。辅助 query 的 identity/E-value 仅作发现来源记录。",
        "",
        f"暂存序列的 TamI-reference identity 范围为 {min(identities):.3f}%–{max(identities):.3f}%，中位数 {median(identities):.3f}%；qcovs 范围为 {min(qcovs):.1f}%–{max(qcovs):.1f}%，中位数 {median(qcovs):.1f}%。",
        "",
        *identity_table,
        "",
        *qcov_table,
        "",
        "逐条判定见 `04_uniref_tami_qc/uniref90_TamI_reference_QC_audit.tsv`；119 条暂存 FASTA 见 `04_uniref_tami_qc/uniref90_retained_after_TamI_QC.faa`。`combined_round2_plus_uniref_provisional.faa` 为 191 条当前输入加 119 条 UniRef 暂存序列的 310 条临时合集，不是最终 MSA 输入。",
        "",
        "## ClusteredNR 待办",
        "",
        "NCBI 远程 BLAST 对 ClusteredNR 的自动调用不可用或长期无返回。本目录已生成可直接上传的 multi-FASTA 与手动下载说明：",
        "",
        "- `01_clusterednr_manual/TamI_representative_queries_for_ClusteredNR.faa`",
        "- `01_clusterednr_manual/README_手动下载.md`",
        "",
        "收到 ClusteredNR hit table 和 complete-sequence FASTA 后，按同一 TamI-reference QC 合并并继续最终 MSA、motif 审查和 Neff/权重计算。",
    ]
    (RUN / "supplemental_search_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest_rows: list[str] = []
    for path in sorted(RUN.rglob("*")):
        if not path.is_file() or path.name == "SUPPLEMENTAL_SHA256SUMS.txt":
            continue
        manifest_rows.append(f"{digest(path)}  {path.relative_to(RUN).as_posix()}")
    (RUN / "SUPPLEMENTAL_SHA256SUMS.txt").write_text("\n".join(manifest_rows) + "\n", encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
