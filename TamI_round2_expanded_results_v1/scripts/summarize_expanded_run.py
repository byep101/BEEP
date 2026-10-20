from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "TamI_round2_expanded_results_v1"
INTAKE = PACKAGE / "02_intake"
MSA = PACKAGE / "03_msa" / "TamI_round2_expanded250_curated_v2_full_MSA.faa"
PREPARE = PACKAGE / "04_prepare" / "final_curated_v2"
QC = PACKAGE / "05_qc"
FINAL_INPUT = INTAKE / "combined_curated_expanded_v2.faa"
CANDIDATE_AUDIT = INTAKE / "expanded_candidate_accession_audit.tsv"
CURATION_AUDIT = ROOT / "TamI_round2_results" / "02_intake" / "legacy_curation_v2" / "curation_audit.tsv"
CONTROL_QC = ROOT / "TamI_round2_results" / "05_qc" / "baseline_curated_prepare_control_v2" / "input_qc.json"
PREVIOUS_SUMMARY = ROOT / "TamI_round2_results" / "05_qc" / "round2_manual_summary.json"
MOTIF = QC / "TamI_round2_expanded250_curated_v2_motif_alignment_audit.tsv"
REFERENCE_ID = "TamI_REF_ADC79647.1"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                rows.append((header.split()[0], "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ValueError(f"Sequence before header in {path}")
            chunks.append("".join(line.split()))
    if header is not None:
        rows.append((header.split()[0], "".join(chunks).upper()))
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    input_rows = read_fasta(FINAL_INPUT)
    msa_rows = read_fasta(MSA)
    focus_rows = read_fasta(PREPARE / "focus_all.faa")
    model_rows = read_fasta(PREPARE / "focus_cov70.faa")
    weights_rows = read_tsv(PREPARE / "sequence_weights.tsv")
    sequence_audit = read_tsv(PREPARE / "sequence_audit.tsv")
    candidate_rows = read_tsv(CANDIDATE_AUDIT)
    curation_rows = read_tsv(CURATION_AUDIT)
    motif_rows = read_tsv(MOTIF)
    column_rows = read_tsv(PREPARE / "column_map.tsv")
    input_qc = json.loads((PREPARE / "input_qc.json").read_text(encoding="utf-8"))
    control_qc = json.loads(CONTROL_QC.read_text(encoding="utf-8"))
    previous = json.loads(PREVIOUS_SUMMARY.read_text(encoding="utf-8"))
    intake_summary = json.loads((INTAKE / "expanded_intake_summary.json").read_text(encoding="utf-8"))

    if len(input_rows) != 191 or len(msa_rows) != 191:
        raise AssertionError(f"Expected final input/MSA N=191/191, found {len(input_rows)}/{len(msa_rows)}")
    if len({len(sequence) for _, sequence in msa_rows}) != 1:
        raise AssertionError("Final full MSA is not rectangular")
    if {sequence_id: sequence.replace("-", "") for sequence_id, sequence in msa_rows} != dict(input_rows):
        raise AssertionError("Final full MSA does not reproduce all input sequences after gap removal")
    if len(model_rows) != 187 or len(weights_rows) != 187:
        raise AssertionError(f"Expected model/weights N=187/187, found {len(model_rows)}/{len(weights_rows)}")
    weights_by_id = {row["sequence_id"]: row for row in weights_rows}
    model_ids = set(weights_by_id)
    neff = sum(float(row["weight"]) for row in weights_rows)
    if not math.isclose(neff, float(input_qc["neff"]), rel_tol=1e-12):
        raise AssertionError("Weight sum does not equal input_qc Neff")
    if any(not 0 < float(row["weight"]) <= 1 for row in weights_rows):
        raise AssertionError("Invalid model weight")
    if any(row["three_anchor_match"] != "True" for row in motif_rows):
        raise AssertionError("Final model contains a three-anchor motif failure")

    reference_aligned = dict(msa_rows)[REFERENCE_ID]
    focus_indexes = [index for index, aa in enumerate(reference_aligned) if aa not in "-."]
    aligned_by_id = dict(msa_rows)
    projected = {
        sequence_id: "".join(sequence[index] for index in focus_indexes)
        for sequence_id, sequence in msa_rows
    }
    ids_by_projection: dict[str, list[str]] = defaultdict(list)
    for sequence_id, sequence in projected.items():
        ids_by_projection[sequence].append(sequence_id)

    excluded_prepare = [row for row in sequence_audit if row["decision"] != "retained"]
    if len(excluded_prepare) != 4 or {row["decision"] for row in excluded_prepare} != {"excluded_exact_focus_duplicate"}:
        raise AssertionError(f"Unexpected final prepare exclusions: {excluded_prepare}")
    exclusion_rows: list[dict[str, object]] = []
    for row in excluded_prepare:
        sequence_id = row["sequence_id"]
        representatives = [candidate for candidate in ids_by_projection[projected[sequence_id]] if candidate in model_ids and candidate != sequence_id]
        if not representatives:
            raise AssertionError(f"No retained focus representative for {sequence_id}")
        exclusion_rows.append(
            {
                "sequence_id": sequence_id,
                "decision": row["decision"],
                "focus_coverage": row["focus_coverage"],
                "retained_focus_representatives": ";".join(representatives),
                "reason": "exact sequence duplicate after projection to TamI reference non-gap positions",
            }
        )
    write_tsv(QC / "final_prepare_exclusions.tsv", exclusion_rows)

    candidate_by_id = {row["accession"]: row for row in candidate_rows}
    curation_by_id = {row["sequence_id"]: row for row in curation_rows}
    legacy_ids = set(curation_by_id)
    weight_detail: list[dict[str, object]] = []
    for row in weights_rows:
        sequence_id = row["sequence_id"]
        if sequence_id == REFERENCE_ID:
            origin = "reference"
            selected_source = "legacy_reference"
            source_support = "legacy_reference"
            identity = curation_by_id[sequence_id]["identity_to_reference_percent"]
            qcovs = ""
        elif sequence_id in legacy_ids:
            origin = "curated_legacy"
            selected_source = "legacy_round1"
            source_support = "legacy_round1"
            identity = curation_by_id[sequence_id]["identity_to_reference_percent"]
            qcovs = ""
        else:
            candidate = candidate_by_id[sequence_id]
            origin = "new_RefSeq" if candidate["selected_source"] == "RefSeq" else "new_clusteredNR"
            selected_source = candidate["selected_source"]
            source_support = candidate["source_support"]
            identity = candidate["pident"]
            qcovs = candidate["reported_qcovs_percent"]
        weight_detail.append(
            {
                "sequence_id": sequence_id,
                "origin": origin,
                "selected_source": selected_source,
                "source_support": source_support,
                "identity_percent": identity,
                "reported_qcovs_percent": qcovs,
                "weight": row["weight"],
                "neighbors_including_self": row["neighbors_including_self"],
            }
        )
    write_tsv(QC / "final_sequence_weights_with_origin.tsv", weight_detail)

    origin_rows: list[dict[str, object]] = []
    for origin in ("reference", "curated_legacy", "new_RefSeq", "new_clusteredNR"):
        members = [row for row in weight_detail if row["origin"] == origin]
        values = [float(row["weight"]) for row in members]
        neighbors = [int(row["neighbors_including_self"]) for row in members]
        origin_rows.append(
            {
                "origin": origin,
                "N": len(members),
                "weight_sum": sum(values),
                "weight_min": min(values),
                "weight_median": statistics.median(values),
                "weight_max": max(values),
                "neighbor_min": min(neighbors),
                "neighbor_median": statistics.median(neighbors),
                "neighbor_max": max(neighbors),
            }
        )
    write_tsv(QC / "final_weight_summary_by_origin.tsv", origin_rows)

    new_weight_sum = sum(float(row["weight"]) for row in weight_detail if str(row["origin"]).startswith("new_"))
    old_weight_sum = neff - new_weight_sum
    control_neff = float(control_qc["neff"])
    net_increment = neff - control_neff
    neighbor_offset = old_weight_sum - control_neff
    previous_neff = float(previous["primary"]["Neff"])
    previous_ratio = float(previous["primary"]["Neff_over_L"])
    current_ratio = float(input_qc["neff_over_model_length"])
    v1_neff = float(previous["frozen_v1"]["Neff"])
    v1_ratio = float(previous["frozen_v1"]["Neff_over_L"])

    benchmark = input_qc["benchmark"]
    for pos, residue in (("101", "L"), ("244", "L"), ("295", "L")):
        if benchmark[pos]["reference_aa"] != residue or benchmark[pos]["status"] != "retained":
            raise AssertionError(f"Benchmark {pos} failed")
    reference_seq = dict(input_rows)[REFERENCE_ID]
    if reference_seq[286] != "E" or reference_seq[289] != "R" or reference_seq[304] != "E" or reference_seq[358] != "C":
        raise AssertionError("Reference motif coordinates are inconsistent")

    comparison_rows = [
        {
            "dataset": "v1_frozen_broad_scope",
            "input_or_model_N": 69,
            "model_L": 384,
            "Neff": v1_neff,
            "Neff_over_L": v1_ratio,
            "notes": "historical broad P450 baseline; includes explicit other-function proteins",
        },
        {
            "dataset": "curated_legacy_control",
            "input_or_model_N": 24,
            "model_L": int(control_qc["model_length"]),
            "Neff": control_neff,
            "Neff_over_L": float(control_qc["neff_over_model_length"]),
            "notes": "TamI-focused legacy set only",
        },
        {
            "dataset": "previous_round2_60_new",
            "input_or_model_N": int(previous["primary"]["model_N"]),
            "model_L": int(previous["primary"]["model_L"]),
            "Neff": previous_neff,
            "Neff_over_L": previous_ratio,
            "notes": "reference + 23 curated legacy + 60 new model sequences",
        },
        {
            "dataset": "expanded_final_curated_v2",
            "input_or_model_N": int(input_qc["model_rows"]),
            "model_L": int(input_qc["model_length"]),
            "Neff": neff,
            "Neff_over_L": current_ratio,
            "notes": "reference + 23 curated legacy + 163 new model sequences",
        },
    ]
    write_tsv(QC / "expanded_final_comparison.tsv", comparison_rows)

    input_flow_rows = [
        {"stage": "raw source rows", "N": 500, "excluded_at_stage": 0, "notes": "250 RefSeq + 250 clusteredNR"},
        {"stage": "unique accessions", "N": 439, "excluded_at_stage": 61, "notes": "same-accession cross-source merge; RefSeq precedence"},
        {"stage": "fixed numeric QC pass", "N": 169, "excluded_at_stage": 270, "notes": "unchanged identity/coverage/length/ambiguity rules"},
        {"stage": "combined with curated legacy", "N": 193, "excluded_at_stage": 0, "notes": "169 new + 24 legacy/reference"},
        {"stage": "post-MSA motif-curated input", "N": 191, "excluded_at_stage": 2, "notes": "excluded WP_453160342.1 and MDT7786182.1"},
        {"stage": "final model after exact focus dedup", "N": 187, "excluded_at_stage": 4, "notes": "four exact focus duplicates"},
    ]
    write_tsv(QC / "expanded_intake_flow.tsv", input_flow_rows)

    summary = {
        "analysis_id": "TamI_round2_expanded250_curated_v2",
        "stage_completed": "sequence merge, full-length MSA, model masking, project-defined weights and Neff",
        "raw_sources": {"RefSeq": 250, "clusteredNR_batch01": 250},
        "unique_accessions": 439,
        "cross_source_accession_overlap": 61,
        "fixed_numeric_qc_pass_new": 169,
        "post_msa_motif_excluded": ["WP_453160342.1", "MDT7786182.1"],
        "final_combined_input_N": 191,
        "prepare_focus_duplicates_excluded": exclusion_rows,
        "final_model_N": int(input_qc["model_rows"]),
        "full_msa_columns": int(input_qc["full_msa_columns"]),
        "model_L": int(input_qc["model_length"]),
        "Neff": neff,
        "Neff_over_L": current_ratio,
        "weights_by_origin": origin_rows,
        "new_model_weight_sum": new_weight_sum,
        "net_Neff_increment_over_curated_legacy_control": net_increment,
        "legacy_reference_weight_change_after_new_neighbors": neighbor_offset,
        "change_vs_previous_round2": {
            "model_N_delta": int(input_qc["model_rows"]) - int(previous["primary"]["model_N"]),
            "Neff_delta": neff - previous_neff,
            "Neff_over_L_relative_change_percent": (current_ratio / previous_ratio - 1.0) * 100.0,
        },
        "motif_anchor_fraction": 1.0,
        "reference_hash_match": bool(input_qc["baseline_reference_hash_match"]),
        "benchmark_positions_retained": [101, 244, 295, 287, 290, 305, 359],
        "planning_target_0.3_met": current_ratio >= 0.3,
        "manual_stop_condition_met": current_ratio <= v1_ratio,
        "stop_reason": "Neff/L remains below the frozen broad-scope V1 value; do not continue to MIp/DCA under the manual",
        "MIp_DCA_run": False,
        "functional_limit": intake_summary["functional_limit"],
        "final_input_sha256": digest(FINAL_INPUT),
        "final_full_msa_sha256": digest(MSA),
        "final_focus_model_sha256": digest(PREPARE / "focus_cov70.faa"),
        "final_weights_sha256": digest(PREPARE / "weights.txt"),
    }
    (QC / "expanded_final_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    origin_by_name = {row["origin"]: row for row in origin_rows}
    duplicate_lines = [
        f"- `{row['sequence_id']}` → `{row['retained_focus_representatives']}`"
        for row in exclusion_rows
    ]
    report = [
        "# TamI 第二轮扩充（RefSeq 250 + clusteredNR 250）",
        "",
        "## 最终结果",
        "",
        "已完成序列验证、固定阈值过滤、跨库合并、完整序列 MSA、模型位点过滤、θ=0.8 权重和 Neff 计算。最终分析使用 `combined_curated_expanded_v2.faa`；初次 MSA 仅作为 motif 审查记录保留。",
        "",
        "| 指标 | 最终值 |",
        "|---|---:|",
        "| 原始来源记录 | 500（RefSeq 250 + clusteredNR 250） |",
        "| 跨库唯一 accession | 439 |",
        "| 固定数值 QC 合格的新序列 | 169 |",
        "| MSA 后 conserved-motif 审查剔除 | 2 |",
        "| 最终 full-length MSA 输入 | 191（旧序列 24 + 新序列 167） |",
        "| exact focus 重复剔除 | 4 |",
        "| 最终模型序列数 | 187（旧序列 24 + 新序列 163） |",
        f"| full MSA 列数 | {input_qc['full_msa_columns']} |",
        f"| model L | {input_qc['model_length']} |",
        f"| Neff | {neff:.12f} |",
        f"| Neff/L | {current_ratio:.12f} |",
        "| 三个 P450 motif anchors | 187/187 |",
        "",
        "## Intake 与人工 QC",
        "",
        "两组各 250 条记录中有 61 个 accession 重叠；这些重叠记录的蛋白序列、比对指标和 QC 判定全部一致。按数据库优先级保留 RefSeq 记录后得到 439 个唯一 accession。",
        "",
        "固定 QC 未改变：identity 25–95%（含端点）、query-span coverage ≥0.70、reported qcovs ≥70%、长度比 0.75–1.25、非标准氨基酸比例 ≤1%，随后按完整蛋白质序列精确去重。439 个 accession 中 169 个通过；270 个被固定数值规则剔除。",
        "",
        "初次 MSA 后剔除两条不完整/结构域截短记录：",
        "",
        "- `WP_453160342.1`：FASTA 注释明确标记为 partial，参考投影中缺失 C359。",
        "- `MDT7786182.1`：参考投影中缺失 C359 及整个 heme motif window。",
        "",
        "最终 prepare 又识别出 4 条 exact focus duplicate：",
        "",
        *duplicate_lines,
        "",
        "## 权重与 Neff",
        "",
        "权重定义保持第一轮项目规则：共同观测比例至少 0.7，pairwise identity ≥0.8 为邻居，包含自身，`w_i = 1 / n_neighbors_i`。",
        "",
        "| 来源 | 模型 N | 权重和 |",
        "|---|---:|---:|",
        f"| reference | {origin_by_name['reference']['N']} | {float(origin_by_name['reference']['weight_sum']):.12f} |",
        f"| 第一轮合格旧序列 | {origin_by_name['curated_legacy']['N']} | {float(origin_by_name['curated_legacy']['weight_sum']):.12f} |",
        f"| 新增 RefSeq | {origin_by_name['new_RefSeq']['N']} | {float(origin_by_name['new_RefSeq']['weight_sum']):.12f} |",
        f"| 新增 clusteredNR | {origin_by_name['new_clusteredNR']['N']} | {float(origin_by_name['new_clusteredNR']['weight_sum']):.12f} |",
        "",
        f"新增 163 条模型序列的权重和为 {new_weight_sum:.12f}。相对 24 条 curated-legacy control（Neff={control_neff:.12f}），最终 Neff 净增加 {net_increment:.12f}。",
        "",
        f"相对上一版 60 条新增模型，模型 N 增加 {int(input_qc['model_rows']) - int(previous['primary']['model_N'])}，Neff 增加 {neff - previous_neff:.12f}，Neff/L 相对提高 {(current_ratio / previous_ratio - 1.0) * 100.0:.2f}%。",
        "",
        "## QC 结论与停止条件",
        "",
        "reference hash 匹配；L101、L244、L295 以及 E287/R290/E305/C359 均保留；最终 187 条模型序列全部通过三 anchor 检查。",
        "",
        f"Neff/L={current_ratio:.6f} 高于上一版 0.045188，但仍低于冻结 V1 的 0.149306，也未达到规划目标 0.3。依照手册停止条件，本次止于 MSA、weights 和 Neff，未运行 MIp/DCA。由于冻结 V1 含 45 条明确异功能 P450，两者范围不同，不能把 V1 更高的 Neff/L 解读为更好的 TamI 功能数据集。",
        "",
        "新增 FASTA 注释仍主要是通用 P450 或 hypothetical protein；这些记录是暂定 homolog，未由现有文件证明为 TamI ortholog，也没有 BGC 邻域证据。",
        "",
        "## 主要文件",
        "",
        "- 最终合并 FASTA：`02_intake/combined_curated_expanded_v2.faa`",
        "- 最终 full MSA：`03_msa/TamI_round2_expanded250_curated_v2_full_MSA.faa`",
        "- 最终 model MSA：`04_prepare/final_curated_v2/focus_cov70.faa`",
        "- 最终权重：`04_prepare/final_curated_v2/weights.txt`、`sequence_weights.tsv`",
        "- 位点映射：`04_prepare/final_curated_v2/column_map.tsv`",
        "- 汇总与审计：`05_qc/expanded_final_summary.json`、`expanded_intake_flow.tsv`、`final_prepare_exclusions.tsv`、`final_sequence_weights_with_origin.tsv`",
        "",
    ]
    (PACKAGE / "README.md").write_text("\n".join(report), encoding="utf-8-sig")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
