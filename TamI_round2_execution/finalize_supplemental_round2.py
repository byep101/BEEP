from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "TamI_round2_supplemental_search_v1"
BASELINE = ROOT / "TamI_round2_results" / "02_intake" / "legacy_curation_v2" / "baseline_curated_unaligned.faa"
EXISTING = ROOT / "TamI_round2_expanded_results_v1" / "02_intake" / "combined_curated_expanded_v2.faa"
UNIREF = RUN / "04_uniref_tami_qc" / "uniref90_retained_after_TamI_QC.faa"
CLUSTER = RUN / "05_clusterednr_tami_qc" / "clusterednr_new_retained_after_TamI_QC.faa"
COMBINED = RUN / "05_clusterednr_tami_qc" / "combined_round2_uniref_clusterednr_provisional.faa"
PREPARE = RUN / "06_final_msa" / "prepare"
MSA_QC = RUN / "06_final_msa" / "qc"
COEV = RUN / "07_coevolution"
CLUSTER_SOURCE = RUN / "01_clusterednr_manual"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_table(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_table(path: Path, rows: list[dict[str, object]], delimiter: str = ",") -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if line.startswith(">"):
            if identifier is not None:
                records.append((identifier, "".join(chunks).upper()))
            identifier = line[1:].split()[0]
            chunks = []
        elif line.strip():
            if identifier is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            chunks.append("".join(line.split()).replace("-", "").replace(".", ""))
    if identifier is not None:
        records.append((identifier, "".join(chunks).upper()))
    if len(records) != len({identifier for identifier, _ in records}):
        raise ValueError(f"Duplicate FASTA ID in {path}")
    return records


def verify_download_column_semantics() -> dict[str, object]:
    downloaded: dict[str, list[str]] = {}
    with (CLUSTER_SOURCE / "TamI_R2_representative_ClusteredNR_hits.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle):
            if len(row) == 13 and row[0].startswith("TamI_branch_01|"):
                downloaded[row[1]] = row
    local: dict[str, list[str]] = {}
    with (CLUSTER_SOURCE / "clusterednr250_column_semantics_check.tsv").open("r", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 11:
                local[row[0]] = row
    shared = sorted(set(downloaded) & set(local))
    result = {
        "branch01_downloaded_rows": len(downloaded),
        "matched_local_alignments": len(shared),
        "identity_matches_within_0.001": sum(abs(float(downloaded[key][2]) - float(local[key][1])) <= 0.001 for key in shared),
        "downloaded_field13_matches_local_ppos_within_0.01": sum(abs(float(downloaded[key][12]) - float(local[key][3])) <= 0.01 for key in shared),
        "downloaded_field13_matches_local_qcovs_within_0.01": sum(abs(float(downloaded[key][12]) - float(local[key][2])) <= 0.01 for key in shared),
        "median_absolute_field13_minus_qcovs": median(abs(float(downloaded[key][12]) - float(local[key][2])) for key in shared),
        "conclusion": "Downloaded column 13 is ppos (percent positives), not qcovs. Intake coverage must use recomputed TamI-reference BLAST qcovs.",
        "thresholds_changed": False,
        "source_files_changed": False,
    }
    if result["branch01_downloaded_rows"] != 250 or result["downloaded_field13_matches_local_ppos_within_0.01"] != 250:
        raise AssertionError(f"Column-semantics proof failed: {result}")
    return result


def main() -> None:
    baseline_ids = {identifier for identifier, _ in read_fasta(BASELINE)}
    existing_ids = {identifier for identifier, _ in read_fasta(EXISTING)}
    uniref_ids = {identifier for identifier, _ in read_fasta(UNIREF)}
    cluster_ids = {identifier for identifier, _ in read_fasta(CLUSTER)}
    combined_ids = [identifier for identifier, _ in read_fasta(COMBINED)]
    if len(baseline_ids) != 24 or len(existing_ids) != 191 or len(uniref_ids) != 119 or len(cluster_ids) != 152:
        raise AssertionError("Unexpected source counts")
    if set(combined_ids) != existing_ids | uniref_ids | cluster_ids or len(combined_ids) != 462:
        raise AssertionError("Combined input membership mismatch")
    if (existing_ids & uniref_ids) or (existing_ids & cluster_ids) or (uniref_ids & cluster_ids):
        raise AssertionError("Source ID groups are not disjoint")

    input_qc = json.loads((PREPARE / "input_qc.json").read_text(encoding="utf-8"))
    sequence_audit = read_table(PREPARE / "sequence_audit.tsv")
    weights = read_table(PREPARE / "sequence_weights.tsv")
    motif = read_table(MSA_QC / "motif_alignment_audit.tsv")
    column_map = read_table(PREPARE / "column_map.tsv")
    if input_qc["input_rows"] != 462 or input_qc["model_rows"] != 457 or input_qc["model_length"] != 401:
        raise AssertionError("Unexpected final MSA/prepare dimensions")
    if len(weights) != 457 or not math.isclose(sum(float(row["weight"]) for row in weights), input_qc["neff"], rel_tol=1e-12):
        raise AssertionError("Weight/Neff mismatch")
    if len(motif) != input_qc["model_rows"] or any(row["three_anchor_match"] != "True" for row in motif):
        raise AssertionError("Final model contains a three-anchor motif failure")
    excluded_focus = [row for row in sequence_audit if row["decision"] != "retained"]
    expected_focus_duplicates = {"WP_364525319.1", "WP_382670764.1", "WP_388324696.1", "WP_427759588.1", "WP_378339705.1"}
    if {row["sequence_id"] for row in excluded_focus} != expected_focus_duplicates:
        raise AssertionError(f"Unexpected prepare exclusions: {excluded_focus}")
    for position, residue in ((101, "L"), (244, "L"), (295, "L"), (287, "E"), (290, "R"), (305, "E"), (359, "C")):
        row = column_map[position - 1]
        if int(row["reference_position"]) != position or row["reference_aa"] != residue or row["status"] != "retained":
            raise AssertionError(f"Reference mapping failure at {position}: {row}")

    def origin(identifier: str) -> str:
        if identifier in baseline_ids:
            return "curated_legacy"
        if identifier in existing_ids:
            return "round2_existing_new"
        if identifier in uniref_ids:
            return "UniRef90_supplement"
        if identifier in cluster_ids:
            return "ClusteredNR_qcov_correction"
        raise KeyError(identifier)

    input_counts = Counter(origin(identifier) for identifier in combined_ids)
    model_counts: Counter[str] = Counter()
    neff_sums: Counter[str] = Counter()
    weight_values: dict[str, list[float]] = {}
    for row in weights:
        group = origin(row["sequence_id"])
        weight = float(row["weight"])
        model_counts[group] += 1
        neff_sums[group] += weight
        weight_values.setdefault(group, []).append(weight)
    source_order = ["curated_legacy", "round2_existing_new", "UniRef90_supplement", "ClusteredNR_qcov_correction"]
    neff_rows: list[dict[str, object]] = []
    for group in source_order:
        neff_rows.append({
            "source_group": group,
            "input_sequences": input_counts[group],
            "model_sequences": model_counts[group],
            "excluded_exact_focus_duplicates": input_counts[group] - model_counts[group],
            "Neff_contribution_sum_weights": neff_sums[group],
            "fraction_of_total_Neff": neff_sums[group] / input_qc["neff"],
            "mean_model_weight": sum(weight_values[group]) / len(weight_values[group]),
            "median_model_weight": median(weight_values[group]),
        })
    write_table(RUN / "final_model_source_neff_contributions.csv", neff_rows)

    semantics = verify_download_column_semantics()
    (RUN / "clusterednr_column_semantics_audit.json").write_text(json.dumps(semantics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    mi_summary = json.loads((COEV / "mi_summary.json").read_text(encoding="utf-8"))
    dca_summary = json.loads((COEV / "dca" / "dca_summary.json").read_text(encoding="utf-8"))
    smoke_summary = json.loads((COEV / "dca_smoke" / "dca_summary.json").read_text(encoding="utf-8"))
    comparison = json.loads((COEV / "comparison_summary.json").read_text(encoding="utf-8"))
    expected_pairs = input_qc["model_length"] * (input_qc["model_length"] - 1) // 2
    if mi_summary["pair_count"] != expected_pairs or dca_summary["pairs"] != expected_pairs or smoke_summary["pairs"] != expected_pairs:
        raise AssertionError("MIp/DCA pair table incomplete")
    if digest(COEV / "MIp_pairs.csv") != mi_summary["pairs_sha256"] or digest(COEV / "dca" / "DCA_pairs.csv") != dca_summary["pairs_sha256"]:
        raise AssertionError("MIp/DCA pair hash mismatch")
    if not dca_summary["outputs_validated"] or dca_summary["alignment_sha256"] != mi_summary["alignment_sha256"] or dca_summary["weights_sha256"] != mi_summary["weights_sha256"]:
        raise AssertionError("DCA structural QC or cross-method input identity failed")

    aa_by_pos = {int(row["reference_position"]): row["reference_aa"] for row in column_map if row["status"] == "retained"}
    top_rows: list[dict[str, object]] = []
    for method, path in (("MIp", COEV / "MIp_top20_pairs.csv"), ("DCA_CN", COEV / "DCA_CN_top20_pairs.csv")):
        for rank, row in enumerate(read_table(path, ","), 1):
            i, j = int(row["reference_i"]), int(row["reference_j"])
            top_rows.append({
                "method": method,
                "rank": rank,
                "reference_i": i,
                "aa_i": aa_by_pos[i],
                "reference_j": j,
                "aa_j": aa_by_pos[j],
                "pair": f"{aa_by_pos[i]}{i}-{aa_by_pos[j]}{j}",
                "score": float(row["score"]),
            })
    write_table(COEV / "top20_pairs_with_residues.csv", top_rows)

    mi_pairs = read_table(COEV / "MIp_pairs.csv", ",")
    eligible = {(int(row["reference_i"]), int(row["reference_j"])) for row in mi_pairs if row["eligible_nonlocal_variable"] == "1"}
    mip_scores = {(int(row["reference_i"]), int(row["reference_j"])): float(row["mip"]) for row in mi_pairs if row["eligible_nonlocal_variable"] == "1"}
    dca_scores = {
        (int(row["reference_i"]), int(row["reference_j"])): float(row["cn"])
        for row in read_table(COEV / "dca" / "DCA_pairs.csv", ",")
        if (int(row["reference_i"]), int(row["reference_j"])) in eligible
    }
    benchmark_by_pos = {int(row["reference_position"]): row for row in comparison["benchmark"]}
    partner_rows: list[dict[str, object]] = []
    for position in (101, 244, 295):
        for method, scores, rank_key in (("MIp", mip_scores, "MIp_rank"), ("DCA_CN", dca_scores, "DCA_CN_rank")):
            pair, score = max(((pair, score) for pair, score in scores.items() if position in pair), key=lambda item: (item[1], item[0]))
            partner = pair[1] if pair[0] == position else pair[0]
            partner_rows.append({
                "benchmark_residue": f"{aa_by_pos[position]}{position}",
                "method": method,
                "residue_rank": benchmark_by_pos[position][rank_key],
                "best_partner": f"{aa_by_pos[partner]}{partner}",
                "best_pair": f"{aa_by_pos[pair[0]]}{pair[0]}-{aa_by_pos[pair[1]]}{pair[1]}",
                "score": score,
            })
    write_table(COEV / "known_site_best_partners.csv", partner_rows)

    cluster_summary = json.loads((RUN / "05_clusterednr_tami_qc" / "clusterednr_representative_query_QC_summary.json").read_text(encoding="utf-8"))
    uniref_summary = json.loads((RUN / "04_uniref_tami_qc" / "uniref90_TamI_reference_QC_summary.json").read_text(encoding="utf-8"))
    e305_counts = Counter(row["E305_manual_checkpoint"] for row in motif)
    summary = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "status": "round2 supplemental MSA, weighting, MIp and DCA complete; DCA optimizer status requires review",
        "representative_queries": 6,
        "input": {
            "existing_round2": 191,
            "UniRef90_new_unique": 119,
            "ClusteredNR_recovered_after_qcov_correction": 152,
            "final_full_length_input": 462,
        },
        "model": {
            "N": input_qc["model_rows"],
            "L": input_qc["model_length"],
            "Neff": input_qc["neff"],
            "Neff_over_L": input_qc["neff_over_model_length"],
            "planning_target_0.3_met": input_qc["neff_over_model_length"] >= 0.3,
            "exact_focus_duplicates_removed": len(excluded_focus),
            "three_anchor_pass": f"{len(motif)}/{len(motif)} final model sequences",
            "E305_distribution": e305_counts,
            "benchmark_positions_retained": [101, 244, 295, 287, 290, 305, 359],
        },
        "field_semantics_correction": semantics,
        "MIp": mi_summary,
        "DCA": {
            "outputs_validated": dca_summary["outputs_validated"],
            "pair_count": dca_summary["pairs"],
            "optimization_status": dca_summary["optimization_status"],
            "optimizer_converged": dca_summary["optimizer_converged"],
            "plmc_effective_samples": dca_summary["plmc_effective_samples"],
            "evcouplings_version": dca_summary["evcouplings_version"],
            "params_sha256": dca_summary["params_sha256"],
            "params_local_archive": "/home/byep101/projects/tami_dca/TamI_round2_supplemental_final/dca/TamI.params",
            "params_copied_to_Windows_package": False,
        },
        "comparison": comparison,
        "source_neff_contributions": neff_rows,
        "source_thresholds_changed": False,
        "orthology_and_BGC_context_established": False,
    }
    (RUN / "final_analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=dict) + "\n", encoding="utf-8")

    (RUN / "FIELD_SEMANTICS_ERRATUM.md").write_text(
        "# ClusteredNR hit-table field correction\n\n"
        "The downloaded 13-column tables were previously interpreted as ending in `qcovs`. Direct reconstruction with local BLAST shows that column 13 is `ppos` (percent positives): for the TamI branch, 250/250 rows match local `ppos` within 0.01, whereas 0/250 match local `qcovs` within 0.01. The median absolute difference from true qcovs is "
        f"{semantics['median_absolute_field13_minus_qcovs']:.3f} percentage points.\n\n"
        "No raw file or QC threshold was changed. True qcovs was recomputed by aligning every complete protein sequence to TamI reference ADC79647.1. This recovers 152 unique sequences that pass the frozen TamI-reference rules and were previously excluded only because percent positives had been mislabeled as coverage. Previously rejected post-MSA motif failures remain blocked from re-entry.\n",
        encoding="utf-8",
    )

    neff_lines = ["| 来源 | 输入 N | 模型 N | Neff 贡献 | 占总 Neff |", "|---|---:|---:|---:|---:|"]
    for row in neff_rows:
        neff_lines.append(
            f"| {row['source_group']} | {row['input_sequences']} | {row['model_sequences']} | {float(row['Neff_contribution_sum_weights']):.6f} | {float(row['fraction_of_total_Neff']):.1%} |"
        )
    top_lines = ["| 方法 | rank | pair | score |", "|---|---:|---|---:|"]
    for row in [row for row in top_rows if int(row["rank"]) <= 10]:
        top_lines.append(f"| {row['method']} | {row['rank']} | {row['pair']} | {float(row['score']):.6g} |")
    benchmark_lines = ["| 位点 | MIp rank | DCA rank |", "|---|---:|---:|"]
    for row in comparison["benchmark"]:
        benchmark_lines.append(f"| {row['residue']} | {row['MIp_rank']} | {row['DCA_CN_rank']} |")

    report = [
        "# TamI 第二轮补充共进化分析",
        "",
        "## 最终数据集",
        "",
        "从最终第二轮输入 191 条开始，代表 query 的 UniRef90 搜索新增 119 条完整序列唯一候选。ClusteredNR 代表 query 的 FASTA 与上一批 250 条字节相同，但字段审计发现旧流程把 percent positives 误当成 qcovs。重新对完整序列做 TamI-reference BLAST 后，恢复 152 条满足冻结阈值的序列。",
        "",
        "最终 full-length MSA 输入 462 条。prepare 去除 5 条 exact focus duplicates，得到 457 条 × 401 位点模型；457/457 条模型序列均通过 E287/R290/C359 三锚点检查。L101、L244、L295、E287、R290、E305 和 C359 均保留。",
        "",
        f"θ=0.8 权重得到 Neff={input_qc['neff']:.6f}，Neff/L={input_qc['neff_over_model_length']:.6f}，超过 0.3 规划目标。",
        "",
        "## Neff 来源贡献",
        "",
        *neff_lines,
        "",
        "贡献按最终模型中的序列权重求和，四组之和等于总 Neff。ClusteredNR 组指因 qcovs 字段纠正而首次进入模型的旧下载序列。",
        "",
        "## 字段纠正",
        "",
        "TamI 分支的 250 行逐行复算显示：下载 CSV 第 13 列与 BLAST `ppos` 250/250 一致，与 `qcovs` 0/250 一致。原始文件和 QC 阈值均未修改；最终判定统一使用本地 TamI-reference BLAST 的真实 qcovs。详情见 `FIELD_SEMANTICS_ERRATUM.md`。",
        "",
        "## MIp 与 DCA",
        "",
        f"MIp/APC 和 plmc DCA 各生成 {expected_pairs:,} 个完整位点对。DCA 使用同一份 457×401 focus alignment 和 external weights；plmc effective samples={dca_summary['plmc_effective_samples']:.1f}。DCA 文件结构、pair 完整性、有限值、坐标和哈希检查通过。",
        "",
        f"100-iteration DCA 的优化器状态为 `{dca_summary['optimization_status']}`，未标记为收敛。该状态不是 `MAXIMUMITERATION`，因此按手册未触发 200 次补跑；结果状态为 `scores_generated_optimization_review`。",
        "",
        f"MIp 与 DCA 的 top20 重叠 {comparison['comparisons'][0]['overlap_pairs']}/20，top50 重叠 {comparison['comparisons'][1]['overlap_pairs']}/50。",
        "",
        *benchmark_lines,
        "",
        "三个已知 Leu 位点均保留并可评分，但均未进入任一方法的 residue top20。排名只作描述性 benchmark，不能单独证明或否定功能重要性。",
        "",
        "### 各方法前 10 位点对",
        "",
        *top_lines,
        "",
        "## 解释范围",
        "",
        "Neff/L 已达到规划目标，使完整 MIp/DCA 技术流程可以运行。新增序列仍主要依据通用 P450 注释、TamI-reference 数值 QC 和 conserved motif；TamI orthology、BGC 邻域与具体底物功能尚未建立。DCA 优化器为 MINIMUMSTEP，因此 pair 排名应作为候选生成结果，并结合结构映射、BGC 证据和稳定性分析复核。",
        "",
        "两个 256 MB `TamI.params` 文件未复制进 Windows 包。正式 DCA 参数文件保留在 WSL，SHA-256 为 `"
        + dca_summary["params_sha256"] + "`；完整 pair 表、原始 couplings、迭代记录和日志均已复制。",
    ]
    (RUN / "README_FINAL.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest: list[str] = []
    for path in sorted(RUN.rglob("*")):
        if not path.is_file() or path.name in {"SUPPLEMENTAL_SHA256SUMS.txt", "FINAL_SHA256SUMS.txt"}:
            continue
        manifest.append(f"{digest(path)}  {path.relative_to(RUN).as_posix()}")
    (RUN / "FINAL_SHA256SUMS.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=dict), flush=True)


if __name__ == "__main__":
    main()
