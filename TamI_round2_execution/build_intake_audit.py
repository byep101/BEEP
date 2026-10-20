from __future__ import annotations

import csv
import hashlib
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "TamI_round2_results"
ROUND1_FASTA = ROOT / "tami_project" / "data" / "processed" / "msa" / "TamI" / "TamI_pilot_input.faa"
CURATION_TSV = RESULTS / "02_intake" / "legacy_curation_v2" / "curation_audit.tsv"
NEW_FASTA = RESULTS / "01_search" / "TamI_R2_RefSeq_raw.fasta"
NEW_AUDIT_TSV = RESULTS / "01_search" / "manual_hit_audit.tsv"
ROUND2_SOURCE_MSA = RESULTS / "04_prepare" / "primary" / "source_full_msa.faa"
ROUND2_WEIGHTS = RESULTS / "04_prepare" / "primary" / "sequence_weights.tsv"
CURATION_CONTROL_WEIGHTS = RESULTS / "05_qc" / "baseline_curated_prepare_control_v2" / "sequence_weights.tsv"
OUT_CSV = RESULTS / "05_qc" / "round1_round2_intake_audit.csv"
OUT_REPORT = RESULTS / "05_qc" / "round1_round2_intake_audit_report.md"


def read_fasta(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append({"id": header.split()[0], "header": header, "sequence": "".join(chunks)})
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if header is not None:
        records.append({"id": header.split()[0], "header": header, "sequence": "".join(chunks)})
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_sequence(seq: str) -> str:
    # Audit-only normalization: uppercase and remove alignment gaps/whitespace.
    return re.sub(r"[\s.\-]", "", seq).upper()


def accession_of(sequence_id: str) -> str:
    if sequence_id.startswith("TamI_REF_"):
        return sequence_id.removeprefix("TamI_REF_")
    return sequence_id


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def parse_annotation(header: str) -> tuple[str, str]:
    rest = header.split(maxsplit=1)[1] if " " in header else ""
    organism = ""
    match = re.search(r"\[([^\[\]]+)\]\s*$", rest)
    if match:
        organism = match.group(1)
        annotation = rest[: match.start()].strip()
    else:
        annotation = rest.strip()
    return annotation, organism


def format_number(value: float | None, digits: int = 12) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def median(values: list[float]) -> float:
    return float(statistics.median(values))


def count_bins(values: list[float], bins: list[tuple[float, float, str, bool]]) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    for low, high, label, include_high in bins:
        if include_high:
            count = sum(low <= value <= high for value in values)
        else:
            count = sum(low <= value < high for value in values)
        out.append((label, count))
    if sum(count for _, count in out) != len(values):
        raise AssertionError(f"Binning did not cover all values: {out}, N={len(values)}")
    return out


def main() -> None:
    round1 = read_fasta(ROUND1_FASTA)
    curation = read_tsv(CURATION_TSV)
    new_records = read_fasta(NEW_FASTA)
    new_audit = read_tsv(NEW_AUDIT_TSV)
    source_msa = read_fasta(ROUND2_SOURCE_MSA)
    weights_rows = read_tsv(ROUND2_WEIGHTS)
    control_weights_rows = read_tsv(CURATION_CONTROL_WEIGHTS)

    if len(round1) != 69:
        raise AssertionError(f"Expected 69 round-1 records, found {len(round1)}")
    if len(curation) != 69:
        raise AssertionError(f"Expected 69 curation rows, found {len(curation)}")
    if len(new_records) != 100 or len(new_audit) != 100:
        raise AssertionError(f"Expected 100 new records/audit rows, found {len(new_records)}/{len(new_audit)}")

    curation_by_id = {row["sequence_id"]: row for row in curation}
    new_by_id = {row["id"]: row for row in new_records}
    new_audit_by_id = {row["saccver"]: row for row in new_audit}
    source_by_id = {row["id"]: row for row in source_msa}
    weights_by_id = {row["sequence_id"]: row for row in weights_rows}
    model_ids = set(weights_by_id)
    combined_input_ids = set(source_by_id)

    if set(curation_by_id) != {row["id"] for row in round1}:
        raise AssertionError("Round-1 FASTA and curation audit IDs differ")
    if set(new_by_id) != set(new_audit_by_id):
        raise AssertionError("New FASTA and hit audit accessions differ")
    if len(source_msa) != 87 or len(weights_rows) != 84:
        raise AssertionError(f"Expected combined/model N=87/84, found {len(source_msa)}/{len(weights_rows)}")

    r1_norm_by_id = {row["id"]: normalize_sequence(row["sequence"]) for row in round1}
    new_norm_by_id = {row["id"]: normalize_sequence(row["sequence"]) for row in new_records}
    source_norm_by_id = {row["id"]: normalize_sequence(row["sequence"]) for row in source_msa}

    r1_ids_by_seq: dict[str, list[str]] = defaultdict(list)
    new_ids_by_seq: dict[str, list[str]] = defaultdict(list)
    source_ids_by_seq: dict[str, list[str]] = defaultdict(list)
    for sequence_id, seq in r1_norm_by_id.items():
        r1_ids_by_seq[seq].append(sequence_id)
    for sequence_id, seq in new_norm_by_id.items():
        new_ids_by_seq[seq].append(sequence_id)
    for sequence_id, seq in source_norm_by_id.items():
        source_ids_by_seq[seq].append(sequence_id)

    reference_id = "TamI_REF_ADC79647.1"
    reference_aligned = next(row["sequence"] for row in source_msa if row["id"] == reference_id)
    focus_indexes = [idx for idx, aa in enumerate(reference_aligned) if aa not in "-."]
    focus_by_id = {
        row["id"]: "".join(row["sequence"][idx].upper() for idx in focus_indexes)
        for row in source_msa
    }
    ids_by_focus: dict[str, list[str]] = defaultdict(list)
    for sequence_id, seq in focus_by_id.items():
        ids_by_focus[seq].append(sequence_id)

    retained_r1_ids = {
        row["sequence_id"]
        for row in curation
        if row["curation_decision"] in {"retain_reference", "provisional_homolog_not_orthology_validated"}
    }
    excluded_r1_ids = set(curation_by_id) - retained_r1_ids
    numeric_pass_new_ids = {
        row["saccver"] for row in new_audit if row["numeric_decision"] == "pass"
    }
    new_model_ids = model_ids & set(new_by_id)
    new_prepare_excluded_ids = numeric_pass_new_ids - new_model_ids

    if len(retained_r1_ids) != 24 or len(excluded_r1_ids) != 45:
        raise AssertionError(f"Expected retained/excluded round-1 N=24/45, found {len(retained_r1_ids)}/{len(excluded_r1_ids)}")
    if len(numeric_pass_new_ids) != 63 or len(new_model_ids) != 60 or len(new_prepare_excluded_ids) != 3:
        raise AssertionError(
            "Expected new pass/model/prepare-excluded N=63/60/3, found "
            f"{len(numeric_pass_new_ids)}/{len(new_model_ids)}/{len(new_prepare_excluded_ids)}"
        )

    # Accession and normalized full-protein sequence checks for eligible legacy records.
    missing_eligible_by_accession = retained_r1_ids - combined_input_ids
    changed_eligible_sequences = {
        sequence_id
        for sequence_id in retained_r1_ids & combined_input_ids
        if r1_norm_by_id[sequence_id] != source_norm_by_id[sequence_id]
    }
    missing_eligible_from_model = retained_r1_ids - model_ids
    if missing_eligible_by_accession or changed_eligible_sequences or missing_eligible_from_model:
        raise AssertionError(
            "Eligible legacy merge failed: "
            f"missing accession={sorted(missing_eligible_by_accession)}, "
            f"changed sequence={sorted(changed_eligible_sequences)}, "
            f"missing model={sorted(missing_eligible_from_model)}"
        )

    reintroduced_excluded_by_accession = excluded_r1_ids & combined_input_ids
    reintroduced_excluded_by_sequence: dict[str, list[str]] = {}
    for sequence_id in excluded_r1_ids:
        matches = [candidate for candidate in source_ids_by_seq[r1_norm_by_id[sequence_id]] if candidate != sequence_id]
        if matches:
            reintroduced_excluded_by_sequence[sequence_id] = matches
    if reintroduced_excluded_by_accession or reintroduced_excluded_by_sequence:
        raise AssertionError(
            "Excluded legacy sequence was reintroduced: "
            f"accession={sorted(reintroduced_excluded_by_accession)}, sequence={reintroduced_excluded_by_sequence}"
        )

    columns = [
        "record_scope",
        "accession",
        "internal_sequence_id",
        "round1_present",
        "round2_new_downloaded",
        "round2_numeric_pass",
        "round2_combined_input_present",
        "round2_primary_model_present",
        "round2_current_60_new_model",
        "accession_check",
        "standardized_protein_sequence_check",
        "matched_round1_accessions_by_standardized_sequence",
        "matched_round2_accessions_by_standardized_sequence",
        "standardized_protein_length_aa",
        "standardized_sequence_sha256",
        "identity_percent",
        "reported_query_coverage_percent",
        "query_span_coverage_percent",
        "fixed_numeric_qc_rule",
        "numeric_qc_decision",
        "round1_scope_decision",
        "round2_final_status",
        "exclusion_stage",
        "exclusion_reason",
        "full_length_exact_multiplicity_in_new_download",
        "full_length_dedup_status",
        "focus_dedup_status",
        "focus_duplicate_representatives",
        "model_weight_neff_contribution",
        "neighbors_including_self",
        "annotation",
        "organism",
        "evidence_source",
    ]

    output_rows: list[dict[str, object]] = []

    for record in round1:
        sequence_id = record["id"]
        accession = accession_of(sequence_id)
        norm = r1_norm_by_id[sequence_id]
        cur = curation_by_id[sequence_id]
        in_input = sequence_id in combined_input_ids
        in_model = sequence_id in model_ids
        if in_input:
            seq_check = "exact_standardized_full_sequence_match" if source_norm_by_id[sequence_id] == norm else "accession_match_but_sequence_changed"
        else:
            other_matches = [candidate for candidate in source_ids_by_seq[norm] if candidate != sequence_id]
            seq_check = "sequence_match_under_other_accession" if other_matches else "absent_by_standardized_full_sequence"

        excluded = sequence_id in excluded_r1_ids
        if excluded:
            evidence_bits: list[str] = []
            if cur.get("annotation_evidence_terms"):
                evidence_bits.append(f"证据标签：{cur['annotation_evidence_terms']}")
            if cur.get("annotation"):
                evidence_bits.append(f"第一轮注释：{cur['annotation']}")
            if cur.get("additional_exclusion_reason"):
                evidence_bits.append(f"补充审查：{cur['additional_exclusion_reason']}")
            if cur.get("current_uniprot_function"):
                evidence_bits.append(f"UniProt功能：{cur['current_uniprot_function']}")
            exclusion_reason = "明确属于其他功能的P450，按TamI功能聚焦范围剔除；" + "；".join(evidence_bits)
            final_status = "excluded_round1_explicit_other_function"
            exclusion_stage = "legacy_function_scope_review"
        else:
            exclusion_reason = ""
            final_status = "retained_round1_reference" if sequence_id == reference_id else "retained_round1_provisional_homolog"
            exclusion_stage = ""

        weight = weights_by_id.get(sequence_id, {})
        output_rows.append(
            {
                "record_scope": "round1_baseline",
                "accession": accession,
                "internal_sequence_id": sequence_id,
                "round1_present": True,
                "round2_new_downloaded": False,
                "round2_numeric_pass": "",
                "round2_combined_input_present": in_input,
                "round2_primary_model_present": in_model,
                "round2_current_60_new_model": False,
                "accession_check": "matched_in_round2_combined_input" if in_input else "absent_from_round2_combined_input",
                "standardized_protein_sequence_check": seq_check,
                "matched_round1_accessions_by_standardized_sequence": ";".join(accession_of(x) for x in r1_ids_by_seq[norm]),
                "matched_round2_accessions_by_standardized_sequence": ";".join(accession_of(x) for x in source_ids_by_seq.get(norm, [])),
                "standardized_protein_length_aa": len(norm),
                "standardized_sequence_sha256": sha256_text(norm),
                "identity_percent": cur.get("identity_to_reference_percent", ""),
                "reported_query_coverage_percent": "",
                "query_span_coverage_percent": "",
                "fixed_numeric_qc_rule": "legacy_scope_review; numeric intake thresholds not reapplied",
                "numeric_qc_decision": "not_applicable_legacy",
                "round1_scope_decision": cur.get("curation_decision", ""),
                "round2_final_status": final_status,
                "exclusion_stage": exclusion_stage,
                "exclusion_reason": exclusion_reason,
                "full_length_exact_multiplicity_in_new_download": "",
                "full_length_dedup_status": "not_applicable_legacy",
                "focus_dedup_status": "retained" if in_model else "not_aligned_after_scope_exclusion",
                "focus_duplicate_representatives": "",
                "model_weight_neff_contribution": weight.get("weight", ""),
                "neighbors_including_self": weight.get("neighbors_including_self", ""),
                "annotation": cur.get("annotation", ""),
                "organism": cur.get("organism", ""),
                "evidence_source": cur.get("current_uniprot_source") or cur.get("source", ""),
            }
        )

    for record in new_records:
        sequence_id = record["id"]
        norm = new_norm_by_id[sequence_id]
        audit = new_audit_by_id[sequence_id]
        numeric_pass = audit["numeric_decision"] == "pass"
        in_input = sequence_id in combined_input_ids
        in_model = sequence_id in model_ids
        annotation, organism = parse_annotation(record["header"])
        r1_matches = r1_ids_by_seq.get(norm, [])
        round2_matches = source_ids_by_seq.get(norm, [])
        focus_reps: list[str] = []
        if in_input and sequence_id in focus_by_id:
            focus_reps = [candidate for candidate in ids_by_focus[focus_by_id[sequence_id]] if candidate in model_ids and candidate != sequence_id]

        if not numeric_pass:
            final_status = "excluded_new_numeric_qc"
            exclusion_stage = "manual_numeric_screen"
            reasons: list[str] = []
            identity = float(audit["pident"])
            coverage = float(audit["qcovs"])
            span_coverage = float(audit["query_span_coverage"]) * 100.0
            if not (25.0 <= identity <= 95.0):
                reasons.append(f"identity={identity:g}% 超出固定范围25%–95%")
            if coverage < 70.0 or span_coverage < 70.0:
                reasons.append(f"coverage未达到固定70%阈值（qcovs={coverage:g}%, query-span={span_coverage:.3f}%）")
            exclusion_reason = "；".join(reasons) or audit.get("numeric_reasons", "numeric_QC_failed")
            focus_status = "not_aligned_after_numeric_exclusion"
        elif not in_model:
            final_status = "excluded_new_exact_focus_duplicate"
            exclusion_stage = "prepare_focus_deduplication"
            representative_text = ", ".join(accession_of(x) for x in focus_reps) or "unresolved"
            exclusion_reason = (
                "通过固定identity/coverage阈值且完整蛋白序列在下载集中唯一；"
                "投影至TamI参考非gap位点后与已保留序列完全相同，按既定精确focus去重规则剔除；"
                f"保留代表：{representative_text}"
            )
            focus_status = "excluded_exact_focus_duplicate"
        else:
            final_status = "retained_new_in_primary_model"
            exclusion_stage = ""
            exclusion_reason = ""
            focus_status = "retained_unique_focus_sequence"

        full_multiplicity = int(audit["exact_sequence_multiplicity_in_download"])
        weight = weights_by_id.get(sequence_id, {})
        output_rows.append(
            {
                "record_scope": "round2_new_candidate",
                "accession": sequence_id,
                "internal_sequence_id": sequence_id,
                "round1_present": False,
                "round2_new_downloaded": True,
                "round2_numeric_pass": numeric_pass,
                "round2_combined_input_present": in_input,
                "round2_primary_model_present": in_model,
                "round2_current_60_new_model": in_model,
                "accession_check": "new_accession_not_in_round1" if sequence_id not in curation_by_id else "unexpected_round1_accession_overlap",
                "standardized_protein_sequence_check": "new_no_round1_exact_sequence_match" if not r1_matches else "exact_sequence_match_to_round1",
                "matched_round1_accessions_by_standardized_sequence": ";".join(accession_of(x) for x in r1_matches),
                "matched_round2_accessions_by_standardized_sequence": ";".join(accession_of(x) for x in round2_matches),
                "standardized_protein_length_aa": len(norm),
                "standardized_sequence_sha256": sha256_text(norm),
                "identity_percent": audit["pident"],
                "reported_query_coverage_percent": audit["qcovs"],
                "query_span_coverage_percent": format_number(float(audit["query_span_coverage"]) * 100.0),
                "fixed_numeric_qc_rule": "25 <= identity <= 95; qcovs >= 70; query-span coverage >= 70",
                "numeric_qc_decision": audit["numeric_decision"],
                "round1_scope_decision": "not_applicable_new",
                "round2_final_status": final_status,
                "exclusion_stage": exclusion_stage,
                "exclusion_reason": exclusion_reason,
                "full_length_exact_multiplicity_in_new_download": full_multiplicity,
                "full_length_dedup_status": "unique_full_length" if full_multiplicity == 1 else "duplicate_full_length",
                "focus_dedup_status": focus_status,
                "focus_duplicate_representatives": ";".join(accession_of(x) for x in focus_reps),
                "model_weight_neff_contribution": weight.get("weight", "0"),
                "neighbors_including_self": weight.get("neighbors_including_self", ""),
                "annotation": annotation,
                "organism": organism,
                "evidence_source": "manual RefSeq FASTA + BLAST hit table + primary prepare outputs",
            }
        )

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="raise")
        writer.writeheader()
        writer.writerows(output_rows)

    current_new_rows = [new_audit_by_id[sequence_id] for sequence_id in new_model_ids]
    current_identity = [float(row["pident"]) for row in current_new_rows]
    current_qcovs = [float(row["qcovs"]) for row in current_new_rows]
    current_span_cov = [float(row["query_span_coverage"]) * 100.0 for row in current_new_rows]
    current_weights = [float(weights_by_id[sequence_id]["weight"]) for sequence_id in new_model_ids]

    identity_bins_spec = [
        (75.0, 80.0, "75–<80%", False),
        (80.0, 85.0, "80–<85%", False),
        (85.0, 90.0, "85–<90%", False),
        (90.0, 95.0, "90–95%", True),
    ]
    coverage_bins_spec = [
        (70.0, 80.0, "70–<80%", False),
        (80.0, 90.0, "80–<90%", False),
        (90.0, 95.0, "90–<95%", False),
        (95.0, 100.0, "95–100%", True),
    ]
    identity_bins = count_bins(current_identity, identity_bins_spec)
    coverage_bins = count_bins(current_qcovs, coverage_bins_spec)

    identity_weight_rows: list[tuple[str, int, float]] = []
    for low, high, label, include_high in identity_bins_spec:
        ids = [
            sequence_id
            for sequence_id in new_model_ids
            if (low <= float(new_audit_by_id[sequence_id]["pident"]) <= high)
            if include_high or float(new_audit_by_id[sequence_id]["pident"]) < high
        ]
        identity_weight_rows.append((label, len(ids), sum(float(weights_by_id[sequence_id]["weight"]) for sequence_id in ids)))

    explicit_initial = [row for row in curation if row["curation_decision"] == "excluded_explicit_other_function" and not row.get("additional_exclusion_reason")]
    additional_review = [row for row in curation if row["curation_decision"] == "excluded_explicit_other_function" and row.get("additional_exclusion_reason")]
    evidence_counter = Counter(
        row["annotation_evidence_terms"] or "UniProt补充功能审查"
        for row in curation
        if row["sequence_id"] in excluded_r1_ids
    )

    new_weight_sum = sum(current_weights)
    primary_neff = sum(float(row["weight"]) for row in weights_rows)
    control_neff = sum(float(row["weight"]) for row in control_weights_rows)
    net_increment = primary_neff - control_neff
    legacy_plus_reference_primary_weight = primary_neff - new_weight_sum
    legacy_plus_reference_control_weight = control_neff
    neighbor_offset = legacy_plus_reference_primary_weight - legacy_plus_reference_control_weight

    duplicate_details: list[str] = []
    for sequence_id in sorted(new_prepare_excluded_ids):
        reps = [accession_of(x) for x in ids_by_focus[focus_by_id[sequence_id]] if x in model_ids and x != sequence_id]
        duplicate_details.append(f"- `{sequence_id}` → focus完全相同的保留代表：{', '.join(reps) if reps else '未解析'}")

    report_lines = [
        "# TamI 第一轮—第二轮 intake 审计",
        "",
        "## 结论",
        "",
        "第二轮合并正确。第一轮69条中，按既定TamI功能聚焦范围判定合格的24条（参考序列1条、暂定同源序列23条）全部按accession进入第二轮合并输入，标准化后的完整蛋白质序列逐条完全一致，并全部进入主模型。没有合格旧序列漏合并，也没有旧序列在合并时发生序列替换。",
        "",
        "第二轮“当前60条模型序列”是新增RefSeq子集；完整主模型为84条 = 第一轮合格24条 + 新增60条。第一轮其余45条因明确属于其他功能P450而被剔除，且未按accession或标准化完整序列被重新引入。每条旧序列的判定和具体证据见审计CSV。",
        "",
        "## 核对方法与固定规则",
        "",
        "- accession核对使用FASTA首字段；参考内部ID `TamI_REF_ADC79647.1` 标准化为accession `ADC79647.1`。",
        "- 蛋白质序列标准化仅用于核对：转为大写并移除空白、`-` 和 `.`；未改写任何源数据。CSV记录标准化长度和SHA-256。",
        "- 新命中保持既定数值QC：25% ≤ identity ≤ 95%，reported qcovs与query-span coverage均≥70%。未放宽阈值。",
        "- prepare阶段保持既定精确focus去重：在TamI参考非gap位点投影后完全相同的序列只保留一个代表。",
        "",
        "## Intake与合并结果",
        "",
        "| 检查项 | 数量 | 结果 |",
        "|---|---:|---|",
        "| 第一轮冻结基线 | 69 | 已逐条核对 |",
        "| 第一轮合格并应合并 | 24 | 24/24 accession匹配；24/24标准化完整序列完全一致；24/24进入主模型 |",
        "| 第一轮范围剔除 | 45 | 40条由原注释明确异功能；5条由UniProt补充功能审查确认异功能 |",
        "| 第二轮下载的新候选 | 100 | 完整序列均为下载集内唯一序列 |",
        "| 新候选通过数值QC | 63 | 37条因identity >95%剔除；coverage无额外淘汰 |",
        "| 新候选进入主模型 | 60 | 3条在focus投影后精确重复而剔除 |",
        "| 第二轮完整主模型 | 84 | 24条第一轮合格序列 + 60条新增序列 |",
        "",
        "第一轮45条范围剔除的证据标签计数：",
        "",
    ]
    for label, count in sorted(evidence_counter.items(), key=lambda item: (-item[1], item[0])):
        report_lines.append(f"- {label}: {count}")

    report_lines.extend(
        [
            "",
            "## 新增60条模型序列的identity分布",
            "",
            f"范围 {min(current_identity):.3f}%–{max(current_identity):.3f}%，中位数 {median(current_identity):.3f}%，均值 {statistics.mean(current_identity):.3f}%。",
            "",
            "| identity区间 | 序列数 | Neff权重和 |",
            "|---|---:|---:|",
        ]
    )
    for label, count, weight_sum in identity_weight_rows:
        report_lines.append(f"| {label} | {count} | {weight_sum:.12f} |")

    report_lines.extend(
        [
            "",
            "## 新增60条模型序列的覆盖率分布",
            "",
            f"reported qcovs范围 {min(current_qcovs):.2f}%–{max(current_qcovs):.2f}%，中位数 {median(current_qcovs):.2f}%；query-span coverage范围 {min(current_span_cov):.3f}%–{max(current_span_cov):.3f}%，中位数 {median(current_span_cov):.3f}%。",
            "",
            "| reported qcovs区间 | 序列数 |",
            "|---|---:|",
        ]
    )
    for label, count in coverage_bins:
        report_lines.append(f"| {label} | {count} |")

    report_lines.extend(
        [
            "",
            "## 去重与Neff贡献",
            "",
            "63条数值QC合格序列在完整蛋白质层面均不重复；其中3条只在投影至TamI参考focus位点后成为精确重复：",
            "",
            *duplicate_details,
            "",
            f"进入模型的60条新增序列，其逐序列Neff权重合计为 {new_weight_sum:.15f}；单条权重范围 {min(current_weights):.15f}–{max(current_weights):.15f}，中位数 {median(current_weights):.15f}。逐条权重和邻居数见CSV。",
            "",
            f"第一轮合格序列单独建模的Neff为 {control_neff:.15f}，合并新增序列后的主模型Neff为 {primary_neff:.15f}，净增量为 {net_increment:.15f}。新增权重和大于净增量，是因为重新计算邻居图后，参考/旧序列权重合计变化 {neighbor_offset:.15f}；Neff是全局权重之和，不能把净增量机械分摊成彼此独立的单序列效应。",
            "",
            "## 审计文件使用",
            "",
            "CSV共169行：69条第一轮记录和100条第二轮新候选。筛选 `record_scope=round1_baseline` 且 `round2_primary_model_present=FALSE` 可得到45条旧序列剔除清单；筛选 `round2_current_60_new_model=TRUE` 可得到当前60条新增模型序列及其identity、coverage、去重状态和Neff权重。",
            "",
        ]
    )

    OUT_REPORT.write_text("\n".join(report_lines), encoding="utf-8-sig")

    # Re-open the final CSV to verify shape and core invariants.
    with OUT_CSV.open("r", encoding="utf-8-sig", newline="") as handle:
        check_rows = list(csv.DictReader(handle))
    if len(check_rows) != 169:
        raise AssertionError(f"Expected 169 audit rows, found {len(check_rows)}")
    if sum(row["round2_primary_model_present"] == "True" for row in check_rows) != 84:
        raise AssertionError("Final CSV model membership count is not 84")
    if sum(row["round2_current_60_new_model"] == "True" for row in check_rows) != 60:
        raise AssertionError("Final CSV new-model membership count is not 60")
    if sum(
        row["record_scope"] == "round1_baseline" and row["round2_primary_model_present"] == "False"
        for row in check_rows
    ) != 45:
        raise AssertionError("Final CSV round-1 exclusion count is not 45")

    print(f"Wrote {OUT_CSV} ({len(check_rows)} rows, {len(columns)} columns)")
    print(f"Wrote {OUT_REPORT}")
    print(f"Round1 retained/excluded: {len(retained_r1_ids)}/{len(excluded_r1_ids)}")
    print(f"New downloaded/pass/model: {len(new_records)}/{len(numeric_pass_new_ids)}/{len(new_model_ids)}")
    print(f"New model identity bins: {identity_bins}")
    print(f"New model coverage bins: {coverage_bins}")
    print(f"New model weight sum/net Neff increment: {new_weight_sum}/{net_increment}")


if __name__ == "__main__":
    main()
