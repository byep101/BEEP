"""Summarize primary and sensitivity Neff results without running MIp/DCA."""
from __future__ import annotations

import csv
import json
import math
import statistics
from collections import Counter
from pathlib import Path


ROOT = Path.home() / "projects/tami_dca/TamI_round2"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def weight_stats(run: Path, new_ids: set[str]) -> dict[str, object]:
    weights = read_tsv(run / "sequence_weights.tsv")
    audits = read_tsv(run / "sequence_audit.tsv")
    excluded = [row for row in audits if row["decision"] != "retained"]
    rows = []
    for row in weights:
        identifier = row["sequence_id"]
        origin = "manual_RefSeq" if identifier in new_ids else ("reference" if identifier == "TamI_REF_ADC79647.1" else "curated_legacy")
        rows.append({
            "sequence_id": identifier,
            "origin": origin,
            "weight": float(row["weight"]),
            "neighbors_including_self": int(row["neighbors_including_self"]),
        })
    values = [float(row["weight"]) for row in rows]
    by_origin = {}
    for origin in ("reference", "curated_legacy", "manual_RefSeq"):
        subset = [row for row in rows if row["origin"] == origin]
        by_origin[origin] = {
            "N": len(subset),
            "weight_sum": sum(float(row["weight"]) for row in subset),
            "neighbor_min": min((int(row["neighbors_including_self"]) for row in subset), default=None),
            "neighbor_max": max((int(row["neighbors_including_self"]) for row in subset), default=None),
        }
    return {
        "rows": rows,
        "N": len(rows),
        "Neff": sum(values),
        "weight_min": min(values),
        "weight_median": statistics.median(values),
        "weight_max": max(values),
        "by_origin": by_origin,
        "excluded_status_counts": dict(Counter(row["decision"] for row in excluded)),
        "excluded_sequence_ids": [row["sequence_id"] for row in excluded],
    }


def main() -> None:
    qc_dir = ROOT / "05_qc"
    primary_run = ROOT / "04_prepare/TamI_round2_prepare_v1"
    sensitivity_run = ROOT / "04_prepare/TamI_round2_close_refseq_prepare_v1"
    primary_qc = json.loads((primary_run / "input_qc.json").read_text())
    sensitivity_qc = json.loads((sensitivity_run / "input_qc.json").read_text())
    v1 = json.loads((ROOT / "00_baseline/TamI_v1_input_qc.json").read_text())
    control = json.loads((qc_dir / "baseline_curated_prepare_control_v2/input_qc.json").read_text())
    intake = read_tsv(ROOT / "02_intake/manual_intake_v1/manual_intake_audit.tsv")
    new_ids = {row["sequence_id"] for row in intake if row["decision"] == "retain_provisional_close_homolog"}
    primary_stats = weight_stats(primary_run, new_ids)
    sensitivity_stats = weight_stats(sensitivity_run, new_ids)
    assert math.isclose(primary_stats["Neff"], primary_qc["neff"], rel_tol=1e-12)
    assert math.isclose(sensitivity_stats["Neff"], sensitivity_qc["neff"], rel_tol=1e-12)

    comparison = [
        {"dataset": "V1_frozen", "input_N": v1["N_input"], "model_N": v1["N_model"], "model_L": v1["L_model"], "Neff": v1["Neff"], "Neff_over_L": v1["Neff_over_L_model"], "interpretation": "historical broad-P450 baseline"},
        {"dataset": "V1_TamI_focused_curation_control", "input_N": control["input_rows"], "model_N": control["model_rows"], "model_L": control["model_length"], "Neff": control["neff"], "Neff_over_L": control["neff_over_model_length"], "interpretation": "45 explicit other-function records removed; inherited alignment"},
        {"dataset": "R2_primary", "input_N": primary_qc["input_rows"], "model_N": primary_qc["model_rows"], "model_L": primary_qc["model_length"], "Neff": primary_qc["neff"], "Neff_over_L": primary_qc["neff_over_model_length"], "interpretation": "curated legacy plus 63 numeric-pass manual RefSeq homologs"},
        {"dataset": "R2_close_RefSeq_sensitivity", "input_N": sensitivity_qc["input_rows"], "model_N": sensitivity_qc["model_rows"], "model_L": sensitivity_qc["model_length"], "Neff": sensitivity_qc["neff"], "Neff_over_L": sensitivity_qc["neff_over_model_length"], "interpretation": "reference plus 63 numeric-pass manual RefSeq homologs"},
    ]
    write_tsv(qc_dir / "round2_final_comparison.tsv", comparison)
    write_tsv(qc_dir / "round2_primary_weights_by_origin.tsv", primary_stats["rows"])
    duplicates = []
    for analysis, stats in (("primary", primary_stats), ("close_refseq", sensitivity_stats)):
        for identifier in stats["excluded_sequence_ids"]:
            duplicates.append({"analysis": analysis, "sequence_id": identifier, "prepare_decision": "excluded_exact_focus_duplicate_or_low_coverage"})
    write_tsv(qc_dir / "round2_prepare_exclusions.tsv", duplicates)

    v1_ratio = v1["Neff_over_L_model"]
    primary_ratio = primary_qc["neff_over_model_length"]
    control_ratio = control["neff_over_model_length"]
    summary = {
        "scope": "homolog expansion, full-length MSA, project-defined weights and Neff only",
        "manual_input": {"downloaded_N": 100, "numeric_pass_N": 63, "identity_gt_95_excluded_N": 37},
        "primary": {
            "input_N": primary_qc["input_rows"],
            "model_N": primary_qc["model_rows"],
            "model_L": primary_qc["model_length"],
            "Neff": primary_qc["neff"],
            "Neff_over_L": primary_ratio,
            "motif_anchor_fraction": 1.0,
            "weight_details": {key: value for key, value in primary_stats.items() if key != "rows"},
        },
        "close_refseq_sensitivity": {
            "input_N": sensitivity_qc["input_rows"],
            "model_N": sensitivity_qc["model_rows"],
            "model_L": sensitivity_qc["model_length"],
            "Neff": sensitivity_qc["neff"],
            "Neff_over_L": sensitivity_qc["neff_over_model_length"],
            "weight_details": {key: value for key, value in sensitivity_stats.items() if key != "rows"},
        },
        "frozen_v1": {"N": v1["N_model"], "L": v1["L_model"], "Neff": v1["Neff"], "Neff_over_L": v1_ratio},
        "curation_only_control": {"N": control["model_rows"], "L": control["model_length"], "Neff": control["neff"], "Neff_over_L": control_ratio},
        "primary_Neff_increment_over_curation_control": primary_qc["neff"] - control["neff"],
        "primary_Neff_over_L_change_vs_frozen_v1_percent": 100 * (primary_ratio / v1_ratio - 1),
        "primary_Neff_over_L_change_vs_curation_control_percent": 100 * (primary_ratio / control_ratio - 1),
        "planning_target_0.3_met": primary_ratio >= 0.3,
        "manual_stop_condition_met": primary_ratio <= v1_ratio,
        "stop_reason": "Neff/L did not improve relative to frozen V1; new close homologs are strongly redundant at theta=0.8",
        "MIp_DCA_structure_run": False,
        "functional_limit": "Supplied annotations are generic P450 labels; TamI orthology and BGC context remain unverified",
    }
    (qc_dir / "round2_manual_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
