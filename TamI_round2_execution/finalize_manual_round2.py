"""Finalize the stopped round-2 checkpoint and a reproducibility manifest."""
from __future__ import annotations

import datetime
import hashlib
import json
from pathlib import Path


ROOT = Path.home() / "projects/tami_dca/TamI_round2"
QC = ROOT / "05_qc"
summary = json.loads((QC / "round2_manual_summary.json").read_text())
status = {
    "recorded_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "stage": "round2_MSA_Neff_completed_stopped_by_manual_rule",
    "scope": "homolog intake, full-length MSA, project-defined sequence weights and Neff",
    "manual_inputs_received": True,
    "manual_input_N": 100,
    "numeric_pass_N": 63,
    "primary_input_N": summary["primary"]["input_N"],
    "primary_model_N": summary["primary"]["model_N"],
    "primary_model_L": summary["primary"]["model_L"],
    "primary_Neff": summary["primary"]["Neff"],
    "primary_Neff_over_L": summary["primary"]["Neff_over_L"],
    "stop_condition_met": True,
    "stop_reason": summary["stop_reason"],
    "planning_target_0.3_met": False,
    "reference_and_benchmark_QC_pass": True,
    "motif_anchor_QC_pass": True,
    "MIp_DCA_structure_run": False,
    "functional_scope_limit": summary["functional_limit"],
    "manual_followup_data": "GenBank genomic neighborhoods for the 63 numeric-pass accessions, preferably +/-10 genes or about +/-20 kb, with BGC membership review",
}
(QC / "run_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

files = [
    ROOT / "01_search/manual_download_v1/TamI_R2_RefSeq_raw.fasta",
    ROOT / "01_search/manual_download_v1/TamI_R2_RefSeq_hits.csv",
    ROOT / "02_intake/manual_intake_v1/combined_curated_provisional_v1.faa",
    ROOT / "02_intake/manual_intake_v1/close_refseq_sensitivity_v1.faa",
    ROOT / "03_msa/TamI_round2_full_MSA.faa",
    ROOT / "03_msa/TamI_round2_close_refseq_full_MSA.faa",
    ROOT / "04_prepare/TamI_round2_prepare_v1/focus_cov70.faa",
    ROOT / "04_prepare/TamI_round2_prepare_v1/weights.txt",
    ROOT / "04_prepare/TamI_round2_prepare_v1/column_map.tsv",
    ROOT / "04_prepare/TamI_round2_prepare_v1/input_qc.json",
    ROOT / "04_prepare/TamI_round2_close_refseq_prepare_v1/focus_cov70.faa",
    ROOT / "04_prepare/TamI_round2_close_refseq_prepare_v1/weights.txt",
    ROOT / "04_prepare/TamI_round2_close_refseq_prepare_v1/column_map.tsv",
    ROOT / "04_prepare/TamI_round2_close_refseq_prepare_v1/input_qc.json",
    QC / "round2_manual_summary.json",
    QC / "round2_final_comparison.tsv",
    ROOT / "reports/TamI_round2_MSA_Neff_report_v1.md",
]
missing = [str(path) for path in files if not path.is_file()]
if missing:
    raise FileNotFoundError("Missing final files: " + ", ".join(missing))
manifest = []
for path in files:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest.append(f"{digest}  {path.relative_to(ROOT)}")
(QC / "ROUND2_FINAL_SHA256SUMS.txt").write_text("\n".join(manifest) + "\n", encoding="ascii")
print(json.dumps(status, ensure_ascii=False, indent=2))
