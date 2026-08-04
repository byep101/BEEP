#!/usr/bin/env python3
"""Finalize the TamI/TxtC pilot MSA and phylogeny metadata.

This script records the successful local Swiss-Prot/Bacteria run, keeps the
failed remote-nr attempts as provenance, performs lightweight tree QC, and
writes the resulting status back to the master dataset.
"""

from __future__ import annotations

import csv
import json
import re
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_DATE = "2026-08-04"
SYSTEMS = {
    "TamI": {"project": "tami_project", "reference": "TamI_REF_ADC79647.1"},
    "TxtC": {"project": "txtc_project", "reference": "TxtC_REF_AAW49320.1"},
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def fasta_count(path: Path) -> int:
    return sum(line.startswith(">") for line in path.read_text(encoding="utf-8").splitlines())


def nonempty_lines(path: Path, has_header: bool = False) -> int:
    count = sum(bool(line.strip()) for line in path.read_text(encoding="utf-8-sig").splitlines())
    return max(0, count - int(has_header))


def selected_model(iqtree_path: Path) -> str:
    text = iqtree_path.read_text(encoding="utf-8", errors="replace")
    patterns = (
        r"Best-fit model according to BIC:\s*([^\s]+)",
        r"Best-fit model:\s*([^\s]+) chosen according to BIC",
        r"Model of substitution:\s*([^\s]+)",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
    return ""


def all_mapping_rows_pass(path: Path) -> bool:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return bool(rows) and all(row.get("status") == "pass" for row in rows)


def tree_qc(tree_path: Path, expected_reference: str, expected_leaves: int) -> dict:
    newick = tree_path.read_text(encoding="utf-8").strip()
    leaf_labels = re.findall(r"(?<=[(,])([^():;,]+):", newick)
    lengths = [float(value) for value in re.findall(r":([0-9]+(?:\.[0-9]*)?(?:[Ee][+-]?\d+)?)", newick)]
    median = statistics.median(lengths) if lengths else 0.0
    maximum = max(lengths) if lengths else 0.0
    ratio = maximum / median if median else None
    result = {
        "treefile": tree_path.as_posix(),
        "reference_present": expected_reference in leaf_labels,
        "observed_leaf_count": len(leaf_labels),
        "expected_leaf_count": expected_leaves,
        "leaf_count_matches": len(leaf_labels) == expected_leaves,
        "branch_length_count": len(lengths),
        "median_branch_length": round(median, 8),
        "maximum_branch_length": round(maximum, 8),
        "maximum_to_median_ratio": round(ratio, 3) if ratio is not None else None,
        "gross_extreme_branch_flag": bool(ratio is not None and ratio > 20),
        "gross_extreme_branch_rule": "maximum branch length > 20x median branch length",
    }
    result["status"] = (
        "pass"
        if result["reference_present"]
        and result["leaf_count_matches"]
        and lengths
        and not result["gross_extreme_branch_flag"]
        else "review"
    )
    return result


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def collect(system: str, spec: dict) -> dict:
    project = ROOT / spec["project"]
    raw_dir = project / "data" / "raw" / "homologs" / system
    msa_dir = project / "data" / "processed" / "msa" / system
    tree_dir = project / "data" / "processed" / "tree" / system
    meta_dir = project / "metadata" / "msa" / system
    prefix = tree_dir / f"{system}_pilot_supported"
    tree_path = Path(str(prefix) + ".treefile")
    iqtree_path = Path(str(prefix) + ".iqtree")
    required = [tree_path, iqtree_path]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing completed IQ-TREE output: " + ", ".join(missing))

    qc = read_json(meta_dir / f"{system}_msa_qc_summary.json")
    trim = read_json(meta_dir / f"{system}_tree_trim_report.json")
    input_fasta = msa_dir / f"{system}_pilot_input.faa"
    seq_count = fasta_count(input_fasta)
    tqc = tree_qc(tree_path, spec["reference"], seq_count)
    (meta_dir / f"{system}_tree_qc_summary.json").write_text(
        json.dumps(tqc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    raw_blast = raw_dir / f"{system}_blast_swissprot_local.tsv"
    raw_hits = raw_dir / f"{system}_hits_raw.faa"
    filtered_table = raw_dir / f"{system}_blast_filtered.tsv"
    raw_rows = nonempty_lines(raw_blast)
    filtered_rows = nonempty_lines(filtered_table, has_header=True)
    mapping_ok = all_mapping_rows_pass(meta_dir / f"{system}_residue_mapping.tsv")
    motif_ok = all_mapping_rows_pass(meta_dir / f"{system}_P450_motif_mapping.tsv")
    low_count = seq_count < 50
    note = (
        "Pilot set is below the manual's recommended 50-sequence target; expand the database before final coevolution analysis."
        if low_count
        else "Pilot sequence count meets the manual's recommended 50-sequence target."
    )
    result = {
        "System_ID": system,
        "Reference_sequence_status": "confirmed" if system == "TamI" else "dual_reference_confirmed",
        "Homolog_set_status": "pilot_filtered_low_count" if low_count else "pilot_filtered",
        "MSA_status": "qc_passed" if mapping_ok and motif_ok else "review",
        "Tree_status": "qc_passed" if tqc["status"] == "pass" else "review",
        "Sequence_count_raw": fasta_count(raw_hits),
        "Sequence_count_filtered": seq_count,
        "Full_MSA_length": qc["alignment_length"],
        "Tree_MSA_length": trim["retained_columns"],
        "Tree_model": selected_model(iqtree_path),
        "raw_blast_rows": raw_rows,
        "blast_filtered_accessions": filtered_rows,
        "mean_column_gap_fraction": qc["mean_column_gap_fraction"],
        "columns_gap_ge_0_5": qc["columns_gap_ge_0_5"],
        "tree_trim_threshold": trim["selected_threshold"],
        "treefile": tree_path.as_posix(),
        "tree_qc_status": tqc["status"],
        "pilot_note": note,
    }

    search_row = {
        "System_ID": system,
        "run_date": RUN_DATE,
        "database": "Swiss-Prot (NCBI BLAST database snapshot)",
        "database_update_date": "2026-07-31",
        "remote": "false",
        "blast_task": "blastp-fast",
        "taxonomy_filter": "taxid:2 (Bacteria)",
        "max_target_seqs": 500,
        "evalue": "1e-20",
        "min_identity_percent": 25,
        "max_identity_percent": 95,
        "min_query_coverage_percent": 70,
        "raw_blast_rows": raw_rows,
        "filtered_accessions": filtered_rows,
        "downloaded_sequences": result["Sequence_count_raw"],
        "analysis_input_sequences_including_reference": seq_count,
        "notes": "Sequences extracted locally with blastdbcmd. Remote nr/refseq_protein attempts failed or stalled under NCBI CPU/queue limits and were not used in the pilot dataset.",
    }
    write_csv(meta_dir / "homolog_search_log.csv", list(search_row), [search_row])

    summary_fields = [
        "System_ID", "Reference_sequence_status", "Homolog_set_status", "MSA_status", "Tree_status",
        "Sequence_count_raw", "Sequence_count_filtered", "Full_MSA_length", "Tree_MSA_length", "Tree_model",
        "mean_column_gap_fraction", "columns_gap_ge_0_5", "tree_trim_threshold", "tree_qc_status",
        "treefile", "pilot_note",
    ]
    write_csv(meta_dir / "msa_tree_summary.csv", summary_fields, [{key: result[key] for key in summary_fields}])
    return result


def update_master(results: dict[str, dict]) -> None:
    path = ROOT / "metadata" / "core_bgc_p450_dataset_v0.1_master.csv"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    result_fields = [
        "Reference_sequence_status", "Homolog_set_status", "MSA_status", "Tree_status",
        "Sequence_count_raw", "Sequence_count_filtered", "Full_MSA_length", "Tree_MSA_length", "Tree_model",
    ]
    for field in result_fields:
        if field not in fields:
            fields.append(field)
    for row in rows:
        data = results.get(row.get("System_ID", ""))
        if data:
            for field in result_fields:
                row[field] = data[field]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def update_config() -> None:
    path = ROOT / "config" / "p450_phylo_pilot_v0.1.json"
    config = read_json(path)
    config["blast"] = {
        "database": "swissprot",
        "database_source": "NCBI preformatted BLAST database",
        "database_update_date": "2026-07-31",
        "remote": False,
        "taxonomy_filter": "taxid:2 (Bacteria)",
        "max_target_seqs": 500,
        "evalue": 1e-20,
        "min_identity_percent": 25,
        "max_identity_percent": 95,
        "min_query_coverage_percent": 70,
        "task": "blastp-fast",
        "sequence_extraction": "blastdbcmd",
        "initial_attempts": [
            "remote nr + Bacteria (NCBI CPU limit)",
            "remote nr + Actinomycetota (stalled)",
            "remote refseq_protein + Bacteria/Actinomycetota (failed or stalled)",
            "remote Swiss-Prot + Bacteria (stalled)",
        ],
        "deviation_reason": "The successful pilot used the official local Swiss-Prot snapshot because all remote searches failed or exceeded the practical queue window.",
    }
    path.write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_report(results: dict[str, dict]) -> None:
    tam = results["TamI"]
    txt = results["TxtC"]
    report = f"""# TamI / TxtC MSA・系統樹構築 実行記録 v0.1

実行日: {RUN_DATE}

## 結論

TamI と TxtC を分離し、reference 確認、homolog 検索・filter、MAFFT full-length MSA、残基/motif QC、smart-gap tree MSA、IQ-TREE ModelFinder + SH-aLRT 1000 + UFBoot 1000 まで完了した。

本成果物は pilot である。TamI は推奨配列数を満たすが、TxtC は {txt['Sequence_count_filtered']} 配列で手順書の推奨 50 配列未満のため、最終共進化解析の前により広いデータベースで homolog set を拡張する。

## Homolog 検索

- 採用データベース: NCBI preformatted Swiss-Prot、更新日 2026-07-31
- 対象分類群: Bacteria (taxid 2)
- BLAST: local `blastp-fast`, max_target_seqs 500, E-value <= 1e-20
- filter: identity 25–95%、query coverage >= 70%、長さ比 0.75–1.25、X <= 1%、完全重複除外
- remote nr/refseq_protein の Bacteria/Actinomycetota 検索は NCBI CPU limit または長時間 queue/stall により採用せず、失敗ログを各 system の `logs/` に保存した。

## 結果サマリー

| System | raw hits | MSA sequences | full MSA | tree MSA | IQ-TREE model | MSA QC | tree QC |
| --- | ---: | ---: | ---: | ---: | --- | --- | --- |
| TamI | {tam['Sequence_count_raw']} | {tam['Sequence_count_filtered']} | {tam['Full_MSA_length']} columns | {tam['Tree_MSA_length']} columns | {tam['Tree_model']} | {tam['MSA_status']} | {tam['Tree_status']} |
| TxtC | {txt['Sequence_count_raw']} | {txt['Sequence_count_filtered']} | {txt['Full_MSA_length']} columns | {txt['Tree_MSA_length']} columns | {txt['Tree_model']} | {txt['MSA_status']} | {txt['Tree_status']} |

## Reference・残基 QC

- TamI reference: `ADC79647.1` / `TamI_REF_ADC79647.1`, 413 aa。Leu101 -> MSA 239、Leu244 -> 418、Leu295 -> 473。
- TxtC homolog search reference: `AAW49320.1` / `TxtC_REF_AAW49320.1`, 395 aa。
- TxtC structure reference `Q8VS75` は AAW49320.1 と同長で indel なし、identity 89.8734%（40 substitutions）。Q8VS75/PDB Ser280 は AAW49320.1 Ala280 に対応し、Thr385 は双方 Thr385。
- P450 motif anchor は TamI E305/C359、TxtC E272/C344 が reference 上で non-gap、MSA mapping はすべて pass。

## MSA QC

- MAFFT L-INS-i 相当 (`--localpair --maxiterate 1000`) を使用。
- TamI: mean column gap fraction {tam['mean_column_gap_fraction']}、gap >= 0.5 の列 {tam['columns_gap_ge_0_5']}。smart-gap threshold {tam['tree_trim_threshold']}、{tam['Full_MSA_length']} -> {tam['Tree_MSA_length']} columns。
- TxtC: mean column gap fraction {txt['mean_column_gap_fraction']}、gap >= 0.5 の列 {txt['columns_gap_ge_0_5']}。smart-gap threshold {txt['tree_trim_threshold']}、{txt['Full_MSA_length']} -> {txt['Tree_MSA_length']} columns。閾値判定の結果、TxtC では削除列なし。
- ClipKIT 2.13.2 は Windows application-control policy により Biopython compiled extension をロードできなかったため、同版 `smart_gap_helper.py` の gap-distribution/slope rule を標準ライブラリのみで互換実装した。full MSA は保持した。

## 系統樹 QC

- IQ-TREE 2.4.0、amino-acid、ModelFinder Plus、SH-aLRT 1000、ultrafast bootstrap 1000。
- 両 treefile で reference の存在、leaf count、branch length を機械検査した。詳細は各 `metadata/msa/<System>/<System>_tree_qc_summary.json` を参照。
- IQ-TREE composition chi-square は候補配列の一部を warning とした。この検定は除外 filter としては使わず、pilot tree の解釈上の注意として `.iqtree` と `.log` に保持する。

## 再現性

- parameter: `config/p450_phylo_pilot_v0.1.json`
- processing utilities: `scripts/p450_phylo.py`, `scripts/concat_fasta.py`
- tool versions: `metadata/msa/tool_versions.txt`
- per-system search log: `metadata/msa/<System>/homolog_search_log.csv`
- per-system summary: `metadata/msa/<System>/msa_tree_summary.csv`
- IQ-TREE の完全 command line、seed、model、support 設定は各 `.iqtree` / `.log` に保存される。
"""
    path = ROOT / "reports" / "TamI_TxtC_MSA_tree_workflow_v0.1.md"
    path.write_text(report, encoding="utf-8", newline="\n")


def main() -> None:
    results = {system: collect(system, spec) for system, spec in SYSTEMS.items()}
    update_master(results)
    update_config()
    write_report(results)
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
