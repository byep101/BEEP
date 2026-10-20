"""Build versioned primary and close-homolog sensitivity FASTAs from manual data."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from import_manual_data import read_fasta, write_tsv


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def main() -> None:
    root = Path("C:/Users/31441/Desktop/work/BEEP/tami_project/coevolution/round2/TamI_round2_manual_input")
    out = root / "intake_v1"
    out.mkdir(exist_ok=False)
    raw = {identifier: (header, sequence) for identifier, header, sequence in read_fasta(root / "TamI_R2_RefSeq_raw.fasta")}
    baseline_rows = read_fasta(root / "baseline_curated_v2_unaligned.faa")
    baseline = [(identifier, sequence) for identifier, _, sequence in baseline_rows]
    if baseline[0][0] != "TamI_REF_ADC79647.1" or len(baseline[0][1]) != 413:
        raise ValueError("Curated baseline reference is missing, not first, or not 413 aa")

    with (root / "validation/manual_hit_audit.tsv").open(encoding="utf-8", newline="") as handle:
        hits = list(csv.DictReader(handle, delimiter="\t"))
    passed = [row for row in hits if row["numeric_decision"] == "pass"]
    if len(passed) != 63:
        raise ValueError(f"Expected 63 manual-filter pass records, got {len(passed)}")

    primary = list(baseline)
    close_only = [baseline[0]]
    seen_primary = {sequence: identifier for identifier, sequence in primary}
    seen_close = {baseline[0][1]: baseline[0][0]}
    audit: list[dict[str, object]] = []
    for row in hits:
        identifier = row["saccver"]
        header, sequence = raw[identifier]
        organism_match = re.search(r"\[([^\]]+)\]$", header)
        organism = organism_match.group(1) if organism_match else ""
        annotation = header[len(identifier):].strip()
        duplicate_primary = seen_primary.get(sequence, "")
        decision = "exclude_numeric_filter"
        if row["numeric_decision"] == "pass":
            decision = "exclude_exact_sequence_duplicate" if duplicate_primary else "retain_provisional_close_homolog"
            if not duplicate_primary:
                primary.append((identifier, sequence))
                seen_primary[sequence] = identifier
            if sequence not in seen_close:
                close_only.append((identifier, sequence))
                seen_close[sequence] = identifier
        audit.append({
            "sequence_id": identifier,
            "decision": decision,
            "reason": row["numeric_reasons"],
            "identity_percent": row["pident"],
            "query_coverage_percent_reported": row["qcovs"],
            "query_span_coverage": row["query_span_coverage"],
            "length": len(sequence),
            "organism": organism,
            "annotation": annotation,
            "duplicate_of_primary": duplicate_primary,
            "orthology_status": "not established; high-similarity provisional homolog",
            "BGC_status": "not checked; genomic-neighborhood data not supplied",
        })

    if len(primary) > 1000 or len(close_only) > 1000:
        raise ValueError("Prepare limit exceeded; do not truncate")
    write_fasta(out / "combined_curated_provisional_v1.faa", primary)
    write_fasta(out / "close_refseq_sensitivity_v1.faa", close_only)
    write_fasta(out / "manual_numeric_pass_v1.faa", [(row["saccver"], raw[row["saccver"]][1]) for row in passed])
    write_tsv(out / "manual_intake_audit.tsv", audit)
    identity_values = [float(row["pident"]) for row in passed]
    summary = {
        "manual_download_records": len(raw),
        "manual_numeric_pass": len(passed),
        "manual_numeric_excluded_identity_gt_95": sum(float(row["pident"]) > 95 for row in hits),
        "manual_pass_identity_min": min(identity_values),
        "manual_pass_identity_max": max(identity_values),
        "manual_pass_annotations": "generic cytochrome P450 family protein",
        "manual_pass_taxonomic_header_scope": "Streptomyces",
        "primary_baseline_rows": len(baseline),
        "primary_combined_rows": len(primary),
        "close_refseq_sensitivity_rows": len(close_only),
        "new_exact_duplicates_of_primary": sum(row["decision"] == "exclude_exact_sequence_duplicate" for row in audit),
        "functional_scope": "TamI-focused provisional homologs",
        "orthology_BGC_limit": "Functional identity and BGC neighborhood are not established by the supplied files",
        "primary_definition": "curated legacy v2 plus manual RefSeq hits passing fixed numeric filters",
        "sensitivity_definition": "reference plus manual RefSeq hits passing fixed numeric filters; excludes unresolved legacy homologs",
    }
    (out / "intake_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    note = (
        "Manual review note\n"
        "All 100 supplied FASTA headers describe cytochrome P450 family proteins from Streptomyces.\n"
        "The supplied files contain no BGC neighborhood or experimentally validated function field.\n"
        "Records passing the manual's identity, query coverage, length, ambiguity, and duplicate rules are retained as provisional high-similarity homologs.\n"
        "They are not claimed to be confirmed TamI orthologs. The close-only set is analyzed separately to expose sensitivity to unresolved legacy membership.\n"
    )
    (out / "manual_review_note.txt").write_text(note, encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
