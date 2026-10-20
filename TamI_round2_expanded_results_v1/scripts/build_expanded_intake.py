from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "TamI_round2_expanded_input_v1"
OUT = ROOT / "TamI_round2_expanded_results_v1" / "02_intake"
BASELINE = ROOT / "TamI_round2_results" / "02_intake" / "legacy_curation_v2" / "baseline_curated_unaligned.faa"
EXCLUDED_BASELINE = ROOT / "TamI_round2_results" / "02_intake" / "legacy_curation_v2" / "excluded_baseline_full_length.faa"
REFERENCE_ID = "TamI_REF_ADC79647.1"
REFERENCE_LENGTH = 413
AA = set("ACDEFGHIKLMNPQRSTVWY")
FIELDS = [
    "qseqid", "saccver", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs",
]
SOURCES = [
    {
        "name": "RefSeq",
        "hits": INPUT / "TamI_R2_RefSeq_hits_250.csv",
        "fasta": INPUT / "TamI_R2_RefSeq_raw_250.fasta",
    },
    {
        "name": "clusteredNR_batch01",
        "hits": INPUT / "TamI_R2_clusteredNR_batch01_hits_250.csv",
        "fasta": INPUT / "TamI_R2_clusteredNR_batch01_250.fasta",
    },
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                rows.append({"id": header.split()[0], "header": header, "sequence": "".join(chunks).upper()})
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Empty FASTA header at {path}:{line_number}")
            chunks = []
        else:
            if header is None:
                raise ValueError(f"Sequence before FASTA header at {path}:{line_number}")
            chunks.append("".join(line.split()))
    if header is not None:
        rows.append({"id": header.split()[0], "header": header, "sequence": "".join(chunks).upper()})
    ids = [row["id"] for row in rows]
    duplicated = [key for key, count in Counter(ids).items() if count > 1]
    if duplicated:
        raise ValueError(f"Duplicate FASTA IDs in {path}: {duplicated}")
    return rows


def read_hits(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line_number, values in enumerate(csv.reader(handle), 1):
            if not values or not any(value.strip() for value in values):
                continue
            if len(values) != len(FIELDS):
                raise ValueError(f"Expected 13 fields at {path}:{line_number}, got {len(values)}")
            row: dict[str, object] = dict(zip(FIELDS, values))
            for key in ("pident", "evalue", "bitscore", "qcovs"):
                row[key] = float(row[key])
            for key in ("length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"):
                row[key] = int(row[key])
            row["source_line"] = line_number
            rows.append(row)
    accessions = [str(row["saccver"]) for row in rows]
    duplicated = [key for key, count in Counter(accessions).items() if count > 1]
    if duplicated:
        raise ValueError(f"Multiple hit rows for an accession in {path}: {duplicated}")
    return rows


def write_fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, sequence in rows:
            handle.write(f">{sequence_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def annotation_parts(header: str, accession: str) -> tuple[str, str]:
    rest = header[len(accession):].strip()
    match = re.search(r"\[([^\[\]]+)\]\s*$", rest)
    if match:
        return rest[:match.start()].strip(), match.group(1)
    return rest, ""


def numeric_failures(hit: dict[str, object], sequence: str, known_excluded: dict[str, str]) -> list[str]:
    failures: list[str] = []
    identity = float(hit["pident"])
    query_span = (int(hit["qend"]) - int(hit["qstart"]) + 1) / REFERENCE_LENGTH
    qcovs = float(hit["qcovs"]) / 100.0
    length_ratio = len(sequence) / REFERENCE_LENGTH
    nonstandard = set(sequence) - AA
    if not 25.0 <= identity <= 95.0:
        failures.append("identity_outside_fixed_25_95")
    if query_span < 0.70:
        failures.append("query_span_coverage_lt_fixed_0.70")
    if qcovs < 0.70:
        failures.append("reported_qcovs_lt_fixed_0.70")
    if not 0.75 <= length_ratio <= 1.25:
        failures.append("length_ratio_outside_fixed_0.75_1.25")
    if nonstandard and len(nonstandard) / len(sequence) > 0.01:
        failures.append("nonstandard_fraction_gt_fixed_0.01")
    if sequence in known_excluded:
        failures.append("exact_match_to_known_explicit_other_function_baseline")
    return failures


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=False)
    baseline_records = read_fasta(BASELINE)
    if len(baseline_records) != 24:
        raise AssertionError(f"Expected 24 curated legacy records, found {len(baseline_records)}")
    if baseline_records[0]["id"] != REFERENCE_ID or len(baseline_records[0]["sequence"]) != REFERENCE_LENGTH:
        raise AssertionError("Curated reference is missing, not first, or not 413 aa")
    excluded_records = read_fasta(EXCLUDED_BASELINE)
    known_excluded = {row["sequence"]: row["id"] for row in excluded_records}

    source_audit: list[dict[str, object]] = []
    source_records_by_accession: dict[str, list[dict[str, object]]] = defaultdict(list)
    source_counts: dict[str, dict[str, int]] = {}
    source_hashes: dict[str, dict[str, str]] = {}

    for source_index, source in enumerate(SOURCES):
        source_name = str(source["name"])
        fasta_path = Path(source["fasta"])
        hits_path = Path(source["hits"])
        fasta_rows = read_fasta(fasta_path)
        hit_rows = read_hits(hits_path)
        if len(fasta_rows) != 250 or len(hit_rows) != 250:
            raise AssertionError(f"Expected 250 FASTA/hit rows for {source_name}, found {len(fasta_rows)}/{len(hit_rows)}")
        fasta_by_id = {row["id"]: row for row in fasta_rows}
        hit_ids = {str(row["saccver"]) for row in hit_rows}
        if set(fasta_by_id) != hit_ids:
            raise AssertionError(
                f"FASTA/hit accessions differ for {source_name}: "
                f"fasta_only={sorted(set(fasta_by_id) - hit_ids)}, hit_only={sorted(hit_ids - set(fasta_by_id))}"
            )
        sequence_counts = Counter(row["sequence"] for row in fasta_rows)
        source_counts[source_name] = {
            "fasta_records": len(fasta_rows),
            "hit_rows": len(hit_rows),
            "exact_duplicate_sequence_groups_within_source": sum(count > 1 for count in sequence_counts.values()),
        }
        source_hashes[source_name] = {"fasta_sha256": digest(fasta_path), "hits_sha256": digest(hits_path)}

        for hit in hit_rows:
            accession = str(hit["saccver"])
            fasta_row = fasta_by_id[accession]
            sequence = str(fasta_row["sequence"])
            header = str(fasta_row["header"])
            annotation, organism = annotation_parts(header, accession)
            failures = numeric_failures(hit, sequence, known_excluded)
            query_span = (int(hit["qend"]) - int(hit["qstart"]) + 1) / REFERENCE_LENGTH
            nonstandard = "".join(sorted(set(sequence) - AA))
            record: dict[str, object] = {
                "source": source_name,
                "source_priority": source_index,
                "source_line": hit["source_line"],
                "accession": accession,
                "pident": hit["pident"],
                "alignment_length": hit["length"],
                "qstart": hit["qstart"],
                "qend": hit["qend"],
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "reported_qcovs_percent": hit["qcovs"],
                "query_span_coverage": query_span,
                "fasta_length": len(sequence),
                "length_ratio": len(sequence) / REFERENCE_LENGTH,
                "nonstandard_symbols": nonstandard,
                "exact_sequence_multiplicity_within_source": sequence_counts[sequence],
                "known_excluded_baseline_match": known_excluded.get(sequence, ""),
                "numeric_decision": "exclude" if failures else "pass",
                "numeric_reasons": ";".join(failures),
                "annotation": annotation,
                "organism": organism,
                "annotation_scope_review": "provisional_generic_or_uncharacterized; no explicit alternative function in supplied header",
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "sequence": sequence,
            }
            source_records_by_accession[accession].append(record)
            source_audit.append({key: value for key, value in record.items() if key != "sequence"})

    accessions_in_order: list[str] = []
    seen_accessions: set[str] = set()
    for record in source_audit:
        accession = str(record["accession"])
        if accession not in seen_accessions:
            seen_accessions.add(accession)
            accessions_in_order.append(accession)

    baseline_ids = {row["id"] for row in baseline_records}
    seen_sequences = {row["sequence"]: row["id"] for row in baseline_records}
    retained_new: list[tuple[str, str]] = []
    candidate_audit: list[dict[str, object]] = []
    overlap_metric_discordance = 0
    overlap_qc_discordance = 0

    for accession in accessions_in_order:
        records = sorted(
            source_records_by_accession[accession],
            key=lambda row: (int(row["source_priority"]), int(row["source_line"])),
        )
        sequences = {str(row["sequence"]) for row in records}
        if len(sequences) != 1:
            raise AssertionError(f"Same accession has inconsistent protein sequence across sources: {accession}")
        selected = records[0]
        metric_tuples = {
            (
                float(row["pident"]), float(row["reported_qcovs_percent"]),
                int(row["qstart"]), int(row["qend"]), int(row["alignment_length"]),
            )
            for row in records
        }
        if len(metric_tuples) > 1:
            overlap_metric_discordance += 1
        qc_decisions = {str(row["numeric_decision"]) for row in records}
        if len(qc_decisions) > 1:
            overlap_qc_discordance += 1
        sequence = str(selected["sequence"])
        duplicate_of = seen_sequences.get(sequence, "")
        if accession in baseline_ids:
            final_decision = "exclude_accession_duplicate_of_curated_legacy"
            final_reason = "accession_already_present_in_curated_legacy"
        elif selected["numeric_decision"] != "pass":
            final_decision = "exclude_fixed_numeric_qc"
            final_reason = str(selected["numeric_reasons"])
        elif duplicate_of:
            final_decision = "exclude_exact_full_length_sequence_duplicate"
            final_reason = f"exact_full_length_duplicate_of:{duplicate_of}"
        else:
            final_decision = "retain_provisional_homolog"
            final_reason = ""
            retained_new.append((accession, sequence))
            seen_sequences[sequence] = accession

        candidate_audit.append(
            {
                "accession": accession,
                "source_support": ";".join(str(row["source"]) for row in records),
                "selected_source": selected["source"],
                "selected_source_line": selected["source_line"],
                "same_accession_sequence_consistent_across_sources": True,
                "source_metric_record_count": len(metric_tuples),
                "source_numeric_decisions": ";".join(str(row["numeric_decision"]) for row in records),
                "source_qc_decision_discordant": len(qc_decisions) > 1,
                "pident": selected["pident"],
                "reported_qcovs_percent": selected["reported_qcovs_percent"],
                "query_span_coverage": selected["query_span_coverage"],
                "fasta_length": selected["fasta_length"],
                "length_ratio": selected["length_ratio"],
                "numeric_decision": selected["numeric_decision"],
                "numeric_reasons": selected["numeric_reasons"],
                "annotation": selected["annotation"],
                "organism": selected["organism"],
                "annotation_scope_review": selected["annotation_scope_review"],
                "sequence_sha256": selected["sequence_sha256"],
                "final_intake_decision": final_decision,
                "final_intake_reason": final_reason,
                "duplicate_of_retained_sequence": duplicate_of,
            }
        )

    combined = [(str(row["id"]), str(row["sequence"])) for row in baseline_records] + retained_new
    if len(combined) > 1000:
        raise AssertionError(f"Combined N={len(combined)} exceeds fixed prepare limit 1000; do not truncate")
    if combined[0][0] != REFERENCE_ID:
        raise AssertionError("Reference is not first after merging")
    combined_ids = [sequence_id for sequence_id, _ in combined]
    combined_sequences = [sequence for _, sequence in combined]
    if len(combined_ids) != len(set(combined_ids)):
        raise AssertionError("Combined FASTA has duplicate accessions")
    if len(combined_sequences) != len(set(combined_sequences)):
        raise AssertionError("Combined FASTA has duplicate full-length sequences")

    write_tsv(OUT / "expanded_source_record_audit.tsv", source_audit)
    write_tsv(OUT / "expanded_candidate_accession_audit.tsv", candidate_audit)
    write_fasta(OUT / "combined_curated_expanded_v1.faa", combined)
    write_fasta(OUT / "new_retained_expanded_v1.faa", retained_new)
    (OUT / "curation_note.txt").write_text(
        "Expanded TamI round-2 intake curation note\n"
        "Scope choice carried forward from the user: retain the curated TamI-focused legacy set (reference plus 23 unresolved/generic homologs), excluding 45 legacy proteins with explicit alternative functions.\n"
        "New candidates were screened with the unchanged manual thresholds: identity 25-95% inclusive, query span coverage >=0.70, reported qcovs >=70%, length ratio 0.75-1.25, nonstandard amino acids <=1%, and exact full-length sequence deduplication.\n"
        "RefSeq takes precedence for accessions present in both downloads. Same-accession sequences were required to be identical across sources.\n"
        "Supplied new headers are generic cytochrome P450 or hypothetical-protein descriptions; no explicit alternative physiological function is asserted. TamI orthology and BGC context remain unverified.\n",
        encoding="utf-8",
    )

    retained_candidates = [row for row in candidate_audit if row["final_intake_decision"] == "retain_provisional_homolog"]
    summary = {
        "version": "expanded_v1",
        "fixed_qc": {
            "identity_percent_inclusive": [25, 95],
            "query_span_coverage_min": 0.70,
            "reported_qcovs_percent_min": 70,
            "length_ratio_inclusive": [0.75, 1.25],
            "nonstandard_fraction_max": 0.01,
            "exact_full_length_deduplication": True,
            "thresholds_changed": False,
        },
        "sources": source_counts,
        "source_hashes": source_hashes,
        "raw_source_rows_total": len(source_audit),
        "unique_accessions_after_cross_source_merge": len(candidate_audit),
        "accessions_supported_by_both_sources": sum(";" in str(row["source_support"]) for row in candidate_audit),
        "same_accession_metric_discordance_count": overlap_metric_discordance,
        "same_accession_numeric_qc_discordance_count": overlap_qc_discordance,
        "selected_numeric_pass_accessions": sum(row["numeric_decision"] == "pass" for row in candidate_audit),
        "selected_numeric_excluded_accessions": sum(row["numeric_decision"] == "exclude" for row in candidate_audit),
        "exact_full_length_duplicates_removed_after_numeric_qc": sum(row["final_intake_decision"] == "exclude_exact_full_length_sequence_duplicate" for row in candidate_audit),
        "retained_new_unique_full_length_sequences": len(retained_new),
        "curated_legacy_sequences": len(baseline_records),
        "combined_unaligned_sequences": len(combined),
        "combined_below_prepare_limit_1000": len(combined) <= 1000,
        "new_retained_identity_min": min(float(row["pident"]) for row in retained_candidates),
        "new_retained_identity_max": max(float(row["pident"]) for row in retained_candidates),
        "new_retained_reported_qcovs_min": min(float(row["reported_qcovs_percent"]) for row in retained_candidates),
        "new_retained_reported_qcovs_max": max(float(row["reported_qcovs_percent"]) for row in retained_candidates),
        "annotation_review": "generic cytochrome P450/cytochrome/hypothetical headers; no explicit alternative function in supplied headers",
        "functional_limit": "TamI orthology and BGC context are not established by the supplied files",
    }
    (OUT / "expanded_intake_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUT / "INTAKE_SHA256SUMS.txt").write_text(
        "".join(
            f"{digest(path)}  {path.name}\n"
            for path in [
                OUT / "combined_curated_expanded_v1.faa",
                OUT / "new_retained_expanded_v1.faa",
                OUT / "expanded_source_record_audit.tsv",
                OUT / "expanded_candidate_accession_audit.tsv",
                OUT / "expanded_intake_summary.json",
                OUT / "curation_note.txt",
            ]
        ),
        encoding="ascii",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
