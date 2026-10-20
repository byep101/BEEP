"""Validate the manually exported RefSeq FASTA and BLAST hit table."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


AA = set("ACDEFGHIKLMNPQRSTVWY")
FIELDS = [
    "qseqid", "saccver", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    header = None
    chunks: list[str] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header.split()[0], header, "".join(chunks).upper()))
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Empty FASTA header at line {number}")
            chunks = []
        else:
            if header is None:
                raise ValueError(f"Sequence before first header at line {number}")
            chunks.append("".join(line.split()))
    if header is not None:
        records.append((header.split()[0], header, "".join(chunks).upper()))
    ids = [row[0] for row in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate FASTA IDs: " + ", ".join(k for k, v in Counter(ids).items() if v > 1))
    return records


def read_hits(path: Path) -> list[dict[str, object]]:
    hits = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle), 1):
            if not row or not any(value.strip() for value in row):
                continue
            if len(row) != 13:
                raise ValueError(f"Expected 13 CSV fields at line {line_number}, got {len(row)}")
            hit: dict[str, object] = dict(zip(FIELDS, row))
            for key in ("pident", "evalue", "bitscore", "qcovs"):
                hit[key] = float(hit[key])
            for key in ("length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"):
                hit[key] = int(hit[key])
            hit["source_line"] = line_number
            hits.append(hit)
    return hits


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    workspace = Path("C:/Users/31441/Desktop/work/BEEP")
    source = workspace / "TamI_round2_manual_input"
    fasta_path = source / "TamI_R2_RefSeq_raw.fasta"
    hits_path = source / "TamI_R2_RefSeq_hits.csv"
    output = source / "validation"
    output.mkdir(exist_ok=True)

    fasta = read_fasta(fasta_path)
    hits = read_hits(hits_path)
    by_id = {identifier: (header, sequence) for identifier, header, sequence in fasta}
    hit_counts = Counter(str(hit["saccver"]) for hit in hits)
    if any(count > 1 for count in hit_counts.values()):
        raise ValueError("Hit table has multiple rows for an accession; best-HSP selection must be explicit")
    fasta_ids = set(by_id)
    hit_ids = set(hit_counts)

    excluded_path = workspace / "TamI_round2_manual_input" / "known_excluded_baseline_full_length.faa"
    excluded_sequences: dict[str, str] = {}
    if excluded_path.exists():
        excluded_sequences = {sequence: identifier for identifier, _, sequence in read_fasta(excluded_path)}

    audit: list[dict[str, object]] = []
    sequence_counts = Counter(sequence for _, _, sequence in fasta)
    for hit in hits:
        accession = str(hit["saccver"])
        header, sequence = by_id[accession]
        query_span = (int(hit["qend"]) - int(hit["qstart"]) + 1) / 413
        nonstandard = sorted(set(sequence) - AA)
        numeric_failures = []
        if not 25 <= float(hit["pident"]) <= 95:
            numeric_failures.append("identity_outside_manual_25_95")
        if query_span < 0.7:
            numeric_failures.append("query_span_coverage_lt_0.7")
        if float(hit["qcovs"]) < 70:
            numeric_failures.append("reported_qcovs_lt_70")
        if not 0.75 <= len(sequence) / 413 <= 1.25:
            numeric_failures.append("length_ratio_outside_0.75_1.25")
        if len(nonstandard) / len(sequence) > 0.01:
            numeric_failures.append("nonstandard_fraction_gt_0.01")
        duplicate_of_excluded = excluded_sequences.get(sequence, "")
        audit.append({
            **hit,
            "fasta_length": len(sequence),
            "length_ratio": len(sequence) / 413,
            "query_span_coverage": query_span,
            "reported_qcovs_fraction": float(hit["qcovs"]) / 100,
            "nonstandard_symbols": "".join(nonstandard),
            "exact_sequence_multiplicity_in_download": sequence_counts[sequence],
            "exact_match_to_known_excluded_baseline": duplicate_of_excluded,
            "numeric_decision": "exclude" if numeric_failures or duplicate_of_excluded else "pass",
            "numeric_reasons": ";".join(numeric_failures + (["matches_known_excluded_baseline"] if duplicate_of_excluded else [])),
            "fasta_header": header,
        })

    identities = [float(row["pident"]) for row in hits]
    lengths = [len(sequence) for _, _, sequence in fasta]
    summary = {
        "fasta_sha256": sha256(fasta_path),
        "hit_table_sha256": sha256(hits_path),
        "fasta_records": len(fasta),
        "nonempty_hit_rows": len(hits),
        "unique_hit_accessions": len(hit_ids),
        "fasta_only_accessions": sorted(fasta_ids - hit_ids),
        "hit_only_accessions": sorted(hit_ids - fasta_ids),
        "all_accessions_match_one_to_one": fasta_ids == hit_ids and len(fasta) == len(hits),
        "query_ids": sorted(set(str(row["qseqid"]) for row in hits)),
        "query_coordinate_convention": "1-based inclusive",
        "hit_table_format": FIELDS,
        "identity_min": min(identities),
        "identity_max": max(identities),
        "identity_gt_95_count": sum(value > 95 for value in identities),
        "identity_25_to_95_count": sum(25 <= value <= 95 for value in identities),
        "sequence_length_min": min(lengths),
        "sequence_length_max": max(lengths),
        "exact_duplicate_sequence_groups": sum(count > 1 for count in sequence_counts.values()),
        "numeric_pass_count": sum(row["numeric_decision"] == "pass" for row in audit),
        "numeric_exclude_count": sum(row["numeric_decision"] == "exclude" for row in audit),
        "known_excluded_sequence_reference_loaded": bool(excluded_sequences),
    }
    if not summary["all_accessions_match_one_to_one"]:
        raise ValueError("FASTA and hit table accessions do not match one-to-one")
    write_tsv(output / "manual_hit_audit.tsv", audit)
    (output / "manual_input_validation.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
