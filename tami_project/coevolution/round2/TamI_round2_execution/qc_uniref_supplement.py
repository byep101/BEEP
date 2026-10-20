from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "TamI_round2_supplemental_search_v1"
INPUT = RUN / "03_uniref_candidates"
OUT = RUN / "04_uniref_tami_qc"
BLAST = INPUT / "uniref90_vs_TamI_reference.tsv"
RAW_FASTA = INPUT / "uniref90_representatives_raw.faa"
SEARCH_AUDIT = INPUT / "uniref90_auxiliary_search_hits.tsv"
EXISTING = ROOT / "TamI_round2_expanded_results_v1" / "02_intake" / "combined_curated_expanded_v2.faa"
EXCLUDED = ROOT / "TamI_round2_results" / "02_intake" / "legacy_curation_v2" / "excluded_baseline_full_length.faa"
REFERENCE_LENGTH = 413
AA = set("ACDEFGHIKLMNPQRSTVWY")
BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore", "qcovs", "qlen", "slen",
]
OTHER_FUNCTION_TERMS = [
    "pikc", "vitamin d", "deoxyerythronolide", "pkss", "mevastatin", "nocardicin",
    "mycinamicin", "biotin biosynthesis", "pentalenic", "neocarzinostatin", "cyclomarin",
    "biflaviolin", "tylh1", "steroid", "vancomycin", "pentalenolactone", "erythromycin",
    "epothilone", "methyl-branched lipid", "dihydromenaquinone", "linalool",
    "methoxybenzoate", "pulcherriminic", "azinomycin", "doxa", "gfsf", "fd-891",
    "fd-892", "sulfonylurea", "terpineol", "thiocarbamate", "pravastatin",
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append({"id": header.split()[0], "header": header, "sequence": "".join(chunks).upper()})
            header = line[1:].strip()
            chunks = []
        elif line.strip():
            if header is None:
                raise ValueError(f"Sequence before header in {path}")
            chunks.append("".join(line.split()).replace("-", "").replace(".", ""))
    if header is not None:
        records.append({"id": header.split()[0], "header": header, "sequence": "".join(chunks).upper()})
    ids = [row["id"] for row in records]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate FASTA IDs in {path}")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for accession, sequence in records:
            handle.write(f">{accession}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_best_blast(path: Path) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        values = line.split("\t")
        if len(values) != len(BLAST_FIELDS):
            raise ValueError(f"Expected {len(BLAST_FIELDS)} BLAST fields at line {line_number}, got {len(values)}")
        row: dict[str, object] = dict(zip(BLAST_FIELDS, values))
        for key in ("pident", "evalue", "bitscore", "qcovs"):
            row[key] = float(row[key])
        for key in ("length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qlen", "slen"):
            row[key] = int(row[key])
        if int(row["qlen"]) != REFERENCE_LENGTH:
            raise ValueError(f"Unexpected query length for {row['sseqid']}: {row['qlen']}")
        accession = str(row["sseqid"])
        if accession not in best or float(row["bitscore"]) > float(best[accession]["bitscore"]):
            best[accession] = row
    return best


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    candidates = read_fasta(RAW_FASTA)
    existing = read_fasta(EXISTING)
    excluded = read_fasta(EXCLUDED)
    blast = read_best_blast(BLAST)
    if not set(blast).issubset({row["id"] for row in candidates}):
        raise ValueError("BLAST output contains an unknown candidate ID")

    descriptions: dict[str, list[str]] = {}
    query_support: dict[str, list[str]] = {}
    with SEARCH_AUDIT.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["uniref_accession"]
            descriptions.setdefault(accession, []).append(row["aux_query_description"])
            query_support.setdefault(accession, []).append(row["query_label"])

    existing_ids = {row["id"] for row in existing}
    existing_sequences = {row["sequence"]: row["id"] for row in existing}
    excluded_sequences = {row["sequence"]: row["id"] for row in excluded}
    candidate_multiplicity = Counter(row["sequence"] for row in candidates)
    seen_retained = dict(existing_sequences)
    audit: list[dict[str, object]] = []
    retained: list[tuple[str, str]] = []

    for candidate in candidates:
        accession = candidate["id"]
        sequence = candidate["sequence"]
        hit = blast.get(accession)
        annotation = " | ".join(dict.fromkeys(descriptions.get(accession, [candidate["header"]])))
        matched_terms = [term for term in OTHER_FUNCTION_TERMS if term in annotation.lower()]
        failures: list[str] = []
        if hit is None:
            identity = qcovs = query_span = evalue = bitscore = ""
            failures.append("no_TamI_reference_alignment")
        else:
            identity = float(hit["pident"])
            qcovs = float(hit["qcovs"])
            query_span = (int(hit["qend"]) - int(hit["qstart"]) + 1) / REFERENCE_LENGTH
            evalue = float(hit["evalue"])
            bitscore = float(hit["bitscore"])
            if not 25.0 <= identity <= 95.0:
                failures.append("identity_outside_fixed_25_95")
            if query_span < 0.70:
                failures.append("query_span_coverage_lt_fixed_0.70")
            if qcovs < 70.0:
                failures.append("reported_qcovs_lt_fixed_70_percent")
            if evalue > 1e-20:
                failures.append("TamI_reference_evalue_gt_fixed_1e-20")
        length_ratio = len(sequence) / REFERENCE_LENGTH
        if not 0.75 <= length_ratio <= 1.25:
            failures.append("length_ratio_outside_fixed_0.75_1.25")
        nonstandard_symbols = sorted(set(sequence) - AA)
        nonstandard_fraction = sum(character not in AA for character in sequence) / len(sequence)
        if nonstandard_fraction > 0.01:
            failures.append("nonstandard_fraction_gt_fixed_0.01")
        if matched_terms:
            failures.append("explicit_other_function_annotation")
        if sequence in excluded_sequences:
            failures.append("exact_match_to_known_explicit_other_function_baseline")

        duplicate_of = seen_retained.get(sequence, "")
        if accession in existing_ids:
            failures.append("accession_already_in_round2_input")
        if duplicate_of:
            failures.append(f"exact_full_length_duplicate_of:{duplicate_of}")
        decision = "exclude" if failures else "retain_provisional_homolog"
        if not failures:
            retained.append((accession, sequence))
            seen_retained[sequence] = accession

        audit.append({
            "uniref_accession": accession,
            "query_support_count": len(set(query_support.get(accession, []))),
            "query_support": ";".join(dict.fromkeys(query_support.get(accession, []))),
            "TamI_reference_identity_percent": identity,
            "TamI_reference_qcovs_percent": qcovs,
            "TamI_reference_query_span_coverage": query_span,
            "TamI_reference_evalue": evalue,
            "TamI_reference_bitscore": bitscore,
            "sequence_length": len(sequence),
            "length_ratio": length_ratio,
            "nonstandard_symbols": "".join(nonstandard_symbols),
            "nonstandard_fraction": nonstandard_fraction,
            "exact_sequence_multiplicity_in_uniref_candidates": candidate_multiplicity[sequence],
            "duplicate_of_existing_or_prior_retained": duplicate_of,
            "known_excluded_baseline_match": excluded_sequences.get(sequence, ""),
            "explicit_other_function_terms": ";".join(matched_terms),
            "annotation": annotation,
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "final_decision": decision,
            "exclusion_reasons": ";".join(failures),
            "orthology_status": "not_established",
            "BGC_status": "not_checked",
        })

    merged = [(row["id"], row["sequence"]) for row in existing] + retained
    if len({identifier for identifier, _ in merged}) != len(merged):
        raise AssertionError("Duplicate IDs in merged input")
    if len({sequence for _, sequence in merged}) != len(merged):
        raise AssertionError("Duplicate full-length sequences in merged input")
    if len(merged) > 1000:
        raise AssertionError("Merged input exceeds fixed 1000-sequence limit; do not truncate")

    write_tsv(OUT / "uniref90_TamI_reference_QC_audit.tsv", audit)
    write_fasta(OUT / "uniref90_retained_after_TamI_QC.faa", retained)
    write_fasta(OUT / "combined_round2_plus_uniref_provisional.faa", merged)
    decisions = Counter(str(row["final_decision"]) for row in audit)
    reasons = Counter(
        reason
        for row in audit
        for reason in str(row["exclusion_reasons"]).split(";")
        if reason
    )
    retained_rows = [row for row in audit if row["final_decision"] == "retain_provisional_homolog"]
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_clusters": len(candidates),
        "TamI_reference_alignments": len(blast),
        "existing_round2_full_length_sequences": len(existing),
        "retained_new_unique_sequences": len(retained),
        "provisional_merged_sequences": len(merged),
        "decision_counts": decisions,
        "exclusion_reason_counts_nonexclusive": reasons,
        "retained_identity_percent_range": [
            min(float(row["TamI_reference_identity_percent"]) for row in retained_rows),
            max(float(row["TamI_reference_identity_percent"]) for row in retained_rows),
        ] if retained_rows else [],
        "retained_qcovs_percent_range": [
            min(float(row["TamI_reference_qcovs_percent"]) for row in retained_rows),
            max(float(row["TamI_reference_qcovs_percent"]) for row in retained_rows),
        ] if retained_rows else [],
        "fixed_qc": {
            "TamI_reference_evalue_max": 1e-20,
            "identity_percent_inclusive": [25, 95],
            "query_span_coverage_min": 0.70,
            "qcovs_percent_min": 70,
            "length_ratio_inclusive": [0.75, 1.25],
            "nonstandard_fraction_max": 0.01,
            "exact_full_length_deduplication": True,
            "thresholds_changed": False,
        },
        "input_sha256": {
            "uniref_candidates": digest(RAW_FASTA),
            "TamI_reference_blast": digest(BLAST),
            "existing_round2_input": digest(EXISTING),
            "known_excluded_baseline": digest(EXCLUDED),
        },
        "output_sha256": {
            "audit": digest(OUT / "uniref90_TamI_reference_QC_audit.tsv"),
            "retained_fasta": digest(OUT / "uniref90_retained_after_TamI_QC.faa"),
            "provisional_merged_fasta": digest(OUT / "combined_round2_plus_uniref_provisional.faa"),
        },
        "downstream_status": "Provisional intake only. Do not run final MSA/weights until the representative-query ClusteredNR branch is supplied and processed.",
    }
    (OUT / "uniref90_TamI_reference_QC_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=dict) + "\n", encoding="utf-8")
    (OUT / "curation_note.txt").write_text(
        "UniRef90 representative-query supplemental intake. Every candidate was re-aligned to TamI reference ADC79647.1. Auxiliary-query identity was not used for acceptance. Frozen search/QC limits and exact full-length deduplication were applied. Retained entries remain provisional homologs; TamI orthology and BGC context are not established. Final MSA/weighting is deferred until the matching ClusteredNR supplemental branch is available.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=dict), flush=True)


if __name__ == "__main__":
    main()
