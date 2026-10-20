from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "TamI_round2_supplemental_search_v1"
SOURCE = RUN / "01_clusterednr_manual"
OUT = RUN / "05_clusterednr_tami_qc"
FASTA = SOURCE / "TamI_R2_representative_ClusteredNR_complete.txt"
HITS = SOURCE / "TamI_R2_representative_ClusteredNR_hits.csv"
BLAST = SOURCE / "clusterednr250_vs_TamI_reference.tsv"
EXISTING = ROOT / "TamI_round2_expanded_results_v1" / "02_intake" / "combined_curated_expanded_v2.faa"
UNIREF = RUN / "04_uniref_tami_qc" / "uniref90_retained_after_TamI_QC.faa"
EXCLUDED_BASELINE = ROOT / "TamI_round2_results" / "02_intake" / "legacy_curation_v2" / "excluded_baseline_full_length.faa"
PRE_MOTIF_INPUT = ROOT / "TamI_round2_expanded_results_v1" / "02_intake" / "combined_curated_expanded_v1.faa"
KNOWN_MOTIF_REJECT_IDS = {"WP_453160342.1", "MDT7786182.1"}
REFERENCE_LENGTH = 413
AA = set("ACDEFGHIKLMNPQRSTVWY")
AUX_FIELDS = [
    "query_id", "subject_accession", "aux_pident", "alignment_length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "aux_last_percent_field",
]
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
    identifiers = [row["id"] for row in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate FASTA IDs in {path}")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_aux_hits(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line_number, values in enumerate(csv.reader(handle), 1):
            if not values:
                continue
            if len(values) != len(AUX_FIELDS):
                raise ValueError(f"Expected 13 fields at {path}:{line_number}, got {len(values)}")
            row: dict[str, object] = dict(zip(AUX_FIELDS, values))
            for key in ("aux_pident", "evalue", "bitscore", "aux_last_percent_field"):
                row[key] = float(row[key])
            for key in ("alignment_length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"):
                row[key] = int(row[key])
            row["source_line"] = line_number
            rows.append(row)
    return rows


def read_best_blast(path: Path) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        values = line.split("\t")
        if len(values) != len(BLAST_FIELDS):
            raise ValueError(f"Expected 15 BLAST fields at line {line_number}, got {len(values)}")
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
    candidates = read_fasta(FASTA)
    hits = read_aux_hits(HITS)
    blast = read_best_blast(BLAST)
    candidate_ids = {row["id"] for row in candidates}
    hit_ids = {str(row["subject_accession"]) for row in hits}
    if candidate_ids != hit_ids:
        raise ValueError(f"FASTA/hit accession mismatch: fasta_only={candidate_ids-hit_ids}, hit_only={hit_ids-candidate_ids}")
    if set(blast) != candidate_ids:
        raise ValueError(f"TamI-reference BLAST did not cover all candidates: {candidate_ids-set(blast)}")
    queries = list(dict.fromkeys(str(row["query_id"]) for row in hits))
    expected_queries = {f"TamI_branch_{index:02d}|" for index in range(1, 7)}
    if {query.split("TamI_REF")[0] if "TamI_REF" in query else query.rsplit("|", 1)[0] + "|" for query in queries} != expected_queries:
        labels = {query.split("|", 1)[0] for query in queries}
        if labels != {f"TamI_branch_{index:02d}" for index in range(1, 7)}:
            raise ValueError(f"Unexpected auxiliary query set: {queries}")

    support: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in hits:
        support[str(row["subject_accession"])].append(row)
    existing = read_fasta(EXISTING)
    uniref = read_fasta(UNIREF)
    excluded_baseline = read_fasta(EXCLUDED_BASELINE)
    pre_motif = {row["id"]: row["sequence"] for row in read_fasta(PRE_MOTIF_INPUT)}
    known_motif_reject_sequences = {pre_motif[identifier]: identifier for identifier in KNOWN_MOTIF_REJECT_IDS}
    excluded_baseline_sequences = {row["sequence"]: row["id"] for row in excluded_baseline}
    existing_ids = {row["id"] for row in existing}
    seen_sequences = {row["sequence"]: row["id"] for row in existing + uniref}
    seen_ids = existing_ids | {row["id"] for row in uniref}
    sequence_counts = Counter(row["sequence"] for row in candidates)
    audit: list[dict[str, object]] = []
    retained: list[tuple[str, str]] = []

    for candidate in candidates:
        accession = candidate["id"]
        sequence = candidate["sequence"]
        hit = blast[accession]
        annotation = candidate["header"][len(accession):].strip()
        aux_rows = support[accession]
        query_labels = list(dict.fromkeys(str(row["query_id"]).split("|", 1)[0] for row in aux_rows))
        failures: list[str] = []
        identity = float(hit["pident"])
        query_span = (int(hit["qend"]) - int(hit["qstart"]) + 1) / REFERENCE_LENGTH
        qcovs = float(hit["qcovs"])
        evalue = float(hit["evalue"])
        if evalue > 1e-20:
            failures.append("TamI_reference_evalue_gt_fixed_1e-20")
        if not 25.0 <= identity <= 95.0:
            failures.append("identity_outside_fixed_25_95")
        if query_span < 0.70:
            failures.append("query_span_coverage_lt_fixed_0.70")
        if qcovs < 70.0:
            failures.append("reported_qcovs_lt_fixed_70_percent")
        length_ratio = len(sequence) / REFERENCE_LENGTH
        if not 0.75 <= length_ratio <= 1.25:
            failures.append("length_ratio_outside_fixed_0.75_1.25")
        nonstandard_fraction = sum(character not in AA for character in sequence) / len(sequence)
        if nonstandard_fraction > 0.01:
            failures.append("nonstandard_fraction_gt_fixed_0.01")
        matched_terms = [term for term in OTHER_FUNCTION_TERMS if term in annotation.lower()]
        if matched_terms:
            failures.append("explicit_other_function_annotation")
        if sequence in excluded_baseline_sequences:
            failures.append("exact_match_to_known_explicit_other_function_baseline")
        if sequence in known_motif_reject_sequences:
            failures.append(f"known_post_MSA_motif_reject:{known_motif_reject_sequences[sequence]}")
        duplicate_of = seen_sequences.get(sequence, "")
        if accession in seen_ids:
            failures.append("accession_already_in_existing_or_uniref_input")
        if duplicate_of:
            failures.append(f"exact_full_length_duplicate_of:{duplicate_of}")
        decision = "exclude" if failures else "retain_provisional_homolog"
        if not failures:
            retained.append((accession, sequence))
            seen_ids.add(accession)
            seen_sequences[sequence] = accession

        audit.append({
            "accession": accession,
            "aux_query_support_count": len(query_labels),
            "aux_query_support": ";".join(query_labels),
            "best_aux_evalue": min(float(row["evalue"]) for row in aux_rows),
            "max_aux_identity_percent": max(float(row["aux_pident"]) for row in aux_rows),
            "aux_last_percent_field_note": "decimal percent from downloaded hit table; inferred as positives, not used for intake QC",
            "TamI_reference_identity_percent": identity,
            "TamI_reference_qcovs_percent": qcovs,
            "TamI_reference_query_span_coverage": query_span,
            "TamI_reference_evalue": evalue,
            "TamI_reference_bitscore": hit["bitscore"],
            "sequence_length": len(sequence),
            "length_ratio": length_ratio,
            "nonstandard_fraction": nonstandard_fraction,
            "exact_sequence_multiplicity_in_download": sequence_counts[sequence],
            "duplicate_of_existing_or_uniref": duplicate_of,
            "known_motif_reject_match": known_motif_reject_sequences.get(sequence, ""),
            "known_explicit_other_function_match": excluded_baseline_sequences.get(sequence, ""),
            "explicit_other_function_terms": ";".join(matched_terms),
            "annotation": annotation,
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "final_decision": decision,
            "exclusion_reasons": ";".join(failures),
        })

    merged = [(row["id"], row["sequence"]) for row in existing + uniref] + retained
    if len({identifier for identifier, _ in merged}) != len(merged):
        raise AssertionError("Duplicate IDs in merged input")
    if len({sequence for _, sequence in merged}) != len(merged):
        raise AssertionError("Duplicate full-length sequences in merged input")
    write_tsv(OUT / "clusterednr_representative_query_TamI_QC_audit.tsv", audit)
    write_fasta(OUT / "clusterednr_new_retained_after_TamI_QC.faa", retained)
    write_fasta(OUT / "combined_round2_uniref_clusterednr_provisional.faa", merged)
    reasons = Counter(reason for row in audit for reason in str(row["exclusion_reasons"]).split(";") if reason)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "downloaded_hit_rows": len(hits),
        "representative_queries": len(queries),
        "unique_candidate_accessions": len(candidates),
        "candidate_fasta_byte_identical_to_prior_clusteredNR_batch01": digest(FASTA) == digest(ROOT / "TamI_round2_expanded_input_v1" / "TamI_R2_clusteredNR_batch01_250.fasta"),
        "TamI_reference_alignments": len(blast),
        "new_after_existing_and_uniref_dedup": len(retained),
        "combined_provisional_sequences": len(merged),
        "decision_counts": Counter(str(row["final_decision"]) for row in audit),
        "exclusion_reason_counts_nonexclusive": reasons,
        "query_hit_counts": Counter(str(row["query_id"]).split("|", 1)[0] for row in hits),
        "fixed_qc_thresholds_changed": False,
        "input_sha256": {"hits": digest(HITS), "fasta": digest(FASTA), "TamI_reference_blast": digest(BLAST)},
        "output_sha256": {
            "audit": digest(OUT / "clusterednr_representative_query_TamI_QC_audit.tsv"),
            "retained_fasta": digest(OUT / "clusterednr_new_retained_after_TamI_QC.faa"),
            "combined_provisional_fasta": digest(OUT / "combined_round2_uniref_clusterednr_provisional.faa"),
        },
    }
    (OUT / "clusterednr_representative_query_QC_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=dict) + "\n", encoding="utf-8")
    (OUT / "curation_note.txt").write_text(
        "Representative-query ClusteredNR candidates were all re-aligned to TamI reference ADC79647.1 and screened with the frozen thresholds. Auxiliary-query metrics were retained only as discovery provenance. Previously rejected post-MSA motif failures were explicitly blocked from re-entry.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=dict), flush=True)


if __name__ == "__main__":
    main()
