from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "TamI_round2_expanded_results_v1"
IN_FASTA = RESULTS / "02_intake" / "combined_curated_expanded_v1.faa"
MOTIF_AUDIT = RESULTS / "05_qc" / "TamI_round2_expanded250_v1_motif_alignment_audit.tsv"
CANDIDATE_AUDIT = RESULTS / "02_intake" / "expanded_candidate_accession_audit.tsv"
OUT_FASTA = RESULTS / "02_intake" / "combined_curated_expanded_v2.faa"
OUT_AUDIT = RESULTS / "02_intake" / "post_msa_motif_review.tsv"
OUT_NOTE = RESULTS / "02_intake" / "curation_note_v2.txt"
EXPECTED_EXCLUSIONS = {"WP_453160342.1", "MDT7786182.1"}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                rows.append((header.split()[0], "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ValueError("Sequence before FASTA header")
            chunks.append("".join(line.split()))
    if header is not None:
        rows.append((header.split()[0], "".join(chunks).upper()))
    return rows


def write_fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, sequence in rows:
            handle.write(f">{sequence_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    input_rows = read_fasta(IN_FASTA)
    motif_rows = read_tsv(MOTIF_AUDIT)
    candidate_rows = read_tsv(CANDIDATE_AUDIT)
    candidate_by_id = {row["accession"]: row for row in candidate_rows}
    failed = {row["sequence_id"] for row in motif_rows if row["three_anchor_match"] != "True"}
    if failed != EXPECTED_EXCLUSIONS:
        raise AssertionError(f"Unexpected motif-review set: {sorted(failed)}")
    motif_by_id = {row["sequence_id"]: row for row in motif_rows}

    review_rows: list[dict[str, object]] = []
    for sequence_id in sorted(failed):
        candidate = candidate_by_id[sequence_id]
        motif = motif_by_id[sequence_id]
        if sequence_id == "WP_453160342.1":
            reason = "supplied_header_marks_partial; C359_missing_in_reference_projection; incomplete_C_terminal_heme_region"
        else:
            reason = "C359_and_entire_heme_window_missing_in_reference_projection; truncated_or_domain_mismatched_sequence"
        review_rows.append(
            {
                "sequence_id": sequence_id,
                "decision": "exclude_after_post_msa_motif_review",
                "reason": reason,
                "selected_source": candidate["selected_source"],
                "pident": candidate["pident"],
                "reported_qcovs_percent": candidate["reported_qcovs_percent"],
                "query_span_coverage": candidate["query_span_coverage"],
                "fasta_length": candidate["fasta_length"],
                "annotation": candidate["annotation"],
                "organism": candidate["organism"],
                "E287": motif["E287"],
                "R290": motif["R290"],
                "E305_manual_checkpoint": motif["E305_manual_checkpoint"],
                "C359": motif["C359"],
                "ExxR_window": motif["ExxR_window"],
                "heme_window": motif["heme_window"],
            }
        )

    curated = [(sequence_id, sequence) for sequence_id, sequence in input_rows if sequence_id not in failed]
    if len(input_rows) != 193 or len(curated) != 191:
        raise AssertionError(f"Expected input/curated N=193/191, found {len(input_rows)}/{len(curated)}")
    if curated[0][0] != "TamI_REF_ADC79647.1":
        raise AssertionError("Reference not first after post-MSA curation")
    write_fasta(OUT_FASTA, curated)

    with OUT_AUDIT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(review_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(review_rows)

    OUT_NOTE.write_text(
        "Expanded TamI round-2 curated-v2 note\n"
        "This version starts from combined_curated_expanded_v1.faa and preserves all fixed numeric thresholds.\n"
        "Post-MSA conserved-motif review excluded exactly two incomplete/domain-truncated records: WP_453160342.1 and MDT7786182.1.\n"
        "WP_453160342.1 is labeled partial in the supplied header and lacks C359 in the reference projection.\n"
        "MDT7786182.1 lacks C359 and the entire projected heme window.\n"
        "All other 189 provisional-input rows in the first alignment retained the three operational anchors E287/R290/C359.\n"
        "TamI orthology and BGC context remain unverified.\n",
        encoding="utf-8",
    )
    hashes = {
        "combined_curated_expanded_v2.faa": hashlib.sha256(OUT_FASTA.read_bytes()).hexdigest(),
        "post_msa_motif_review.tsv": hashlib.sha256(OUT_AUDIT.read_bytes()).hexdigest(),
        "curation_note_v2.txt": hashlib.sha256(OUT_NOTE.read_bytes()).hexdigest(),
    }
    print(f"Input N={len(input_rows)}; curated-v2 N={len(curated)}; excluded={sorted(failed)}")
    for name, value in hashes.items():
        print(f"{value}  {name}")


if __name__ == "__main__":
    main()
