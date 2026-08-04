#!/usr/bin/env python3
"""Reproducible helpers for the TamI/TxtC pilot MSA and tree workflow.

The implementation intentionally uses only the Python standard library so the
QC and metadata steps remain runnable on locked-down workstations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path


AA_ALLOWED = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-?")
GAP_CHARS = set("Xx-?*")


def read_fasta(path: str | Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with Path(path).open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append({"header": header, "id": header.split()[0], "seq": "".join(chunks).upper()})
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence found before FASTA header in {path}")
                chunks.append(re.sub(r"\s+", "", line))
    if header is not None:
        records.append({"header": header, "id": header.split()[0], "seq": "".join(chunks).upper()})
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def write_fasta(records: list[dict[str, str]], path: str | Path, width: int = 80) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(f">{record['header']}\n")
            seq = record["seq"]
            for start in range(0, len(seq), width):
                handle.write(seq[start : start + width] + "\n")


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def global_align(a: str, b: str) -> tuple[str, str]:
    match, mismatch, gap = 2, -1, -2
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0], trace[i][0] = i * gap, 1
    for j in range(1, m + 1):
        score[0][j], trace[0][j] = j * gap, 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            values = (score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch),
                      score[i - 1][j] + gap, score[i][j - 1] + gap)
            best = max(range(3), key=values.__getitem__)
            score[i][j], trace[i][j] = values[best], best
    out_a: list[str] = []
    out_b: list[str] = []
    i, j = n, m
    while i or j:
        direction = trace[i][j]
        if direction == 0:
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i -= 1; j -= 1
        elif direction == 1:
            out_a.append(a[i - 1]); out_b.append("-"); i -= 1
        else:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1
    return "".join(reversed(out_a)), "".join(reversed(out_b))


def cmd_compare(args: argparse.Namespace) -> None:
    a = read_fasta(args.fasta_a)[0]
    b = read_fasta(args.fasta_b)[0]
    aa, bb = global_align(a["seq"], b["seq"])
    rows: list[dict[str, object]] = []
    pos_a = pos_b = matches = 0
    for col, (ra, rb) in enumerate(zip(aa, bb), 1):
        if ra != "-": pos_a += 1
        if rb != "-": pos_b += 1
        if ra == rb and ra != "-": matches += 1
        if ra != rb:
            rows.append({"alignment_column": col, "position_a": pos_a if ra != "-" else "",
                         "residue_a": ra, "position_b": pos_b if rb != "-" else "", "residue_b": rb})
    summary = {
        "fasta_a": str(args.fasta_a), "id_a": a["id"], "length_a": len(a["seq"]),
        "fasta_b": str(args.fasta_b), "id_b": b["id"], "length_b": len(b["seq"]),
        "alignment_length": len(aa), "identical_positions": matches,
        "identity_fraction": round(matches / len(aa), 6), "difference_count": len(rows),
        "indel_columns": sum(1 for x, y in zip(aa, bb) if "-" in (x, y)),
        "coordinate_systems_share_positions": len(aa) == len(a["seq"]) == len(b["seq"]),
    }
    Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_json).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if args.output_tsv:
        with Path(args.output_tsv).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["alignment_column", "position_a", "residue_a", "position_b", "residue_b"], delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
    print(json.dumps(summary, ensure_ascii=False))


BLAST_FIELDS = ["qseqid", "saccver", "pident", "length", "qcovs", "evalue", "bitscore", "staxids", "sscinames", "stitle"]


def cmd_filter_blast(args: argparse.Namespace) -> None:
    best: dict[str, dict[str, str]] = {}
    with Path(args.input).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, fieldnames=BLAST_FIELDS, delimiter="\t")
        for row in reader:
            try:
                passed = (args.min_identity <= float(row["pident"]) <= args.max_identity
                          and float(row["qcovs"]) >= args.min_qcov
                          and float(row["evalue"]) <= args.max_evalue)
            except (TypeError, ValueError):
                continue
            if not passed:
                continue
            acc = row["saccver"].strip()
            if acc not in best or float(row["bitscore"]) > float(best[acc]["bitscore"]):
                best[acc] = row
    rows = sorted(best.values(), key=lambda row: float(row["bitscore"]), reverse=True)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.output_tsv).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BLAST_FIELDS, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    Path(args.output_accessions).write_text("".join(f"{r['saccver']}\n" for r in rows), encoding="utf-8")
    print(json.dumps({"filtered_hits": len(rows), "output": args.output_tsv}))


def cmd_fetch(args: argparse.Namespace) -> None:
    accessions = [line.strip() for line in Path(args.accessions).read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    target = Path(args.output); target.parent.mkdir(parents=True, exist_ok=True)
    chunks: list[str] = []
    for start in range(0, len(accessions), args.batch_size):
        batch = accessions[start : start + args.batch_size]
        params = {"db": "protein", "id": ",".join(batch), "rettype": "fasta", "retmode": "text", "tool": "BEEP_P450_phylo"}
        if args.email:
            params["email"] = args.email
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
        last_error: Exception | None = None
        for attempt in range(1, 6):
            try:
                with urllib.request.urlopen(url, timeout=120) as response:
                    text = response.read().decode("utf-8")
                if not text.lstrip().startswith(">"):
                    raise RuntimeError(f"Unexpected NCBI response: {text[:200]}")
                chunks.append(text.rstrip() + "\n")
                last_error = None
                break
            except Exception as exc:  # network retry is deliberate
                last_error = exc
                time.sleep(min(30, 2 ** attempt))
        if last_error:
            raise last_error
        time.sleep(0.34 if args.email else 0.5)
    target.write_text("".join(chunks), encoding="utf-8", newline="\n")
    downloaded = read_fasta(target)
    print(json.dumps({"requested": len(accessions), "downloaded": len(downloaded), "sha256": sha256(target)}))


def cmd_filter_fasta(args: argparse.Namespace) -> None:
    reference = read_fasta(args.reference)[0]
    ref_len = len(reference["seq"])
    output = [{"header": f"{args.reference_id} source={reference['id']}", "id": args.reference_id, "seq": reference["seq"]}]
    seen_sequences = {reference["seq"]}
    metadata: list[dict[str, object]] = [{"sequence_id": args.reference_id, "source_id": reference["id"], "length": ref_len,
                                          "length_ratio": 1.0, "x_fraction": 0.0, "status": "reference"}]
    excluded: list[dict[str, object]] = []
    for record in read_fasta(args.hits):
        seq = record["seq"].replace("-", "")
        ratio = len(seq) / ref_len
        x_fraction = seq.count("X") / len(seq) if seq else 1.0
        invalid = sorted(set(seq) - AA_ALLOWED)
        reason = ""
        if not args.min_length_ratio <= ratio <= args.max_length_ratio:
            reason = "length_ratio"
        elif x_fraction > args.max_x_fraction:
            reason = "x_fraction"
        elif invalid:
            reason = "invalid_characters:" + "".join(invalid)
        elif seq in seen_sequences:
            reason = "exact_sequence_duplicate"
        if reason:
            excluded.append({"sequence_id": record["id"], "length": len(seq), "length_ratio": round(ratio, 6),
                             "x_fraction": round(x_fraction, 6), "reason": reason})
            continue
        seen_sequences.add(seq)
        safe_id = re.sub(r"[^A-Za-z0-9_.|:-]", "_", record["id"])
        output.append({"header": safe_id, "id": safe_id, "seq": seq})
        metadata.append({"sequence_id": safe_id, "source_id": record["id"], "length": len(seq),
                         "length_ratio": round(ratio, 6), "x_fraction": round(x_fraction, 6), "status": "included"})
    write_fasta(output, args.output)
    for path, rows, fields in ((args.metadata, metadata, ["sequence_id", "source_id", "length", "length_ratio", "x_fraction", "status"]),
                               (args.excluded, excluded, ["sequence_id", "length", "length_ratio", "x_fraction", "reason"])):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with Path(path).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"included": len(output), "excluded": len(excluded), "reference_length": ref_len}))


def cmd_stats(args: argparse.Namespace) -> None:
    rows = []
    for file_name in args.files:
        records = read_fasta(file_name)
        lengths = [len(r["seq"]) for r in records]
        rows.append({"file": file_name, "num_sequences": len(records), "min_length": min(lengths),
                     "max_length": max(lengths), "mean_length": round(sum(lengths) / len(lengths), 2), "sha256": sha256(file_name)})
    print(json.dumps(rows, ensure_ascii=False, indent=2))


def cmd_qc(args: argparse.Namespace) -> None:
    records = read_fasta(args.alignment)
    lengths = {len(r["seq"]) for r in records}
    if len(lengths) != 1:
        raise ValueError(f"Alignment sequences do not share a length: {sorted(lengths)}")
    aln_len = lengths.pop(); nseq = len(records)
    seq_rows = []
    for r in records:
        gaps = sum(ch in GAP_CHARS for ch in r["seq"])
        seq_rows.append({"sequence_id": r["id"], "alignment_length": aln_len,
                         "ungapped_length": aln_len - gaps, "gap_fraction": round(gaps / aln_len, 6)})
    column_rows = []
    gap_fracs = []
    for col in range(aln_len):
        chars = [r["seq"][col] for r in records]
        gaps = sum(ch in GAP_CHARS for ch in chars)
        frac = gaps / nseq; gap_fracs.append(frac)
        nongap = Counter(ch for ch in chars if ch not in GAP_CHARS)
        column_rows.append({"alignment_column": col + 1, "gap_fraction": round(frac, 6),
                            "nongap_count": nseq - gaps, "consensus": nongap.most_common(1)[0][0] if nongap else "-"})
    summary = {"alignment": str(args.alignment), "sequence_count": nseq, "alignment_length": aln_len,
               "columns_gap_ge_0_5": sum(x >= 0.5 for x in gap_fracs),
               "columns_gap_ge_0_8": sum(x >= 0.8 for x in gap_fracs),
               "mean_column_gap_fraction": round(sum(gap_fracs) / aln_len, 6),
               "sequences_gap_ge_0_5": sum(float(r["gap_fraction"]) >= 0.5 for r in seq_rows),
               "sha256": sha256(args.alignment)}
    Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_json).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for path, rows, fields in ((args.sequence_tsv, seq_rows, ["sequence_id", "alignment_length", "ungapped_length", "gap_fraction"]),
                               (args.column_tsv, column_rows, ["alignment_column", "gap_fraction", "nongap_count", "consensus"])):
        with Path(path).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    print(json.dumps(summary))


def cmd_map(args: argparse.Namespace) -> None:
    records = read_fasta(args.alignment)
    matches = [r for r in records if r["id"] == args.reference_id or r["header"].startswith(args.reference_id + " ")]
    if len(matches) != 1:
        raise ValueError(f"Expected one reference {args.reference_id!r}; found {len(matches)}")
    ref = matches[0]
    expected = {}
    for item in args.expected:
        pos, residue = item.split(":", 1); expected[int(pos)] = residue.upper()
    rows = []
    for position in args.positions:
        ungapped = 0; msa_col = None; observed = None
        for col, residue in enumerate(ref["seq"], 1):
            if residue not in GAP_CHARS:
                ungapped += 1
            if ungapped == position and residue not in GAP_CHARS:
                msa_col, observed = col, residue
                break
        exp = expected.get(position, "")
        status = "pass" if msa_col and (not exp or observed == exp) else ("residue_mismatch" if msa_col else "not_mapped")
        rows.append({"reference_id": args.reference_id, "reference_position": position, "expected_residue": exp,
                     "observed_residue": observed or "", "msa_column": msa_col or "", "status": status})
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
        fields = ["reference_id", "reference_position", "expected_residue", "observed_residue", "msa_column", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"mapped": len(rows), "passed": sum(r["status"] == "pass" for r in rows), "rows": rows}))


def smart_gap_threshold(gap_fractions: list[float]) -> float:
    counts = sorted(Counter(round(x, 4) for x in gap_fractions).items(), reverse=True)
    if len(counts) < 2:
        return 1.0
    alignment_length = len(gap_fractions)
    cumulative = 0.0; slopes: list[float] = []
    for index in range(1, len(counts)):
        cumulative += counts[index - 1][1] / alignment_length
        current_fraction = counts[index][0]; previous_fraction = counts[index - 1][0]
        if current_fraction != previous_fraction:
            slopes.append(abs((counts[index][1] / alignment_length) / (current_fraction - previous_fraction)))
    slopes = slopes[: len(slopes) // 2]
    if not slopes:
        return counts[0][0]
    if len(slopes) == 1:
        return counts[0][0]
    differences = [abs(slopes[i + 1] - slopes[i]) for i in range(len(slopes) - 1)]
    max_index = max(range(len(differences)), key=differences.__getitem__)
    return counts[max_index][0]


def cmd_trim(args: argparse.Namespace) -> None:
    records = read_fasta(args.input)
    lengths = {len(r["seq"]) for r in records}
    if len(lengths) != 1:
        raise ValueError("Input is not an alignment")
    aln_len = lengths.pop(); nseq = len(records)
    fractions = [sum(r["seq"][col] in GAP_CHARS for r in records) / nseq for col in range(aln_len)]
    threshold = smart_gap_threshold(fractions)
    keep = [col for col, frac in enumerate(fractions) if frac <= threshold]
    trimmed = [{"header": r["header"], "id": r["id"], "seq": "".join(r["seq"][i] for i in keep)} for r in records]
    write_fasta(trimmed, args.output)
    report = {"mode": "smart-gap-compatible", "algorithm_source": "ClipKIT 2.13.2 smart_gap_helper.py",
              "gap_characters": "Xx-?*", "selected_threshold": threshold, "input_columns": aln_len,
              "retained_columns": len(keep), "trimmed_columns": aln_len - len(keep), "output_sha256": sha256(args.output)}
    if args.report_json:
        Path(args.report_json).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("compare"); p.add_argument("--fasta-a", required=True); p.add_argument("--fasta-b", required=True)
    p.add_argument("--output-json", required=True); p.add_argument("--output-tsv"); p.set_defaults(func=cmd_compare)
    p = sub.add_parser("filter-blast"); p.add_argument("--input", required=True); p.add_argument("--output-tsv", required=True)
    p.add_argument("--output-accessions", required=True); p.add_argument("--min-identity", type=float, default=25)
    p.add_argument("--max-identity", type=float, default=95); p.add_argument("--min-qcov", type=float, default=70)
    p.add_argument("--max-evalue", type=float, default=1e-20); p.set_defaults(func=cmd_filter_blast)
    p = sub.add_parser("fetch"); p.add_argument("--accessions", required=True); p.add_argument("--output", required=True)
    p.add_argument("--email"); p.add_argument("--batch-size", type=int, default=100); p.set_defaults(func=cmd_fetch)
    p = sub.add_parser("filter-fasta"); p.add_argument("--reference", required=True); p.add_argument("--hits", required=True)
    p.add_argument("--reference-id", required=True); p.add_argument("--output", required=True); p.add_argument("--metadata", required=True)
    p.add_argument("--excluded", required=True); p.add_argument("--min-length-ratio", type=float, default=0.75)
    p.add_argument("--max-length-ratio", type=float, default=1.25); p.add_argument("--max-x-fraction", type=float, default=0.01)
    p.set_defaults(func=cmd_filter_fasta)
    p = sub.add_parser("stats"); p.add_argument("files", nargs="+"); p.set_defaults(func=cmd_stats)
    p = sub.add_parser("qc"); p.add_argument("--alignment", required=True); p.add_argument("--summary-json", required=True)
    p.add_argument("--sequence-tsv", required=True); p.add_argument("--column-tsv", required=True); p.set_defaults(func=cmd_qc)
    p = sub.add_parser("map"); p.add_argument("--alignment", required=True); p.add_argument("--reference-id", required=True)
    p.add_argument("--positions", nargs="+", type=int, required=True); p.add_argument("--expected", nargs="*", default=[])
    p.add_argument("--output", required=True); p.set_defaults(func=cmd_map)
    p = sub.add_parser("trim"); p.add_argument("--input", required=True); p.add_argument("--output", required=True)
    p.add_argument("--report-json"); p.set_defaults(func=cmd_trim)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
