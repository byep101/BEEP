#!/usr/bin/env python3
"""Minimal, auditable TamI MIp/APC pilot workflow.

This implementation follows TamI_MIp_DCA_pilot_Manual_v1.0.  It deliberately
uses only the Python standard library so that the prepare and MIp stages can be
run on a managed Windows workstation.  DCA remains a separate EVcouplings/plmc
stage and is never approximated here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path


AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
EXPECTED_REFERENCE_LENGTH = 413
EXPECTED_REFERENCE_SHA256 = "fbc563a7aeebba572bac52607531eb227351b690037bf8108f176a3f8701027a"
KNOWN_SITES = {101: "L", 244: "L", 295: "L"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def fasta_id(header: str) -> str:
    return header.split(None, 1)[0]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before first header at {path}:{line_number}")
                chunks.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records found: {path}")
    ids = [fasta_id(header) for header, _ in records]
    duplicates = sorted(identifier for identifier, count in Counter(ids).items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate FASTA IDs: {', '.join(duplicates)}")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")


def write_json(path: Path, data: object) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_weights(path: Path) -> list[float]:
    values: list[float] = []
    with path.open(encoding="ascii") as handle:
        for line_number, raw in enumerate(handle, 1):
            value = raw.strip()
            if not value:
                continue
            try:
                weight = float(value)
            except ValueError as exc:
                raise ValueError(f"Invalid weight at {path}:{line_number}: {value}") from exc
            if not math.isfinite(weight) or weight <= 0:
                raise ValueError(f"Weight must be finite and positive at {path}:{line_number}")
            values.append(weight)
    return values


def compute_weights(
    sequences: list[str], identity_threshold: float = 0.8, overlap_fraction: float = 0.7
) -> tuple[list[float], list[int], int]:
    if not sequences:
        raise ValueError("Cannot weight an empty alignment")
    length = len(sequences[0])
    if any(len(sequence) != length for sequence in sequences):
        raise ValueError("Weighting requires equal-length sequences")
    minimum_overlap = math.ceil(overlap_fraction * length)
    neighbors = [1] * len(sequences)
    for first in range(len(sequences)):
        for second in range(first + 1, len(sequences)):
            comparable = 0
            matches = 0
            for residue_first, residue_second in zip(sequences[first], sequences[second]):
                if residue_first in AA and residue_second in AA:
                    comparable += 1
                    if residue_first == residue_second:
                        matches += 1
            if comparable >= minimum_overlap and matches / comparable >= identity_threshold:
                neighbors[first] += 1
                neighbors[second] += 1
    return [1.0 / count for count in neighbors], neighbors, minimum_overlap


def weighted_mi(first: list[str], second: list[str], weights: list[float]) -> tuple[float, float]:
    joint: Counter[tuple[str, str]] = Counter()
    marginal_first: Counter[str] = Counter()
    marginal_second: Counter[str] = Counter()
    observable_weight = 0.0
    for residue_first, residue_second, weight in zip(first, second, weights):
        if residue_first in AA and residue_second in AA:
            joint[(residue_first, residue_second)] += weight
            marginal_first[residue_first] += weight
            marginal_second[residue_second] += weight
            observable_weight += weight
    if observable_weight == 0:
        return 0.0, 0.0
    score = 0.0
    for (residue_first, residue_second), count in joint.items():
        pair_frequency = count / observable_weight
        first_frequency = marginal_first[residue_first] / observable_weight
        second_frequency = marginal_second[residue_second] / observable_weight
        score += pair_frequency * math.log(pair_frequency / (first_frequency * second_frequency))
    return score, observable_weight


def weighted_entropy(column: list[str], weights: list[float]) -> tuple[float, float]:
    counts: Counter[str] = Counter()
    observable_weight = 0.0
    for residue, weight in zip(column, weights):
        if residue in AA:
            counts[residue] += weight
            observable_weight += weight
    if observable_weight == 0:
        return 0.0, 0.0
    entropy = -sum((count / observable_weight) * math.log(count / observable_weight) for count in counts.values())
    return entropy, observable_weight


def prepare(args: argparse.Namespace) -> None:
    msa_path = Path(args.msa).resolve()
    reference_path = Path(args.reference).resolve()
    output = Path(args.out).resolve()
    output.mkdir(parents=True, exist_ok=True)

    alignment = read_fasta(msa_path)
    aligned_lengths = {len(sequence) for _, sequence in alignment}
    if len(aligned_lengths) != 1:
        raise ValueError(f"MSA sequences are not equal length: {sorted(aligned_lengths)}")
    full_length = aligned_lengths.pop()
    by_id = {fasta_id(header): (index, header, sequence) for index, (header, sequence) in enumerate(alignment)}
    if args.ref_id not in by_id:
        raise ValueError(f"Reference ID {args.ref_id!r} is absent from MSA")
    reference_index, _, aligned_reference = by_id[args.ref_id]

    reference_records = read_fasta(reference_path)
    if len(reference_records) != 1:
        raise ValueError("Reference FASTA must contain exactly one sequence")
    reference_sequence = reference_records[0][1].replace("-", "")
    if any(residue not in AA for residue in reference_sequence):
        bad = sorted(set(reference_sequence) - AA)
        raise ValueError(f"Reference FASTA contains non-standard residues: {''.join(bad)}")
    ungapped_alignment_reference = aligned_reference.replace("-", "")
    if ungapped_alignment_reference != reference_sequence:
        raise ValueError("Ungapped MSA reference does not exactly match reference FASTA")
    if len(reference_sequence) != EXPECTED_REFERENCE_LENGTH:
        raise ValueError(
            f"TamI reference length is {len(reference_sequence)}, expected {EXPECTED_REFERENCE_LENGTH}"
        )
    reference_sha = sha256_text(reference_sequence)
    if reference_sha != EXPECTED_REFERENCE_SHA256:
        raise ValueError(
            f"TamI reference SHA256 {reference_sha} does not match recorded baseline "
            f"{EXPECTED_REFERENCE_SHA256}"
        )
    for position, expected in KNOWN_SITES.items():
        observed = reference_sequence[position - 1]
        if observed != expected:
            raise ValueError(f"Reference position {position} is {observed}, expected {expected}")

    focus_columns = [column for column, residue in enumerate(aligned_reference) if residue != "-"]
    focus_records = [
        (fasta_id(header), "".join(sequence[column] for column in focus_columns))
        for header, sequence in alignment
    ]
    write_fasta(output / "focus_all.faa", focus_records)

    sequence_rows: list[dict[str, object]] = []
    coverage_pass: list[tuple[int, str, str]] = []
    for input_order, (identifier, sequence) in enumerate(focus_records, 1):
        standard_count = sum(residue in AA for residue in sequence)
        coverage = standard_count / len(sequence)
        ambiguity_count = sum(residue not in AA and residue != "-" for residue in sequence)
        status = "coverage_pass" if coverage >= args.coverage else "excluded_low_coverage"
        sequence_rows.append(
            {
                "input_order": input_order,
                "model_order": "",
                "sequence_id": identifier,
                "is_reference": int(input_order - 1 == reference_index),
                "reference_coverage": f"{coverage:.9f}",
                "standard_residue_count": standard_count,
                "gap_count": sequence.count("-"),
                "ambiguity_count": ambiguity_count,
                "status": status,
                "duplicate_of": "",
                "neighbor_count": "",
                "weight": "",
            }
        )
        if status == "coverage_pass":
            coverage_pass.append((input_order - 1, identifier, sequence))

    seen: dict[str, str] = {}
    retained: list[tuple[int, str, str]] = []
    for original_index, identifier, sequence in coverage_pass:
        if sequence in seen:
            sequence_rows[original_index]["status"] = "excluded_exact_focus_duplicate"
            sequence_rows[original_index]["duplicate_of"] = seen[sequence]
        else:
            seen[sequence] = identifier
            retained.append((original_index, identifier, sequence))
            sequence_rows[original_index]["status"] = "retained"
    if not retained:
        raise ValueError("No sequences remain after coverage and duplicate filters")
    if not any(original_index == reference_index for original_index, _, _ in retained):
        raise ValueError("Reference sequence was unexpectedly removed")

    model_positions: list[int] = []
    for focus_index in range(len(focus_columns)):
        occupancy = sum(sequence[focus_index] in AA for _, _, sequence in retained) / len(retained)
        if occupancy >= args.occupancy:
            model_positions.append(focus_index)
    if len(model_positions) < 2:
        raise ValueError("Fewer than two model columns remain after occupancy filtering")

    model_records = [
        (identifier, "".join(sequence[index] for index in model_positions))
        for _, identifier, sequence in retained
    ]
    write_fasta(output / "focus_cov70.faa", model_records)
    weights, neighbors, minimum_overlap = compute_weights(
        [sequence for _, sequence in model_records], args.identity, args.overlap
    )
    with (output / "weights.txt").open("w", encoding="ascii", newline="\n") as handle:
        for weight in weights:
            handle.write(f"{weight:.17g}\n")
    for model_order, ((original_index, _, _), neighbor_count, weight) in enumerate(
        zip(retained, neighbors, weights), 1
    ):
        sequence_rows[original_index]["model_order"] = model_order
        sequence_rows[original_index]["neighbor_count"] = neighbor_count
        sequence_rows[original_index]["weight"] = f"{weight:.17g}"
    write_tsv(
        output / "sequence_audit.tsv",
        [
            "input_order",
            "model_order",
            "sequence_id",
            "is_reference",
            "reference_coverage",
            "standard_residue_count",
            "gap_count",
            "ambiguity_count",
            "status",
            "duplicate_of",
            "neighbor_count",
            "weight",
        ],
        sequence_rows,
    )

    focus_to_model = {focus_index: model_index + 1 for model_index, focus_index in enumerate(model_positions)}
    column_rows: list[dict[str, object]] = []
    retained_sequences = [sequence for _, _, sequence in retained]
    for focus_index, full_column in enumerate(focus_columns):
        occupancy = sum(sequence[focus_index] in AA for sequence in retained_sequences) / len(retained_sequences)
        column_rows.append(
            {
                "reference_position": focus_index + 1,
                "reference_residue": reference_sequence[focus_index],
                "full_msa_column": full_column + 1,
                "model_column": focus_to_model.get(focus_index, ""),
                "occupancy": f"{occupancy:.9f}",
                "status": "retained" if focus_index in focus_to_model else "excluded_low_occupancy",
                "known_site": int(focus_index + 1 in KNOWN_SITES),
            }
        )
    write_tsv(
        output / "column_map.tsv",
        [
            "reference_position",
            "reference_residue",
            "full_msa_column",
            "model_column",
            "occupancy",
            "status",
            "known_site",
        ],
        column_rows,
    )

    neff = sum(weights)
    known_site_status = {
        str(position): {
            "expected_residue": residue,
            "observed_residue": reference_sequence[position - 1],
            "retained": (position - 1) in focus_to_model,
            "model_column": focus_to_model.get(position - 1),
            "full_msa_column": focus_columns[position - 1] + 1,
        }
        for position, residue in KNOWN_SITES.items()
    }
    warnings: list[str] = []
    if neff / len(model_positions) < 1.0:
        warnings.append("low_depth: Neff/L_model < 1; treat this as a method pilot, not a reliability guarantee")
    qc = {
        "workflow": "TamI MIp/DCA pilot v1.0 - prepare",
        "input_msa": str(msa_path),
        "input_reference": str(reference_path),
        "reference_id": args.ref_id,
        "input_msa_sha256": sha256_file(msa_path),
        "input_reference_file_sha256": sha256_file(reference_path),
        "reference_sequence_sha256": reference_sha,
        "N_input": len(alignment),
        "L_full_msa": full_length,
        "L_reference": len(reference_sequence),
        "N_coverage_pass": len(coverage_pass),
        "N_model": len(model_records),
        "L_model": len(model_positions),
        "coverage_threshold": args.coverage,
        "occupancy_threshold": args.occupancy,
        "identity_threshold": args.identity,
        "minimum_pairwise_overlap_fraction": args.overlap,
        "minimum_pairwise_overlap_positions": minimum_overlap,
        "Neff": neff,
        "Neff_over_L_model": neff / len(model_positions),
        "weights_sha256": sha256_file(output / "weights.txt"),
        "focus_alignment_sha256": sha256_file(output / "focus_cov70.faa"),
        "known_sites": known_site_status,
        "warnings": warnings,
    }
    write_json(output / "input_qc.json", qc)
    print(
        f"Prepared {output}: N={len(model_records)}, L_model={len(model_positions)}, "
        f"Neff={neff:.6f}, Neff/L={neff / len(model_positions):.6f}"
    )


def load_column_map(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    retained = [row for row in rows if row["status"] == "retained"]
    retained.sort(key=lambda row: int(row["model_column"]))
    for expected, row in enumerate(retained, 1):
        if int(row["model_column"]) != expected:
            raise ValueError("column_map.tsv model columns are not contiguous")
    return retained


def run_mi(args: argparse.Namespace) -> None:
    run_dir = Path(args.run).resolve()
    alignment_path = run_dir / "focus_cov70.faa"
    weights_path = run_dir / "weights.txt"
    mapping_path = run_dir / "column_map.tsv"
    qc_path = run_dir / "input_qc.json"
    for required in (alignment_path, weights_path, mapping_path, qc_path):
        if not required.exists():
            raise FileNotFoundError(f"Required prepare output is missing: {required}")
    alignment = read_fasta(alignment_path)
    sequences = [sequence for _, sequence in alignment]
    weights = load_weights(weights_path)
    mapping = load_column_map(mapping_path)
    if len(sequences) != len(weights):
        raise ValueError("Sequence and weight counts differ")
    lengths = {len(sequence) for sequence in sequences}
    if lengths != {len(mapping)}:
        raise ValueError(f"Alignment/model map length mismatch: alignment={lengths}, map={len(mapping)}")
    length = len(mapping)
    columns = [[sequence[index] for sequence in sequences] for index in range(length)]
    entropies: list[tuple[float, float]] = [weighted_entropy(column, weights) for column in columns]
    variable = [entropy > 1e-12 for entropy, _ in entropies]

    pairs: list[tuple[int, int, float, float]] = []
    row_sums = [0.0] * length
    total_mi = 0.0
    for first in range(length):
        for second in range(first + 1, length):
            mi_score, observable_weight = weighted_mi(columns[first], columns[second], weights)
            if not math.isfinite(mi_score):
                raise ValueError(f"Non-finite MI for model columns {first + 1}, {second + 1}")
            pairs.append((first, second, mi_score, observable_weight))
            row_sums[first] += mi_score
            row_sums[second] += mi_score
            total_mi += mi_score
    row_means = [value / (length - 1) for value in row_sums]
    overall_mean = (2.0 * total_mi) / (length * (length - 1))

    pair_path = run_dir / "MIp_pairs.csv"
    with pair_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "model_column_i",
            "reference_position_i",
            "reference_residue_i",
            "full_msa_column_i",
            "model_column_j",
            "reference_position_j",
            "reference_residue_j",
            "full_msa_column_j",
            "MI",
            "MIp",
            "observable_weight",
            "variable_i",
            "variable_j",
            "reference_separation",
            "eligible_nonlocal_variable",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for first, second, mi_score, observable_weight in pairs:
            correction = 0.0 if overall_mean == 0 else row_means[first] * row_means[second] / overall_mean
            mip = mi_score - correction
            if not math.isfinite(mip):
                raise ValueError(f"Non-finite MIp for model columns {first + 1}, {second + 1}")
            first_map = mapping[first]
            second_map = mapping[second]
            separation = abs(int(first_map["reference_position"]) - int(second_map["reference_position"]))
            eligible = variable[first] and variable[second] and separation >= args.minimum_separation
            writer.writerow(
                {
                    "model_column_i": first + 1,
                    "reference_position_i": first_map["reference_position"],
                    "reference_residue_i": first_map["reference_residue"],
                    "full_msa_column_i": first_map["full_msa_column"],
                    "model_column_j": second + 1,
                    "reference_position_j": second_map["reference_position"],
                    "reference_residue_j": second_map["reference_residue"],
                    "full_msa_column_j": second_map["full_msa_column"],
                    "MI": f"{mi_score:.17g}",
                    "MIp": f"{mip:.17g}",
                    "observable_weight": f"{observable_weight:.17g}",
                    "variable_i": int(variable[first]),
                    "variable_j": int(variable[second]),
                    "reference_separation": separation,
                    "eligible_nonlocal_variable": int(eligible),
                }
            )

    with pair_path.open(encoding="utf-8", newline="") as handle:
        eligible_rows = [
            row for row in csv.DictReader(handle) if row["eligible_nonlocal_variable"] == "1"
        ]
    eligible_rows.sort(
        key=lambda row: (
            -float(row["MIp"]),
            int(row["reference_position_i"]),
            int(row["reference_position_j"]),
        )
    )
    top_pairs_path = run_dir / "MIp_top20_pairs.csv"
    with top_pairs_path.open("w", encoding="utf-8", newline="") as handle:
        top_fieldnames = ["pair_rank"] + fieldnames
        writer = csv.DictWriter(handle, fieldnames=top_fieldnames, lineterminator="\n")
        writer.writeheader()
        for rank, row in enumerate(eligible_rows[:20], 1):
            writer.writerow({"pair_rank": rank, **row})

    site_best: dict[int, tuple[float, int, str, int]] = {}
    for pair_rank, row in enumerate(eligible_rows, 1):
        score = float(row["MIp"])
        for side, other in (("i", "j"), ("j", "i")):
            position = int(row[f"reference_position_{side}"])
            candidate = (
                score,
                int(row[f"reference_position_{other}"]),
                row[f"reference_residue_{other}"],
                pair_rank,
            )
            if position not in site_best or candidate[0] > site_best[position][0]:
                site_best[position] = candidate
    site_ranks = {
        position: 1 + sum(1 for other in site_best.values() if other[0] > best[0])
        for position, best in site_best.items()
    }
    residue_rank_rows: list[dict[str, object]] = []
    for model_index, row in enumerate(mapping, 1):
        position = int(row["reference_position"])
        best = site_best.get(position)
        rank = site_ranks.get(position)
        residue_rank_rows.append(
            {
                "reference_position": position,
                "reference_residue": row["reference_residue"],
                "model_column": model_index,
                "variable": int(variable[model_index - 1]),
                "MIp_residue_rank": rank if rank is not None else "",
                "MIp_max_score": f"{best[0]:.17g}" if best else "",
                "best_partner_position": best[1] if best else "",
                "best_partner_residue": best[2] if best else "",
                "best_pair_rank": best[3] if best else "",
                "MIp_residue_top20": int(rank is not None and rank <= 20),
                "known_site": int(position in KNOWN_SITES),
            }
        )
    residue_ranks_path = run_dir / "MIp_residue_ranks.csv"
    with residue_ranks_path.open("w", encoding="utf-8", newline="") as handle:
        rank_fieldnames = [
            "reference_position",
            "reference_residue",
            "model_column",
            "variable",
            "MIp_residue_rank",
            "MIp_max_score",
            "best_partner_position",
            "best_partner_residue",
            "best_pair_rank",
            "MIp_residue_top20",
            "known_site",
        ]
        writer = csv.DictWriter(handle, fieldnames=rank_fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(residue_rank_rows)

    entropy_rows: list[dict[str, object]] = []
    for index, ((entropy, observable_weight), row) in enumerate(zip(entropies, mapping), 1):
        entropy_rows.append(
            {
                "model_column": index,
                "reference_position": row["reference_position"],
                "reference_residue": row["reference_residue"],
                "full_msa_column": row["full_msa_column"],
                "entropy": f"{entropy:.17g}",
                "observable_weight": f"{observable_weight:.17g}",
                "variable": int(variable[index - 1]),
            }
        )
    write_tsv(
        run_dir / "site_entropy.tsv",
        [
            "model_column",
            "reference_position",
            "reference_residue",
            "full_msa_column",
            "entropy",
            "observable_weight",
            "variable",
        ],
        entropy_rows,
    )
    summary = {
        "workflow": "TamI MIp/DCA pilot v1.0 - MIp/APC",
        "N": len(sequences),
        "L_model": length,
        "pair_count": len(pairs),
        "expected_pair_count": length * (length - 1) // 2,
        "finite_scores": True,
        "mi_log_base": "natural",
        "pseudocount": 0,
        "gap_treatment": "pairwise missing; gaps and non-standard residues excluded, counts renormalized",
        "apc_scope": "all model column pairs; diagonal excluded",
        "minimum_reference_separation": args.minimum_separation,
        "eligible_rule": "both sites variable and absolute reference-position separation >= 6",
        "focus_alignment_sha256": sha256_file(alignment_path),
        "weights_sha256": sha256_file(weights_path),
        "MIp_pairs_sha256": sha256_file(pair_path),
        "overall_mean_MI": overall_mean,
        "variable_site_count": sum(variable),
        "constant_site_count": length - sum(variable),
        "eligible_pair_count": len(eligible_rows),
        "MIp_top20_pairs_sha256": sha256_file(top_pairs_path),
        "MIp_residue_ranks_sha256": sha256_file(residue_ranks_path),
        "known_site_MIp": {
            str(position): {
                "retained": any(
                    int(row["reference_position"]) == position for row in mapping
                ),
                "residue_rank": site_ranks.get(position),
                "max_score": site_best[position][0] if position in site_best else None,
                "best_partner_position": site_best[position][1] if position in site_best else None,
                "best_pair_rank": site_best[position][3] if position in site_best else None,
                "residue_top20": bool(position in site_ranks and site_ranks[position] <= 20),
            }
            for position in KNOWN_SITES
        },
    }
    write_json(run_dir / "mi_summary.json", summary)
    print(
        f"MIp complete: L_model={length}, pairs={len(pairs)}, "
        f"variable_sites={sum(variable)}, all_scores_finite=true"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="TamI minimal prepare and MIp/APC workflow")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="validate inputs, filter columns, and calculate weights")
    prepare_parser.add_argument("--msa", required=True)
    prepare_parser.add_argument("--reference", required=True)
    prepare_parser.add_argument("--ref-id", required=True)
    prepare_parser.add_argument("--out", required=True)
    prepare_parser.add_argument("--coverage", type=float, default=0.7)
    prepare_parser.add_argument("--occupancy", type=float, default=0.7)
    prepare_parser.add_argument("--identity", type=float, default=0.8)
    prepare_parser.add_argument("--overlap", type=float, default=0.7)
    prepare_parser.set_defaults(func=prepare)

    mi_parser = subparsers.add_parser("mi", help="calculate weighted MI and APC-corrected MIp")
    mi_parser.add_argument("--run", required=True)
    mi_parser.add_argument("--minimum-separation", type=int, default=6)
    mi_parser.set_defaults(func=run_mi)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    for name in ("coverage", "occupancy", "identity", "overlap"):
        if hasattr(args, name) and not 0 < getattr(args, name) <= 1:
            parser.error(f"--{name} must be in (0, 1]")
    try:
        args.func(args)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
