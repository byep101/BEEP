from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ROUND2 = ROOT / "TamI_round2_expanded_results_v1"
MODEL_FASTA = ROUND2 / "04_prepare" / "final_curated_v2" / "focus_cov70.faa"
FULL_FASTA = ROUND2 / "02_intake" / "combined_curated_expanded_v2.faa"
WEIGHTS = ROUND2 / "05_qc" / "final_sequence_weights_with_origin.tsv"
CANDIDATES = ROUND2 / "02_intake" / "expanded_candidate_accession_audit.tsv"
OUT = ROOT / "TamI_round2_supplemental_search_v1"
QUERY_DIR = OUT / "00_queries"
REFERENCE_ID = "TamI_REF_ADC79647.1"
AA = set("ACDEFGHIKLMNPQRSTVWY")
TARGET_CLUSTERS = 6


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
                raise ValueError(f"Sequence before header in {path}")
            chunks.append("".join(line.split()))
    if header is not None:
        rows.append((header.split()[0], "".join(chunks).upper()))
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, rows: list[tuple[str, str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for label, accession, sequence in rows:
            handle.write(f">{label}|{accession}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def distance(seq_a: str, seq_b: str) -> tuple[float, float, float]:
    shared = [(a, b) for a, b in zip(seq_a, seq_b) if a in AA and b in AA]
    overlap = len(shared) / len(seq_a)
    if not shared or overlap < 0.70:
        return 1.0, 0.0, overlap
    identity = sum(a == b for a, b in shared) / len(shared)
    return 1.0 - identity, identity, overlap


def key(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a < b else (b, a)


def safe_leaf(sequence_id: str) -> str:
    return sequence_id.replace(":", "_").replace(";", "_").replace("(", "_").replace(")", "_").replace(",", "_")


def main() -> None:
    QUERY_DIR.mkdir(parents=True, exist_ok=False)
    (OUT / "01_clusterednr").mkdir()
    (OUT / "02_uniref").mkdir()
    (OUT / "03_reference_qc").mkdir()
    model_by_id = dict(read_fasta(MODEL_FASTA))
    full_by_id = dict(read_fasta(FULL_FASTA))
    weight_rows = read_tsv(WEIGHTS)
    weight_by_id = {row["sequence_id"]: row for row in weight_rows}
    candidate_rows = read_tsv(CANDIDATES)
    candidate_by_id = {row["accession"]: row for row in candidate_rows}

    pool_ids = [REFERENCE_ID]
    excluded_from_query_pool: list[dict[str, object]] = []
    for row in weight_rows:
        sequence_id = row["sequence_id"]
        if sequence_id == REFERENCE_ID or not row["origin"].startswith("new_"):
            continue
        annotation = candidate_by_id[sequence_id]["annotation"]
        lower = annotation.lower()
        if "cytochrome p450" not in lower:
            excluded_from_query_pool.append(
                {"sequence_id": sequence_id, "reason": "query representative requires supplied cytochrome P450 annotation", "annotation": annotation}
            )
            continue
        pool_ids.append(sequence_id)

    if len(pool_ids) < TARGET_CLUSTERS:
        raise AssertionError("Credible representative pool is too small")
    if any(sequence_id not in model_by_id or sequence_id not in full_by_id for sequence_id in pool_ids):
        raise AssertionError("Representative pool ID missing from model or full-length FASTA")

    leaf_distance: dict[tuple[str, str], float] = {}
    leaf_identity_to_ref: dict[str, float] = {REFERENCE_ID: 1.0}
    leaf_overlap_to_ref: dict[str, float] = {REFERENCE_ID: 1.0}
    for index, sequence_id in enumerate(pool_ids):
        for other_id in pool_ids[index + 1:]:
            dist, identity, overlap = distance(model_by_id[sequence_id], model_by_id[other_id])
            leaf_distance[key(sequence_id, other_id)] = dist
            if sequence_id == REFERENCE_ID:
                leaf_identity_to_ref[other_id] = identity
                leaf_overlap_to_ref[other_id] = overlap
            elif other_id == REFERENCE_ID:
                leaf_identity_to_ref[sequence_id] = identity
                leaf_overlap_to_ref[sequence_id] = overlap

    clusters: dict[str, dict[str, object]] = {
        sequence_id: {
            "members": [sequence_id],
            "size": 1,
            "height": 0.0,
            "newick": safe_leaf(sequence_id),
        }
        for sequence_id in pool_ids
    }
    active = list(pool_ids)
    cluster_dist = dict(leaf_distance)
    merge_index = 0
    six_cluster_members: list[list[str]] | None = None
    merge_rows: list[dict[str, object]] = []

    while len(active) > 1:
        if len(active) == TARGET_CLUSTERS and six_cluster_members is None:
            six_cluster_members = [list(clusters[cluster_id]["members"]) for cluster_id in active]
        pair = min((key(a, b) for idx, a in enumerate(active) for b in active[idx + 1:]), key=lambda pair_key: (cluster_dist[pair_key], pair_key))
        left, right = pair
        dist = cluster_dist[pair]
        merge_index += 1
        merged_id = f"UPGMA_{merge_index:04d}"
        left_cluster = clusters[left]
        right_cluster = clusters[right]
        merged_height = dist / 2.0
        left_branch = max(0.0, merged_height - float(left_cluster["height"]))
        right_branch = max(0.0, merged_height - float(right_cluster["height"]))
        merged_members = list(left_cluster["members"]) + list(right_cluster["members"])
        merged_size = int(left_cluster["size"]) + int(right_cluster["size"])
        clusters[merged_id] = {
            "members": merged_members,
            "size": merged_size,
            "height": merged_height,
            "newick": f"({left_cluster['newick']}:{left_branch:.8f},{right_cluster['newick']}:{right_branch:.8f})",
        }
        merge_rows.append(
            {
                "merge_index": merge_index,
                "left_cluster": left,
                "right_cluster": right,
                "distance": dist,
                "height": merged_height,
                "merged_size": merged_size,
            }
        )
        for other in [item for item in active if item not in {left, right}]:
            d_left = cluster_dist[key(left, other)]
            d_right = cluster_dist[key(right, other)]
            merged_distance = (
                int(left_cluster["size"]) * d_left + int(right_cluster["size"]) * d_right
            ) / merged_size
            cluster_dist[key(merged_id, other)] = merged_distance
        active = [item for item in active if item not in {left, right}] + [merged_id]
        for pair_key in [pair_key for pair_key in cluster_dist if left in pair_key or right in pair_key]:
            del cluster_dist[pair_key]

    if six_cluster_members is None or len(six_cluster_members) != TARGET_CLUSTERS:
        raise AssertionError("Failed to capture six UPGMA groups")
    tree_root = clusters[active[0]]
    (QUERY_DIR / "credible_set_upgma_distance_tree.nwk").write_text(str(tree_root["newick"]) + ";\n", encoding="ascii")
    write_tsv(QUERY_DIR / "upgma_merge_history.tsv", merge_rows)
    if excluded_from_query_pool:
        write_tsv(QUERY_DIR / "representative_pool_exclusions.tsv", excluded_from_query_pool)

    representatives: list[dict[str, object]] = []
    membership_rows: list[dict[str, object]] = []
    cluster_infos: list[tuple[int, list[str]]] = list(enumerate(six_cluster_members, 1))
    cluster_infos.sort(key=lambda item: (REFERENCE_ID not in item[1], -len(item[1]), sorted(item[1])[0]))
    for output_index, (_, members) in enumerate(cluster_infos, 1):
        if REFERENCE_ID in members:
            representative = REFERENCE_ID
            selection_reason = "TamI reference fixed as representative of its UPGMA group"
        else:
            def medoid_score(sequence_id: str) -> tuple[float, float, str]:
                if len(members) == 1:
                    average = 0.0
                else:
                    average = sum(leaf_distance[key(sequence_id, other)] for other in members if other != sequence_id) / (len(members) - 1)
                weight = float(weight_by_id[sequence_id]["weight"])
                return average, -weight, sequence_id

            representative = min(members, key=medoid_score)
            selection_reason = "minimum mean project-distance medoid within six-group UPGMA cut; weight used as tie-breaker"
        label = f"TamI_branch_{output_index:02d}"
        candidate = candidate_by_id.get(representative, {})
        representatives.append(
            {
                "query_label": label,
                "accession": "ADC79647.1" if representative == REFERENCE_ID else representative,
                "internal_sequence_id": representative,
                "cluster_size": len(members),
                "project_identity_to_TamI_percent": 100.0 * leaf_identity_to_ref[representative],
                "project_overlap_to_TamI": leaf_overlap_to_ref[representative],
                "model_weight": weight_by_id[representative]["weight"],
                "neighbors_including_self": weight_by_id[representative]["neighbors_including_self"],
                "origin": weight_by_id[representative]["origin"],
                "selected_source": weight_by_id[representative]["selected_source"],
                "full_length_aa": len(full_by_id[representative]),
                "annotation": "TamI reference ADC79647.1" if representative == REFERENCE_ID else candidate.get("annotation", ""),
                "organism": "Streptomyces sp. 307-9" if representative == REFERENCE_ID else candidate.get("organism", ""),
                "selection_reason": selection_reason,
                "full_sequence_sha256": hashlib.sha256(full_by_id[representative].encode("ascii")).hexdigest(),
            }
        )
        for member in sorted(members):
            membership_rows.append(
                {
                    "query_group": label,
                    "group_size": len(members),
                    "representative_sequence_id": representative,
                    "member_sequence_id": member,
                    "member_identity_to_TamI_percent": 100.0 * leaf_identity_to_ref[member],
                    "member_origin": weight_by_id[member]["origin"],
                }
            )

    write_tsv(QUERY_DIR / "representative_selection.tsv", representatives)
    write_tsv(QUERY_DIR / "representative_group_membership.tsv", membership_rows)
    write_fasta(
        QUERY_DIR / "representative_queries.faa",
        [(str(row["query_label"]), str(row["internal_sequence_id"]), full_by_id[str(row["internal_sequence_id"])]) for row in representatives],
    )
    for row in representatives:
        label = str(row["query_label"])
        sequence_id = str(row["internal_sequence_id"])
        write_fasta(QUERY_DIR / f"{label}.faa", [(label, sequence_id, full_by_id[sequence_id])])

    summary = {
        "method": "six-group UPGMA cut on project-defined pairwise model-sequence distances; medoid representative per group; TamI reference fixed in its group",
        "tree_interpretation": "distance-based query diversification aid, not a calibrated phylogeny",
        "credible_pool_definition": "TamI reference plus final new model sequences with supplied cytochrome P450 annotation; all already passed TamI-reference numeric and final three-anchor QC",
        "credible_pool_N": len(pool_ids),
        "excluded_from_query_pool_N": len(excluded_from_query_pool),
        "representative_N": len(representatives),
        "representative_accessions": [row["accession"] for row in representatives],
        "query_fasta_sha256": hashlib.sha256((QUERY_DIR / "representative_queries.faa").read_bytes()).hexdigest(),
        "downstream_acceptance_rule": "auxiliary-query hit is only a discovery candidate; recompute identity, query coverage, length, sequence quality and P450 motif status against TamI reference before intake; remove known accessions and exact sequences",
        "fixed_TamI_reference_thresholds": {
            "identity_percent_inclusive": [25, 95],
            "query_span_coverage_min": 0.70,
            "reported_or_recomputed_query_coverage_min": 0.70,
            "length_ratio_inclusive": [0.75, 1.25],
            "nonstandard_fraction_max": 0.01,
            "three_anchor_review": ["E287", "R290", "C359"],
        },
    }
    (QUERY_DIR / "representative_selection_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(
        "# TamI representative-query supplemental search\n\n"
        "The representative query set was selected from the final TamI round-2 model after TamI-reference numeric QC, full-length MSA, focus deduplication and three-anchor motif QC.\n\n"
        "The UPGMA distance tree is used only to diversify search queries. It is not treated as a calibrated phylogeny or functional proof.\n\n"
        "Every auxiliary-query hit remains a discovery candidate. Before intake, it must be aligned back to TamI reference ADC79647.1 and pass the unchanged TamI-reference identity, coverage, length, ambiguity, exact-deduplication and motif checks.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(json.dumps(representatives, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
