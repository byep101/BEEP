from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.error
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "TamI_round2_supplemental_search_v1"
SEARCH = RUN / "02_uniref"
OUT = RUN / "03_uniref_candidates"
BASE = "https://rest.uniprot.org/uniref"
USER_AGENT = "TamI-round2-research/1.0 (byep101@gmail.com)"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_one_fasta(data: bytes, expected_id: str) -> tuple[str, str]:
    lines = data.decode("utf-8").splitlines()
    headers = [line[1:].strip() for line in lines if line.startswith(">")]
    if len(headers) != 1:
        raise ValueError(f"Expected one FASTA record for {expected_id}, received {len(headers)}")
    actual_id = headers[0].split()[0]
    if actual_id != expected_id:
        raise ValueError(f"UniRef response ID mismatch: expected {expected_id}, got {actual_id}")
    sequence = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    sequence = sequence.upper().replace("-", "").replace(".", "")
    if not sequence:
        raise ValueError(f"Empty sequence for {expected_id}")
    return headers[0], sequence


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for accession, description, sequence in records:
            handle.write(f">{accession} {description}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def fetch(accession: str, cache_path: Path) -> tuple[bytes, dict[str, str]]:
    if cache_path.exists():
        return cache_path.read_bytes(), {"retrieval_mode": "cache", "url": f"{BASE}/{accession}.fasta"}
    url = f"{BASE}/{accession}.fasta"
    last_error: Exception | None = None
    for attempt in range(1, 5):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "text/plain"})
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                data = response.read()
                metadata = {
                    "retrieval_mode": "network",
                    "url": url,
                    "http_status": str(response.status),
                    "last_modified": response.headers.get("Last-Modified", ""),
                    "etag": response.headers.get("ETag", ""),
                }
            cache_path.write_bytes(data)
            return data, metadata
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            time.sleep(min(2 ** attempt, 15))
    raise RuntimeError(f"Failed to retrieve {accession}: {last_error}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "raw_fasta"
    cache.mkdir(exist_ok=True)
    long_rows: list[dict[str, object]] = []
    by_accession: dict[str, list[dict[str, object]]] = defaultdict(list)
    result_hashes: dict[str, str] = {}
    for branch in range(1, 7):
        label = f"TamI_branch_{branch:02d}"
        path = SEARCH / f"uniref90_{label}" / "result.tsv"
        if not path.exists():
            raise FileNotFoundError(f"Missing completed UniRef90 result: {path}")
        result_hashes[label] = digest(path)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            raise ValueError(f"No hits in {path}")
        for row in rows:
            item: dict[str, object] = {
                "query_label": label,
                "hit_rank": int(row["Hit"]),
                "database": row["DB"],
                "uniref_accession": row["Accession"],
                "aux_query_description": row["Description"],
                "organism": row["Organism"],
                "reported_length": int(row["Length"]),
                "aux_query_score_bits": float(row["Score(Bits)"]),
                "aux_query_identity_percent": float(row["Identities(%)"]),
                "aux_query_positives_percent": float(row["Positives(%)"]),
                "aux_query_evalue": float(row["E()"]),
            }
            long_rows.append(item)
            by_accession[str(item["uniref_accession"])].append(item)

    ordered_accessions = list(dict.fromkeys(str(row["uniref_accession"]) for row in long_rows))
    records: list[tuple[str, str, str]] = []
    retrieval_rows: list[dict[str, object]] = []
    retrieved_utc = datetime.now(timezone.utc).isoformat()
    for index, accession in enumerate(ordered_accessions, 1):
        cache_path = cache / f"{accession}.fasta"
        data, metadata = fetch(accession, cache_path)
        header, sequence = parse_one_fasta(data, accession)
        hits = by_accession[accession]
        description = header[len(accession):].strip()
        records.append((accession, description, sequence))
        retrieval_rows.append({
            "discovery_order": index,
            "uniref_accession": accession,
            "query_support_count": len({str(hit["query_label"]) for hit in hits}),
            "query_support": ";".join(dict.fromkeys(str(hit["query_label"]) for hit in hits)),
            "best_aux_query_evalue": min(float(hit["aux_query_evalue"]) for hit in hits),
            "max_aux_query_identity_percent": max(float(hit["aux_query_identity_percent"]) for hit in hits),
            "representative_header": header,
            "sequence_length": len(sequence),
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "raw_fasta_sha256": hashlib.sha256(data).hexdigest(),
            "retrieved_utc": retrieved_utc,
            **metadata,
        })
        if metadata["retrieval_mode"] == "network":
            time.sleep(0.12)
        if index % 25 == 0 or index == len(ordered_accessions):
            print(f"retrieved {index}/{len(ordered_accessions)} UniRef90 representatives", flush=True)

    write_tsv(OUT / "uniref90_auxiliary_search_hits.tsv", long_rows)
    write_tsv(OUT / "uniref90_retrieval_audit.tsv", retrieval_rows)
    write_fasta(OUT / "uniref90_representatives_raw.faa", records)
    summary = {
        "created_utc": retrieved_utc,
        "database": "UniRef90",
        "queries_completed": 6,
        "raw_hit_rows": len(long_rows),
        "unique_uniref90_clusters": len(ordered_accessions),
        "clusters_supported_by_multiple_queries": sum(int(row["query_support_count"]) > 1 for row in retrieval_rows),
        "result_tsv_sha256": result_hashes,
        "candidate_fasta_sha256": digest(OUT / "uniref90_representatives_raw.faa"),
        "acceptance_warning": "These are discovery candidates only. Auxiliary-query metrics are not intake QC. Re-align every representative to TamI reference ADC79647.1 and apply the frozen TamI-reference thresholds.",
    }
    (OUT / "uniref90_collection_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
