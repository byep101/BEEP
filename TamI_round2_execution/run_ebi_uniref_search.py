from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEARCH_ROOT = ROOT / "TamI_round2_supplemental_search_v1"
QUERY_DIR = SEARCH_ROOT / "00_queries"
OUT_ROOT = SEARCH_ROOT / "02_uniref"
BASE = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
EMAIL = "byep101@gmail.com"


def get(url: str, timeout: int = 60) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "TamI-round2-research/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def post(url: str, fields: dict[str, str], timeout: int = 60) -> bytes:
    data = urllib.parse.urlencode(fields).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=data,
        headers={"User-Agent": "TamI-round2-research/1.0", "Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.read()
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} from {url}: {body}") from exc


def result_type_ids(xml_bytes: bytes) -> list[str]:
    root = ET.fromstring(xml_bytes)
    ids: list[str] = []
    for element in root.iter():
        if element.tag.split("}")[-1] == "identifier" and element.text:
            ids.append(element.text.strip())
    return list(dict.fromkeys(ids))


def run_one(label: str, database: str, max_hits: int, max_wait_seconds: int) -> dict[str, object]:
    query_path = QUERY_DIR / f"{label}.faa"
    if not query_path.exists():
        raise FileNotFoundError(query_path)
    sequence = query_path.read_text(encoding="ascii")
    out_dir = OUT_ROOT / f"{database}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    record_path = out_dir / "job_record.json"
    existing_record: dict[str, object] | None = None
    if record_path.exists():
        existing_record = json.loads(record_path.read_text(encoding="utf-8"))
        if existing_record.get("final_status") == "FINISHED":
            print(f"{label} {database}: already FINISHED ({existing_record.get('job_id')})", flush=True)
            return existing_record

    submitted = datetime.now(timezone.utc).isoformat()
    fields = {
        "email": EMAIL,
        "title": f"TamI_round2_{label}_{database}",
        "program": "blastp",
        "stype": "protein",
        "database": database,
        "matrix": "BLOSUM62",
        "exp": "1e-50",
        "filter": "F",
        "alignments": str(max_hits),
        "scores": str(max_hits),
        "sequence": sequence,
    }
    if existing_record and existing_record.get("job_id") and existing_record.get("final_status") == "TIMEOUT":
        job_id = str(existing_record["job_id"])
        submitted = str(existing_record.get("submitted_utc", submitted))
        print(f"{label} {database}: resuming {job_id}", flush=True)
    else:
        job_id = post(f"{BASE}/run", fields).decode("utf-8").strip()
        if not job_id or "\n" in job_id or " " in job_id:
            raise RuntimeError(f"Unexpected EBI job id: {job_id!r}")
        print(f"{label} {database}: submitted {job_id}", flush=True)
    start = time.monotonic()
    status_history: list[dict[str, object]] = list(existing_record.get("status_history", [])) if existing_record else []
    final_status = "UNKNOWN"
    while True:
        status = get(f"{BASE}/status/{job_id}").decode("utf-8").strip()
        elapsed = time.monotonic() - start
        if not status_history or status_history[-1]["status"] != status:
            print(f"{label} {database}: {status} at {elapsed:.1f}s", flush=True)
        status_history.append({"elapsed_seconds": round(elapsed, 3), "status": status})
        if status in {"FINISHED", "ERROR", "FAILURE", "NOT_FOUND"}:
            final_status = status
            break
        if elapsed > max_wait_seconds:
            final_status = "TIMEOUT"
            break
        time.sleep(5)

    record: dict[str, object] = {
        "query_label": label,
        "query_file": str(query_path),
        "database": database,
        "program": "blastp",
        "stype": "protein",
        "matrix": "BLOSUM62",
        "evalue": 1e-50,
        "filter": "F",
        "max_hits": max_hits,
        "email": EMAIL,
        "job_id": job_id,
        "submitted_utc": submitted,
        "final_status": final_status,
        "status_history": status_history,
    }
    if final_status == "FINISHED":
        types_xml = get(f"{BASE}/resulttypes/{job_id}")
        (out_dir / "resulttypes.xml").write_bytes(types_xml)
        types = result_type_ids(types_xml)
        record["result_types"] = types
        for result_type in types:
            try:
                data = get(f"{BASE}/result/{job_id}/{result_type}", timeout=120)
            except urllib.error.HTTPError as exc:
                record.setdefault("download_errors", []).append({"result_type": result_type, "http_status": exc.code})
                continue
            (out_dir / f"result.{result_type}").write_bytes(data)
    record_path.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{label} {database}: final {final_status}", flush=True)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default="uniref90", choices=["uniref50", "uniref90", "uniref100"])
    parser.add_argument("--query-label", action="append", help="Repeat for selected labels; default is all six")
    parser.add_argument("--max-hits", type=int, default=250)
    parser.add_argument("--max-wait-seconds", type=int, default=600)
    args = parser.parse_args()
    labels = args.query_label or [f"TamI_branch_{index:02d}" for index in range(1, 7)]
    records: list[dict[str, object]] = []
    for label in labels:
        records.append(run_one(label, args.database, args.max_hits, args.max_wait_seconds))
        if records[-1]["final_status"] != "FINISHED":
            break
    summary = {
        "database": args.database,
        "requested_queries": labels,
        "completed_queries": [record["query_label"] for record in records if record["final_status"] == "FINISHED"],
        "records": records,
    }
    (OUT_ROOT / f"{args.database}_run_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
