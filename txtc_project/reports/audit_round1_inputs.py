"""Read-only input integrity audit; no biological modelling or ranking."""
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'txtc_project/reports/round1_input_audit_20261110'
OUT.mkdir(exist_ok=True)
if (OUT / 'input_audit.json').exists():
    raise FileExistsError('Audit output already exists')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def fasta(path):
    records = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            records.append([line[1:].split()[0], ''])
        elif line.strip():
            if not records:
                raise ValueError('Sequence before header')
            records[-1][1] += line.strip()
    return records

base = ROOT / 'txtc_project/data'
ref_a = fasta(base / 'raw/ncbi/AAW49320.1_TxtC_Streptomyces_turgidiscabies_Car8.fasta')[0][1]
ref_b = fasta(base / 'raw/uniprot/Q8VS75_TxtC_UniProt.fasta')[0][1]
paths = {
    'pilot': base / 'processed/msa/TxtC/TxtC_full_msa.faa',
    'prefilter': base / 'processed/homologs/TxtC_expanded/TxtC_expanded_prefilter.faa',
    'filtered': base / 'processed/homologs/TxtC_expanded/TxtC_expanded_filtered.faa',
}
stats = []
for name, path in paths.items():
    records = fasta(path)
    sequences = {s for _, s in records}
    ungapped = {s.replace('-', '') for s in sequences}
    stats.append(dict(dataset=name, path=str(path.relative_to(ROOT)),
        record_count=len(records), unique_sequence_count=len(sequences),
        unique_id_count=len({i for i, _ in records}),
        unexpected_character_records=sum(bool(set(s) - set('ACDEFGHIKLMNPQRSTVWYXBZJUO-')) for _, s in records),
        reference_matching_rule='Exact amino-acid sequence after removing alignment gaps',
        contains_exact_bgc_reference=ref_a in ungapped,
        contains_exact_structure_reference=ref_b in ungapped, sha256=digest(path)))
cluster = ROOT / 'metadata/coevolution/TxtC/TxtC_seqid95_clusters.tsv'
rows = [r.split('\t') for r in cluster.read_text().splitlines() if r.strip()]
assert all(len(r) == 2 for r in rows)
zip_path = Path('C:/Users/31441/Downloads/TxtC_round1_support_v1.zip')
manual = zip_path.with_name('TxtC_round1_operation_manual_v1.docx')
checks = []
with zipfile.ZipFile(zip_path) as z:
    prefix = 'TxtC_round1_support_v1/'
    for row in z.read(prefix + 'SHA256SUMS.txt').decode().splitlines():
        if not row.strip():
            continue
        expected, filename = row.split(maxsplit=1)
        filename = filename.lstrip('*')
        actual = hashlib.sha256(z.read(prefix + filename)).hexdigest()
        checks.append(dict(file=filename, expected_sha256=expected, actual_sha256=actual, matches=expected == actual))
audit = dict(audit_date='2026-11-10', timezone='Asia/Tokyo',
    scope='Non-functional input integrity audit only',
    repository_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
    manual_baseline_commit='042e7ec479d9940839e50f5812efae4c156efcf0',
    datasets=stats, clusters=dict(rows=len(rows), representatives=len({r[0] for r in rows}), sha256=digest(cluster)),
    support_zip_sha256=digest(zip_path), manual_sha256=digest(manual),
    support_manifest_checks=checks, executed_analysis_steps=[],
    completion_status='Input audit complete; coevolution analysis not performed')
(OUT / 'input_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(audit, ensure_ascii=False, indent=2))
