#!/usr/bin/env python3
"""Clean curator-approved TxtC candidate FASTA; do NOT claim orthology from length QC."""
import argparse
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

AA = set('ACDEFGHIKLMNPQRSTVWYBXZJUO')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--candidates', required=True, help='Raw curated candidate FASTA (not raw BLAST table)')
    ap.add_argument('--output', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--min-length-ratio', type=float, default=0.70)
    ap.add_argument('--max-length-ratio', type=float, default=1.30)
    a = ap.parse_args()
    ref_list = list(SeqIO.parse(a.reference, 'fasta'))
    if len(ref_list) != 1:
        raise SystemExit('Reference FASTA must contain exactly one sequence')
    ref = ref_list[0]
    rs = str(ref.seq).upper().strip().removesuffix('*')
    if not rs or set(rs) - AA:
        raise SystemExit('Reference has invalid amino acid symbols; inspect manually')
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    Path(a.manifest).parent.mkdir(parents=True, exist_ok=True)
    records = [SeqRecord(Seq(rs), id='TXT_C_REF_Q8VS75', description='TxtC reference UniProt Q8VS75')]
    seen = {rs: 'TXT_C_REF_Q8VS75'}
    trace = [{'source_id': ref.id, 'assigned_id': 'TXT_C_REF_Q8VS75', 'length': len(rs), 'status': 'reference', 'duplicate_of': ''}]
    for rec in SeqIO.parse(a.candidates, 'fasta'):
        seq = str(rec.seq).upper().replace(' ', '').removesuffix('*')
        rid = rec.id
        reason = 'accepted'
        dup = ''
        if not seq or (set(seq) - AA):
            reason = 'invalid_letters_or_internal_stop'
        elif len(seq) < a.min_length_ratio * len(rs) or len(seq) > a.max_length_ratio * len(rs):
            reason = 'outside_length_window'
        elif seq in seen:
            reason, dup = 'exact_duplicate', seen[seq]
        if reason == 'accepted':
            assigned = f'TXTC_HOM_{len(records):06d}'
            records.append(SeqRecord(Seq(seq), id=assigned, description=''))
            seen[seq] = assigned
        else:
            assigned = ''
        trace.append({'source_id': rid, 'assigned_id': assigned, 'length': len(seq), 'status': reason, 'duplicate_of': dup})
    SeqIO.write(records, a.output, 'fasta')
    with open(a.manifest, 'w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=list(trace[0]))
        writer.writeheader(); writer.writerows(trace)
    print(f'Reference length: {len(rs)} aa; kept including reference: {len(records)}; candidate rejects: {len(trace)-len(records)}')
    print('Note: sequence length filtering does NOT establish TxtC orthology. Curate BLAST hits and BGC context.')

if __name__ == '__main__':
    main()
