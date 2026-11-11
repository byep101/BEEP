#!/usr/bin/env python3
"""Produce TxtC focus alignment, reference-position map, QC and external weights."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ALPHABET = '-ACDEFGHIKLMNPQRSTVWY'
CANONICAL = set(ALPHABET[1:])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--msa', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--focus', default='TXT_C_REF_Q8VS75')
    ap.add_argument('--min-row-coverage', type=float, default=.70)
    ap.add_argument('--max-column-gap', type=float, default=.50)
    ap.add_argument('--identity-threshold', type=float, default=.80)
    ap.add_argument('--anchors', nargs='*', default=[], help='Optional REF_POSITION:AA e.g. 237:T; do not guess numbering')
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(args.msa, 'fasta'))
    if not records:
        raise SystemExit('Empty MSA')
    names = [r.id for r in records]
    if names.count(args.focus) != 1:
        raise SystemExit(f'Expected focus ID {args.focus} exactly once')
    records = [next(r for r in records if r.id == args.focus)] + [r for r in records if r.id != args.focus]
    seqs = [str(r.seq).upper().replace('.', '-') for r in records]
    if len({len(s) for s in seqs}) != 1:
        raise SystemExit('MSA rows have different lengths')
    ref_aln = seqs[0]
    ref_pos = 0; focus_cols = []
    for col, char in enumerate(ref_aln):
        if char != '-':
            ref_pos += 1
            if char not in CANONICAL:
                raise SystemExit(f'Reference contains noncanonical residue at reference position {ref_pos}: {char}')
            focus_cols.append((col, ref_pos, char))
    if not focus_cols:
        raise SystemExit('No reference residues found')
    row_keep = [True]
    for s in seqs[1:]:
        fraction = sum(s[j] in CANONICAL for j, _, _ in focus_cols)/len(focus_cols)
        row_keep.append(fraction >= args.min_row_coverage)
    filtered = [(rec, s) for rec, s, keep in zip(records, seqs, row_keep) if keep]
    if len(filtered) < 2:
        raise SystemExit('No homologs retained after row coverage QC')
    def norm(c):
        return c if c in CANONICAL else '-'
    kept_cols = []
    for c, rpos, aa in focus_cols:
        gap = sum(norm(seq[c]) == '-' for _,seq in filtered) / len(filtered)
        if gap <= args.max_column_gap:
            kept_cols.append((c,rpos,aa,gap))
    if len(kept_cols) < 20:
        raise SystemExit('Too few model columns after gap filtering; inspect the dataset')
    # Exact focus duplicates are removed AFTER defining the final model columns.
    model_records=[]; seen={}
    for rec, s in filtered:
        ms = ''.join(norm(s[c]) for c,_,_,_ in kept_cols)
        if ms in seen:
            if rec.id == args.focus:
                raise SystemExit('Focus was duplicated and ordered incorrectly')
            continue
        seen[ms] = rec.id
        model_records.append(SeqRecord(Seq(ms), id=rec.id, description=''))
    if model_records[0].id != args.focus:
        raise SystemExit('Focus moved from row 1')
    L = len(kept_cols)
    pos_map = {pos:(i+1, aa) for i, (_,pos,aa,_) in enumerate(kept_cols)}
    for x in args.anchors:
        pos_s, aa = x.split(':')
        pos = int(pos_s); aa = aa.upper()
        if pos not in pos_map or pos_map[pos][1] != aa:
            raise SystemExit(f'Anchor failed ({x}). STOP: check reference / PDB numbering and gap mask, do not force.')
    seq_strs = [str(r.seq) for r in model_records]
    arr = np.array([[ALPHABET.index(c) for c in s] for s in seq_strs], dtype=np.int8)
    n = len(arr)
    # External weights: identical non-gap residues / all model columns >= identity threshold.
    # This explicitly defined weighting is used for BOTH MIp and plmc -w.
    weights = np.zeros(n, dtype=float)
    for i in range(n):
        identical = np.sum((arr == arr[i]) & (arr != 0), axis=1)
        nhood = int(np.sum(identical / L >= args.identity_threshold))
        weights[i] = 1.0 / max(nhood, 1)
    np.savetxt(out / 'weights.txt', weights, fmt='%.12g')
    SeqIO.write(model_records, out / 'txtc_focus.fasta', 'fasta')
    with open(out / 'position_map.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f); writer.writerow(['model_position','ref_position','ref_aa','msa_column_1based','gap_fraction'])
        for i,(c,rpos,aa,gf) in enumerate(kept_cols,1): writer.writerow([i,rpos,aa,c+1,f'{gf:.6f}'])
    with open(out / 'site_qc.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f); writer.writerow(['model_position','ref_position','ref_aa','is_variable','non_gap_fraction'])
        for j,(_,rpos,aa,_) in enumerate(kept_cols):
            observed = {row[j] for row in seq_strs if row[j] != '-'}
            nongap = sum(row[j] != '-' for row in seq_strs) / n
            writer.writerow([j+1,rpos,aa,int(len(observed)>1),f'{nongap:.6f}'])
    neff = float(weights.sum())
    metrics = {'n_input': len(records), 'n_post_coverage':len(filtered), 'n_final':n,'L':L,
               'Neff':neff,'Neff_per_L':neff/L,'identity_threshold':args.identity_threshold,
               'weight_definition':'identical non-gap aa / all L model columns >= threshold',
               'min_row_coverage':args.min_row_coverage,'max_column_gap':args.max_column_gap,
               'expected_all_pairs':L*(L-1)//2,'anchors_checked':args.anchors,
               'gate_0_3':'PASS' if neff/L >= .3 else 'LOW_DEPTH_REVIEW'}
    (out / 'qc.json').write_text(json.dumps(metrics, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(metrics, indent=2, ensure_ascii=False))

if __name__ == '__main__':
    main()
