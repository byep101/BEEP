#!/usr/bin/env python3
"""Local, standard-library-only TamI pilot: prepare -> mi -> compare.

This is a workflow pilot, not a validated functional-residue predictor.
MIp: weighted empirical MI with pairwise missing-data deletion, followed by
Dunn-style off-diagonal APC. The weighting definition is explicitly recorded.
No sequence search, network upload, tree inference or automatic master edit.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys

VERSION = '1.0.0'
AA = set('ACDEFGHIKLMNPQRSTVWY')
MISSING = set('XBZJUO?*')
BENCHMARK = {101: 'L', 244: 'L', 295: 'L'}
REF_HASH = 'fbc563a7aeebba572bac52607531eb227351b690037bf8108f176a3f8701027a'


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def save_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')


def table(path: Path, fields, rows, delimiter='\t') -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8-sig', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter=delimiter, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def read_table(path: Path, delimiter='\t'):
    with path.open(encoding='utf-8-sig', newline='') as inp:
        return list(csv.DictReader(inp, delimiter=delimiter))


def read_fasta(path: Path):
    rows = []
    header = None
    seq = []
    with path.open(encoding='utf-8-sig') as inp:
        for line in inp:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    rows.append((header.split()[0], ''.join(seq)))
                header = line[1:].strip()
                if not header:
                    raise ValueError('Empty FASTA header')
                seq = []
            else:
                if header is None:
                    raise ValueError('Sequence before first FASTA header')
                seq.append(''.join(line.split()).upper().replace('.', '-'))
    if header is not None:
        rows.append((header.split()[0], ''.join(seq)))
    if not rows or any(not s for _, s in rows):
        raise ValueError('Empty FASTA or empty sequence')
    ids = [x[0] for x in rows]
    if len(set(ids)) != len(ids):
        raise ValueError('FASTA IDs must be unique')
    return rows


def write_fasta(path: Path, rows) -> None:
    with path.open('w', encoding='ascii', newline='\n') as out:
        for ident, seq in rows:
            out.write('>' + ident + '\n')
            for start in range(0, len(seq), 80):
                out.write(seq[start:start+80] + '\n')


def check_run(run: Path):
    q = load_json(run / 'input_qc.json')
    for name, expected in q['output_hashes'].items():
        if digest(run / name) != expected:
            raise ValueError(f'{name} changed since prepare; create a new run directory')
    return q


def sequence_weights(seqs, theta=0.8, minimum_overlap=0.7):
    """Symmetric, nongap-pair identity; requires shared coverage of model length.

    A neighbor has identity >= theta over shared 20-AA sites and at least
    minimum_overlap * L shared sites. Each sequence counts itself once.
    These are project-defined external weights, not a claim about EV defaults.
    """
    length = len(seqs[0])
    counts = [1] * len(seqs)
    for a in range(len(seqs)):
        for b in range(a):
            shared = [(x, y) for x, y in zip(seqs[a], seqs[b]) if x in AA and y in AA]
            if len(shared) / length >= minimum_overlap:
                ident = sum(x == y for x, y in shared) / len(shared)
                if ident >= theta:
                    counts[a] += 1
                    counts[b] += 1
    return [1.0 / n for n in counts], counts


def prepare(args):
    src, ref_file, out = Path(args.msa), Path(args.reference), Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output directory is not empty; use a new run name (no overwrite)')
    rows = read_fasta(src)
    if len(rows) > 1000:
        raise ValueError('This pilot is limited to <=1000 rows. Check that the input is TamI, not TxtC expansion.')
    if len({len(s) for _, s in rows}) != 1:
        raise ValueError('Input is not a rectangular FASTA MSA; do not use raw FASTA/A3M')
    for ident, seq in rows:
        bad = set(seq) - AA - MISSING - {'-'}
        if bad:
            raise ValueError(f'Unexpected symbols in {ident}: {sorted(bad)}')
    ref_rows = read_fasta(ref_file)
    if len(ref_rows) != 1:
        raise ValueError('Reference FASTA must contain exactly one sequence')
    original_ref = ref_rows[0][1]
    if any(x not in AA for x in original_ref):
        raise ValueError('Reference must be ungapped, standard 20-AA sequence')
    if len(original_ref) != 413:
        raise ValueError('Expected the 413-aa ADC79647.1/D3Y1J3 TamI reference')
    matches = [s for ident, s in rows if ident == args.ref_id]
    if len(matches) != 1:
        raise ValueError(f'Expected one exact reference ID {args.ref_id}; inspect FASTA headers')
    aligned_ref = matches[0]
    if aligned_ref.replace('-', '') != original_ref:
        raise ValueError('Ungapped reference row does not equal reference FASTA; do not force numbering')
    if hashlib.sha256(original_ref.encode('ascii')).hexdigest() != REF_HASH and not args.allow_reference_hash_difference:
        raise ValueError('Reference differs from the recorded TamI v0.1 sequence SHA256. Verify the source; only then use --allow-reference-hash-difference.')
    for pos, aa in BENCHMARK.items():
        if original_ref[pos-1] != aa:
            raise ValueError(f'Reference identity mismatch: expected {aa}{pos}')
    full_cols = [i for i, aa in enumerate(aligned_ref) if aa != '-']
    rows = [(args.ref_id, aligned_ref)] + [(ident, seq) for ident, seq in rows if ident != args.ref_id]
    clean, audit, seen = [], [], set()
    for ident, seq in rows:
        changes = sum(c in MISSING for c in seq)
        seq = ''.join(c if c in AA or c == '-' else '-' for c in seq)
        focus = ''.join(seq[c] for c in full_cols)
        cov = sum(c in AA for c in focus) / len(focus)
        reason = 'retained'
        if cov < 0.7:
            reason = 'excluded_focus_coverage_lt_0.7'
        elif focus in seen:
            reason = 'excluded_exact_focus_duplicate'
        if reason == 'retained':
            clean.append((ident, focus))
            seen.add(focus)
        audit.append({'sequence_id': ident, 'focus_coverage': cov,
                      'ambiguous_to_gap_in_full_msa': changes, 'decision': reason})
    if len(clean) < 2:
        raise ValueError('Fewer than two distinct sequences remain; do not fit a coevolution model')
    occupancy = [sum(s[j] in AA for _, s in clean) / len(clean) for j in range(413)]
    keep = [j for j, x in enumerate(occupancy) if x >= 0.7]
    if len(keep) < 2:
        raise ValueError('Too few occupied columns remain')
    compact = [(ident, ''.join(s[j] for j in keep)) for ident, s in clean]
    # Identical rows created by masking are retained and downweighted, so the
    # previously computed column-occupancy mask remains unchanged.
    model_rows = compact
    if len({s for _, s in model_rows}) < 2:
        raise ValueError('No sequence variation remains after masking')
    seqs = [s for _, s in model_rows]
    weights, neighbors = sequence_weights(seqs)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, out / 'source_full_msa.faa')
    shutil.copyfile(ref_file, out / 'source_reference.faa')
    write_fasta(out / 'focus_all.faa', clean)
    write_fasta(out / 'focus_cov70.faa', model_rows)
    model_map = {c: i+1 for i, c in enumerate(keep)}
    mapping = [{'reference_position': j+1, 'reference_aa': original_ref[j],
                'full_msa_column': full_cols[j]+1, 'column_occupancy_at_mask': occupancy[j],
                'model_column': model_map.get(j, ''),
                'status': 'retained' if j in model_map else 'filtered_low_occupancy'} for j in range(413)]
    table(out / 'column_map.tsv', list(mapping[0]), mapping)
    table(out / 'sequence_audit.tsv', list(audit[0]), audit)
    with (out / 'weights.txt').open('w', encoding='ascii', newline='\n') as f:
        f.write(''.join(f'{x:.17g}\n' for x in weights))
    wrows = [{'sequence_id': r[0], 'weight': w, 'neighbors_including_self': n}
             for r, w, n in zip(model_rows, weights, neighbors)]
    table(out / 'sequence_weights.tsv', list(wrows[0]), wrows)
    q = {'version': VERSION, 'time_utc': datetime.now(timezone.utc).isoformat(),
         'reference_id': args.ref_id, 'input_rows': len(rows), 'full_msa_columns': len(aligned_ref),
         'focus_rows_before_model_dedup': len(clean), 'model_rows': len(model_rows),
         'reference_length': 413, 'model_length': len(keep), 'neff': sum(weights),
         'neff_over_model_length': sum(weights)/len(keep), 'neff_over_reference_length': sum(weights)/413,
         'weighting': {'theta_identity': 0.8, 'minimum_shared_fraction_of_model': 0.7,
                      'identity_denominator': 'pairwise jointly observed 20-AA sites',
                      'self_included': True, 'method': 'project_defined_external_weights'},
         'coverage_cutoffs': {'sequence_focus': 0.7, 'column_at_mask_before_model_dedup': 0.7},
         'model_min_column_occupancy_after_dedup': min(sum(s[j] in AA for s in seqs)/len(seqs) for j in range(len(keep))),
         'reference_sequence_sha256': hashlib.sha256(original_ref.encode('ascii')).hexdigest(),
         'baseline_reference_hash_match': hashlib.sha256(original_ref.encode('ascii')).hexdigest() == REF_HASH,
         'source_msa_sha256': digest(src), 'source_reference_file_sha256': digest(ref_file),
         'benchmark': {str(k): mapping[k-1] for k in BENCHMARK},
         'interpretation': 'method-development pilot; no depth threshold certifies biological validity',
         'output_hashes': {name: digest(out/name) for name in ['focus_cov70.faa','column_map.tsv','weights.txt']}}
    save_json(out / 'input_qc.json', q)
    print(json.dumps({k:q[k] for k in ['input_rows','model_rows','model_length','neff','neff_over_model_length']}, indent=2))


def weighted_mi(columns, weights):
    length = len(columns)
    mi = [[0.0]*length for _ in range(length)]
    support = [[0.0]*length for _ in range(length)]
    entropy = []
    for col in columns:
        c = Counter()
        for aa, w in zip(col, weights):
            if aa in AA: c[aa] += w
        total = sum(c.values())
        entropy.append(-sum((v/total)*math.log(v/total) for v in c.values()) if total else 0.0)
    for i in range(length):
        for j in range(i+1, length):
            counts = Counter()
            for a, b, w in zip(columns[i], columns[j], weights):
                if a in AA and b in AA: counts[(a,b)] += w
            total = sum(counts.values())
            if total <= 0:
                raise ValueError(f'No joint observations at model pair {i+1},{j+1}')
            fa, fb = Counter(), Counter()
            for (a,b), v in counts.items(): fa[a] += v; fb[b] += v
            score = sum((v/total)*math.log((v*total)/(fa[a]*fb[b])) for (a,b),v in counts.items())
            mi[i][j] = mi[j][i] = max(0.0, score)  # numerical negative round-off only
            support[i][j] = support[j][i] = total
    mean = [sum(row)/(length-1) for row in mi]
    grand = sum(mean)/length
    mip = [[0.0]*length for _ in range(length)]
    for i in range(length):
        for j in range(i+1,length):
            mip[i][j] = mip[j][i] = mi[i][j] - (mean[i]*mean[j]/grand if grand > 1e-15 else 0.0)
    return mi, mip, support, entropy


def model_context(run):
    q = check_run(run)
    rows = read_fasta(run / 'focus_cov70.faa')
    maps = sorted([r for r in read_table(run/'column_map.tsv') if r['model_column']], key=lambda r:int(r['model_column']))
    weights = [float(x) for x in (run/'weights.txt').read_text().splitlines() if x.strip()]
    if len(weights) != len(rows) or len(maps) != len(rows[0][1]):
        raise ValueError('Alignment/weight/mapping length mismatch')
    return q, rows, maps, weights


def run_mi(args):
    run = Path(args.run)
    if (run/'mi_summary.json').exists(): raise ValueError('MI outputs already exist; use a new run directory to rerun')
    q, rows, maps, weights = model_context(run)
    columns = list(zip(*(s for _,s in rows)))
    mi, mip, support, entropy = weighted_mi(columns, weights)
    scores = []
    for i in range(len(maps)):
        for j in range(i+1,len(maps)):
            rp_i, rp_j = int(maps[i]['reference_position']), int(maps[j]['reference_position'])
            eligible = rp_j-rp_i >= 6 and entropy[i] > 1e-12 and entropy[j] > 1e-12
            scores.append({'model_i': i+1, 'model_j': j+1, 'reference_i':rp_i, 'reference_j':rp_j,
                           'aa_i':maps[i]['reference_aa'], 'aa_j':maps[j]['reference_aa'],
                           'mi':mi[i][j], 'mip':mip[i][j], 'pair_weight_sum':support[i][j],
                           'eligible_nonlocal_variable':int(eligible)})
    table(run/'MIp_pairs.csv', list(scores[0]), scores, ',')
    site_rows = [{'model_column': i+1, 'reference_position': int(r['reference_position']),
                  'entropy_nats': entropy[i], 'variable': int(entropy[i] > 1e-12)} for i,r in enumerate(maps)]
    table(run/'site_entropy.tsv', list(site_rows[0]), site_rows)
    save_json(run/'mi_summary.json', {'status':'computed_not_validated', 'pair_count':len(scores),
              'eligible_pairs':sum(r['eligible_nonlocal_variable'] for r in scores),
              'frequency_estimator':'weighted empirical, no pseudocount, pairwise delete missing',
              'log_base':'e', 'apc':'off-diagonal row means, diagonal excluded; compute before sequence-distance filtering',
              'min_reference_separation':6, 'alignment_sha256':digest(run/'focus_cov70.faa'),
              'weights_sha256':digest(run/'weights.txt'), 'pairs_sha256':digest(run/'MIp_pairs.csv'),
              'site_entropy_sha256':digest(run/'site_entropy.tsv')})
    print('MI/APC complete. Numerical scores are not functional validation.')


def rank_min(values):
    # Equal scores receive identical minimum rank, avoiding ID-based tie claims.
    sorted_vals = sorted(values.items(), key=lambda x:(-x[1],x[0]))
    ranks, last, rank = {}, None, 0
    for i,(key,value) in enumerate(sorted_vals,1):
        if last is None or value != last: rank = i
        ranks[key] = rank; last = value
    return ranks


def compare(args):
    run = Path(args.run)
    q = check_run(run)
    ms = load_json(run/'mi_summary.json')
    if digest(run/'MIp_pairs.csv') != ms['pairs_sha256'] or digest(run/'site_entropy.tsv') != ms['site_entropy_sha256']:
        raise ValueError('MI outputs have changed')
    mirows = read_table(run/'MIp_pairs.csv',',')
    dca_path = run/args.dca_label/'DCA_pairs.csv'
    if not dca_path.exists():
        raise ValueError('DCA_pairs.csv is missing. MI alone does not complete the combined pilot.')
    ds = load_json(run/args.dca_label/'dca_summary.json')
    if not ds.get('outputs_validated'): raise ValueError('DCA outputs did not pass structural file checks')
    for k in ['alignment_sha256','weights_sha256']:
        if ms[k] != ds[k]: raise ValueError('MI and DCA do not use the same alignment/weights')
    if digest(dca_path) != ds['pairs_sha256']: raise ValueError('DCA output has changed')
    d = {(int(r['reference_i']),int(r['reference_j'])):float(r['cn']) for r in read_table(dca_path,',')}
    all_pairs = {(int(r['reference_i']),int(r['reference_j'])) for r in mirows}
    if set(d) != all_pairs: raise ValueError('MI/DCA pair coordinate universes differ')
    eligible = [(int(r['reference_i']),int(r['reference_j'])) for r in mirows if r['eligible_nonlocal_variable']=='1']
    mip = {(int(r['reference_i']),int(r['reference_j'])):float(r['mip']) for r in mirows}
    eligible_set = set(eligible)
    methods = {'MIp':{p:s for p,s in mip.items() if p in eligible_set}, 'DCA_CN':{p:s for p,s in d.items() if p in eligible_set}}
    site_scores, site_ranks = {},{}
    for method,pairs in methods.items():
        best = {}
        for (i,j),v in pairs.items():
            for pos in [i,j]: best[pos] = max(best.get(pos,-math.inf),v)
        site_scores[method]=best; site_ranks[method]=rank_min(best)
        ranked = sorted(pairs,key=lambda p:(-pairs[p],p))
        table(run/f'{method}_top20_pairs.csv',['reference_i','reference_j','score'],
              [{'reference_i':p[0],'reference_j':p[1],'score':pairs[p]} for p in ranked[:20]],',')
    retained = {int(r['reference_position']) for r in read_table(run/'column_map.tsv') if r['model_column']}
    model_map = {int(r['reference_position']):r['model_column'] for r in read_table(run/'column_map.tsv')}
    report = []
    for pos in BENCHMARK:
        status = 'filtered_low_occupancy' if pos not in retained else ('scored' if pos in site_ranks['MIp'] else 'unscorable_constant_or_no_eligible_partner')
        report.append({'residue':f'L{pos}','reference_position':pos,'model_column':model_map[pos],
                       'status':status,'MIp_rank':site_ranks['MIp'].get(pos,''),
                       'DCA_CN_rank':site_ranks['DCA_CN'].get(pos,''),
                       'MIp_top20':int(site_ranks['MIp'].get(pos,10**9)<=20),
                       'DCA_top20':int(site_ranks['DCA_CN'].get(pos,10**9)<=20)})
    table(run/'benchmark_summary.csv',list(report[0]),report,',')
    ranks_rows=[]
    for pos in sorted(set(site_scores['MIp']) | set(site_scores['DCA_CN'])):
        ranks_rows.append({'reference_position':pos,'MIp_rank':site_ranks['MIp'].get(pos,''),
                           'DCA_CN_rank':site_ranks['DCA_CN'].get(pos,''),
                           'MIp_max_pair_score':site_scores['MIp'].get(pos,''),
                           'DCA_max_pair_score':site_scores['DCA_CN'].get(pos,'')})
    table(run/'residue_ranks.csv',['reference_position','MIp_rank','DCA_CN_rank','MIp_max_pair_score','DCA_max_pair_score'],ranks_rows,',')
    comparisons=[]
    for k in [20,50]:
        actual=min(k,len(eligible))
        sets=[set(sorted(p,key=lambda x:(-p[x],x))[:actual]) for p in methods.values()]
        inter=len(sets[0]&sets[1]); union=len(sets[0]|sets[1])
        comparisons.append({'requested_k':k,'actual_k':actual,'overlap_pairs':inter,
                            'jaccard':inter/union if union else None,'interpretation':'descriptive only; ties cut by coordinate'})
    status='scores_generated_optimization_review' if not ds.get('optimizer_converged') else 'technical_pilot_complete'
    summary={'status':status,'N':q['model_rows'],'L_model':q['model_length'],'Neff':q['neff'],
             'Neff_L':q['neff_over_model_length'],'comparisons':comparisons,
             'residue_score_definition':'maximum eligible pair score for each site; descriptive heuristic',
             'benchmark':report,'biological_conclusion':'not established','dca_optimization_status':ds['optimization_status']}
    save_json(run/'comparison_summary.json',summary)
    text=f"""# TamI MIp/APC + DCA pilot record\n\nStatus: {status}\n\nN = {q['model_rows']}; model L = {q['model_length']}; Neff = {q['neff']:.3f}; Neff/L = {q['neff_over_model_length']:.4f}.\n\nThe two methods use the same focus alignment and external sequence weights.\nDCA optimizer: {ds['optimization_status']}\n\nSee benchmark_summary.csv and residue_ranks.csv. A rank or a recovered known\nsite is not evidence of catalytic causality, nor a performance estimate over all\nfunctional residues. Filtered and invariant sites must not be called true negatives.\n\nPending: biological homolog/BGC-context review; stability/null controls;\nstructure annotation and stronger depth evaluation. Other systems remain paused.\n"""
    (run/'run_record.md').write_text(text,encoding='utf-8')
    print(json.dumps(summary,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    subs=p.add_subparsers(dest='action',required=True)
    a=subs.add_parser('prepare'); a.add_argument('--msa',required=True); a.add_argument('--reference',required=True)
    a.add_argument('--ref-id',default='TamI_REF_ADC79647.1'); a.add_argument('--out',required=True)
    a.add_argument('--allow-reference-hash-difference',action='store_true',help='Only after manually verifying a genuinely updated reference')
    a.set_defaults(func=prepare)
    a=subs.add_parser('mi'); a.add_argument('--run',required=True); a.set_defaults(func=run_mi)
    a=subs.add_parser('compare'); a.add_argument('--run',required=True); a.add_argument('--dca-label',default='dca'); a.set_defaults(func=compare)
    args=p.parse_args()
    try: args.func(args)
    except (ValueError, OSError, KeyError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr); return 1
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
