#!/usr/bin/env python3
"""Strictly validate plmc six-column couplings and compare with MIp on the same pair set."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mip',required=True);ap.add_argument('--plmc',required=True)
    ap.add_argument('--position-map',required=True);ap.add_argument('--site-qc',required=True)
    ap.add_argument('--outdir',required=True);ap.add_argument('--min-separation',type=int,default=6)
    args=ap.parse_args();out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    mp=pd.read_csv(args.position_map);L=len(mp)
    pm=dict(zip(mp.model_position.astype(int),mp.ref_position.astype(int)))
    mip=pd.read_csv(args.mip)
    required={'model_i','model_j','MIp','ref_i','ref_j'}
    if not required.issubset(mip.columns):raise SystemExit('Missing MIp columns')
    dca=[]
    with open(args.plmc,encoding='utf-8') as f:
        for line_no,line in enumerate(f,1):
            line=line.strip()
            if not line or line.startswith('#'):continue
            p=line.split()
            if len(p)!=6:raise SystemExit(f'plmc line {line_no}: expected 6 columns; got {len(p)}')
            try:i,j,score=int(p[0]),int(p[2]),float(p[5])
            except Exception as e:raise SystemExit(f'plmc invalid line {line_no}: {e}')
            if i>j:i,j=j,i
            if not(1<=i<j<=L):raise SystemExit(f'plmc residue indices {i},{j} outside 1..L; STOP: mapping mismatch')
            dca.append((i,j,score))
    dc=pd.DataFrame(dca,columns=['model_i','model_j','CN'])
    exp=L*(L-1)//2
    if len(dc)!=exp or len(mip)!=exp:
        raise SystemExit(f'Not all scores written: expected {exp}, MIp={len(mip)}, DCA={len(dc)}')
    if dc.duplicated(['model_i','model_j']).any():raise SystemExit('Duplicated plmc pair')
    merged=mip.merge(dc,on=['model_i','model_j'],validate='one_to_one')
    if len(merged)!=exp:raise SystemExit('plmc and MIp pair spaces do not match exactly')
    if any(int(row.ref_i)!=pm[int(row.model_i)] or int(row.ref_j)!=pm[int(row.model_j)] for row in merged.itertuples()):
        raise SystemExit('MIp reference mapping mismatch')
    site=pd.read_csv(args.site_qc)
    variable={int(r.model_position) for r in site.itertuples() if int(r.is_variable)==1}
    eligible=merged[(abs(merged.ref_i-merged.ref_j)>=args.min_separation) &
                    merged.model_i.isin(variable) & merged.model_j.isin(variable)].copy()
    if eligible.empty:raise SystemExit('No eligible nonlocal variable pairs')
    a=eligible.sort_values('MIp',ascending=False,kind='stable').reset_index(drop=True)
    b=eligible.sort_values('CN',ascending=False,kind='stable').reset_index(drop=True)
    a['MIp_rank']=np.arange(1,len(a)+1);b['DCA_rank']=np.arange(1,len(b)+1)
    both=a[['model_i','model_j','ref_i','ref_j','MIp','MIp_rank']].merge(
        b[['model_i','model_j','CN','DCA_rank']],on=['model_i','model_j'],validate='one_to_one')
    both.to_csv(out/'ranked_comparison.csv',index=False)
    a.to_csv(out/'mip_ranked.csv',index=False)
    b.to_csv(out/'dca_ranked.csv',index=False)
    scores=[];rho,pvalue=spearmanr(both.MIp,both.CN)
    scores.append(f'All pair scores: {exp}')
    scores.append(f'Eligible nonlocal and variable pairs: {len(eligible)}')
    scores.append(f'Spearman rho on eligible pairs: {rho:.6f}; p={pvalue:.4g}')
    for k in [20,50,100]:
        top_a=set(zip(a.head(k).model_i,a.head(k).model_j))
        top_b=set(zip(b.head(k).model_i,b.head(k).model_j))
        common=top_a&top_b
        both[both[['model_i','model_j']].apply(tuple,axis=1).isin(common)].sort_values('MIp_rank').to_csv(out/f'intersection_top{k}.csv',index=False)
        scores.append(f'Top{k} intersection: {len(common)} / {k} ({len(common)/k:.1%})')
        a.head(k).to_csv(out/f'mip_top{k}.csv',index=False)
        b.head(k).to_csv(out/f'dca_top{k}.csv',index=False)
    # Max-scoring eligible pair is a transparent residue-level benchmark, NOT previous TamI rank definition.
    residues=[]
    for mpos in range(1,L+1):
        aa=a[(a.model_i==mpos)|(a.model_j==mpos)]
        bb=b[(b.model_i==mpos)|(b.model_j==mpos)]
        residues.append({'model_position':mpos,'ref_position':pm[mpos],
                         'MIp_best_pair_rank':int(aa.iloc[0].MIp_rank) if len(aa) else '',
                         'DCA_best_pair_rank':int(bb.iloc[0].DCA_rank) if len(bb) else ''})
    pd.DataFrame(residues).to_csv(out/'residue_best_pair_rank.csv',index=False)
    (out/'comparison_summary.txt').write_text('\n'.join(scores)+'\n',encoding='utf-8')
    print('\n'.join(scores))

if __name__=='__main__':main()
