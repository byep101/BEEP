#!/usr/bin/env python3
"""Weighted 21-state MI with off-diagonal APC; uses SAME external weights as plmc."""
import argparse
import csv
from pathlib import Path
import numpy as np
from Bio import SeqIO

ABC='-ACDEFGHIKLMNPQRSTVWY'

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--msa',required=True);ap.add_argument('--weights',required=True)
    ap.add_argument('--position-map',required=True);ap.add_argument('--output',required=True)
    args=ap.parse_args()
    seqs=[str(x.seq).upper() for x in SeqIO.parse(args.msa,'fasta')]
    if not seqs or len({len(s) for s in seqs})!=1: raise SystemExit('Invalid MSA length')
    L=len(seqs[0]);N=len(seqs)
    w=np.atleast_1d(np.loadtxt(args.weights));
    if len(w)!=N or np.any(w <=0) or not np.isfinite(w).all(): raise SystemExit('Weights N/validity mismatch')
    w=w/w.sum(); code=np.array([[ABC.index(c) for c in seq] for seq in seqs],dtype=np.int16)
    with open(args.position_map,encoding='utf-8',newline='') as f:
        pm={int(r['model_position']):int(r['ref_position']) for r in csv.DictReader(f)}
    if len(pm)!=L:raise SystemExit('Position map length mismatch')
    one=np.array([np.bincount(code[:,j],weights=w,minlength=21) for j in range(L)])
    mi=np.zeros((L,L),dtype=np.float64)
    for i in range(L):
        for j in range(i+1,L):
            joint=np.bincount(code[:,i]*21+code[:,j],weights=w,minlength=441).reshape((21,21))
            independent=np.outer(one[i],one[j])
            nz=joint>0
            mi[i,j]=mi[j,i]=np.sum(joint[nz]*np.log(joint[nz]/independent[nz]))
    mean_i=mi.sum(axis=1)/(L-1)
    global_mean=mean_i.mean()
    apc=np.outer(mean_i,mean_i)/global_mean if global_mean>0 else np.zeros_like(mi)
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    with open(args.output,'w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        writer.writerow(['model_i','model_j','ref_i','ref_j','MI_raw','APC','MIp'])
        for i in range(L):
            for j in range(i+1,L):
                writer.writerow([i+1,j+1,pm[i+1],pm[j+1],f'{mi[i,j]:.10g}',f'{apc[i,j]:.10g}',f'{mi[i,j]-apc[i,j]:.10g}'])
    print(f'Wrote MIp score for {L*(L-1)//2} pairs; N={N}, L={L}, Neff={np.sum(np.loadtxt(args.weights)):.5f}')

if __name__=='__main__':main()
