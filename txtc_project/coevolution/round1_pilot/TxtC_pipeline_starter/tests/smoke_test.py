#!/usr/bin/env python3
"""Synthetic input smoke test only. No biological TxtC claims."""
from pathlib import Path
import csv
import tempfile
import subprocess
import sys
import random
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=ROOT/'scripts'
AA='ACDEFGHIKLMNPQRSTVWY'

def run(*args):
    cmd=[sys.executable,str(SCRIPTS/args[0])]+list(map(str,args[1:]))
    p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    print(p.stdout)
    if p.returncode:raise RuntimeError('Failed: '+repr(cmd))

with tempfile.TemporaryDirectory() as td:
    d=Path(td)
    rnd=random.Random(5)
    ref='ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQ' # 34 positions
    SeqIO.write([SeqRecord(Seq(ref),id='Q8VS75')],d/'ref.fasta','fasta')
    cands=[]
    for n in range(18):
        a=list(ref)
        for j in range(len(a)):
            if rnd.random()<.22:
                a[j]=rnd.choice(AA)
        cands.append(SeqRecord(Seq(''.join(a)),id=f'test{n}'))
    cands.append(SeqRecord(Seq(ref),id='duplicate_ref'))
    SeqIO.write(cands,d/'cand.fasta','fasta')
    run('01_prepare_fasta.py','--reference',d/'ref.fasta','--candidates',d/'cand.fasta',
        '--output',d/'clean.fasta','--manifest',d/'manifest.csv')
    # all inputs have same length for synthetic MSA
    (d/'msa.fasta').write_text((d/'clean.fasta').read_text())
    run('02_focus_weights.py','--msa',d/'msa.fasta','--outdir',d/'focus',
        '--anchors','1:A')
    run('03_mip.py','--msa',d/'focus/txtc_focus.fasta','--weights',d/'focus/weights.txt',
        '--position-map',d/'focus/position_map.csv','--output',d/'mip.csv')
    with open(d/'cn.txt','w') as f:
        for i in range(len(ref)):
            for j in range(i+1,len(ref)):
                f.write(f'{i+1} {ref[i]} {j+1} {ref[j]} 0 {rnd.random():.8f}\n')
    run('04_compare.py','--mip',d/'mip.csv','--plmc',d/'cn.txt',
        '--position-map',d/'focus/position_map.csv','--site-qc',d/'focus/site_qc.csv',
        '--outdir',d/'compare')
    # Minimal synthetic PDB with synthetic coordinates, not any real protein.
    three={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE',
           'K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER',
           'T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
    lines=[];atomnum=0
    for i,c in enumerate(ref,1):
        for name,dx,el in [('CA',0.,'C'),('CB',1.2,'C')]:
            if c=='G' and name=='CB':continue
            atomnum+=1
            x=2.0*i+dx
            lines.append(f'ATOM  {atomnum:5d} {name:4s} {three[c]:3s} A{i:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}          {el:>2s}\n')
    lines+=['TER\n','END\n']
    (d/'fake.pdb').write_text(''.join(lines))
    run('05_structure_contacts.py','--reference',d/'ref.fasta','--position-map',d/'focus/position_map.csv',
        '--pairs',d/'compare/dca_ranked.csv','--structure',d/'fake.pdb','--chain','A',
        '--out-prefix',d/'structure/fake')
    assert (d/'structure/fake_contacts.csv').exists()
    assert sum(1 for _ in open(d/'mip.csv'))==len(ref)*(len(ref)-1)//2+1
    print('SMOKE TEST PASS on synthetic data; no real biological dataset was used.')
