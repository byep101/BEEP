#!/usr/bin/env python3
"""Map PDB modeled protein positions to Q8VS75 and score TxtC pair/ligand distances."""
import argparse
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBParser, is_aa
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

AA='ACDEFGHIKLMNPQRSTVWY'

def heavy_atoms(res):
    return [a for a in res.get_atoms() if getattr(a,'element','').upper() not in ('H','D')]

def min_dist(aa,bb):
    if not aa or not bb:return np.nan
    A=np.array([a.coord for a in aa]);B=np.array([a.coord for a in bb])
    return float(np.sqrt(((A[:,None,:]-B[None,:,:])**2).sum(axis=2)).min())

def representative(res):
    atom='CA' if res.resname.upper()=='GLY' else 'CB'
    if res.has_id(atom):return res[atom]
    return res['CA'] if res.has_id('CA') else None

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference',required=True);ap.add_argument('--position-map',required=True)
    ap.add_argument('--pairs',required=True,help='ranked MIp or DCA CSV with model_i/model_j/ref_i/ref_j')
    ap.add_argument('--structure',required=True,help='PDB .cif or .pdb')
    ap.add_argument('--chain',default='A');ap.add_argument('--out-prefix',required=True)
    ap.add_argument('--top-k',type=int,default=20)
    args=ap.parse_args()
    refrec=list(SeqIO.parse(args.reference,'fasta'))
    if len(refrec)!=1:raise SystemExit('Reference FASTA must have exactly one record')
    ref=str(refrec[0].seq).upper()
    parser=MMCIFParser(QUIET=True) if str(args.structure).lower().endswith('.cif') else PDBParser(QUIET=True)
    model=next(parser.get_structure('txtc',args.structure).get_models())
    if args.chain not in model:raise SystemExit(f'Chain {args.chain} not found; available: {[c.id for c in model]}')
    chain=model[args.chain]
    residues=[];pdb_chars=[]
    for r in chain:
        if not is_aa(r,standard=False):continue
        c=seq1(r.resname,custom_map={'MSE':'M','SEC':'C','PYL':'K'},undef_code='X')
        if c not in AA:continue
        residues.append(r);pdb_chars.append(c)
    seq=''.join(pdb_chars)
    aligner=PairwiseAligner(mode='global',match_score=2,mismatch_score=-4,
                            open_gap_score=-8,extend_gap_score=-0.4)
    match=aligner.align(ref,seq)[0]
    ref_to_res={};n_mapped=0;n_match=0
    for (a0,a1),(b0,b1) in zip(*match.aligned):
        if a1-a0!=b1-b0:raise SystemExit('Unexpected aligned block lengths')
        for a,b in zip(range(a0,a1),range(b0,b1)):
            n_mapped+=1
            if ref[a]==seq[b]:ref_to_res[a+1]=residues[b];n_match+=1
    identity=n_match/max(n_mapped,1)
    if identity<0.90:raise SystemExit(f'Mapped identity too low ({identity:.3f}); wrong protein / chain / indexing')
    pos=pd.read_csv(args.position_map)
    model_to_ref=dict(zip(pos.model_position.astype(int),pos.ref_position.astype(int)))
    pairs=pd.read_csv(args.pairs)
    if not {'model_i','model_j','ref_i','ref_j'}.issubset(pairs.columns):raise SystemExit('Pair CSV lacks model/ref fields')
    Path(args.out_prefix).parent.mkdir(parents=True,exist_ok=True)
    cands=[]
    for idx,r in enumerate(pairs.itertuples(index=False),start=1):
        mi,mj,ri,rj=int(r.model_i),int(r.model_j),int(r.ref_i),int(r.ref_j)
        if model_to_ref.get(mi)!=ri or model_to_ref.get(mj)!=rj:raise SystemExit('Model->reference pair mapping mismatch')
        res_i,res_j=ref_to_res.get(ri),ref_to_res.get(rj)
        if res_i is None or res_j is None:
            d=np.nan;b=np.nan;status='non_evaluable'
        else:
            d=min_dist(heavy_atoms(res_i),heavy_atoms(res_j))
            atom_i,atom_j=representative(res_i),representative(res_j)
            b=float(np.linalg.norm(atom_i.coord-atom_j.coord)) if atom_i is not None and atom_j is not None else np.nan
            status='evaluable' if np.isfinite(d) and np.isfinite(b) else 'non_evaluable'
        cands.append({'rank_in_input':idx,'ref_i':ri,'ref_j':rj,'model_i':mi,'model_j':mj,
                      'min_heavy_A':d,'Cb_Cb_A':b,'status':status,
                      'heavy_contact_le5':(d<=5) if status=='evaluable' else '',
                      'Cb_contact_le8':(b<=8) if status=='evaluable' else ''})
    scored=pd.DataFrame(cands)
    scored.to_csv(args.out_prefix+'_contacts.csv',index=False)
    mapping=[]
    for rpos,res in sorted(ref_to_res.items()):
        icode=res.id[2].strip()
        mapping.append({'ref_position':rpos,'reference_aa':ref[rpos-1],
                        'pdb_chain':args.chain,'pdb_author_resid':res.id[1],
                        'insertion_code':icode,'pdb_aa':seq1(res.resname,undef_code='X')})
    pd.DataFrame(mapping).to_csv(args.out_prefix+'_pdb_mapping.csv',index=False)
    hemes=[r for r in model.get_residues() if r.resname.upper() in ('HEM','HEB')]
    substrate=[r for r in model.get_residues() if r.resname.upper() in ('C8H','C8B')]
    heme_atoms=sum([heavy_atoms(r) for r in hemes],[])
    substrate_atoms=sum([heavy_atoms(r) for r in substrate],[])
    contact_res=[]
    for ri in sorted(set(scored.ref_i)|set(scored.ref_j)):
        res=ref_to_res.get(int(ri))
        contact_res.append({'ref_position':int(ri),'ref_aa':ref[int(ri)-1],
                            'mapped':int(res is not None),
                            'to_HEM_min_heavy_A':min_dist(heavy_atoms(res),heme_atoms) if res else np.nan,
                            'to_substrate_min_heavy_A':min_dist(heavy_atoms(res),substrate_atoms) if res else np.nan})
    pd.DataFrame(contact_res).to_csv(args.out_prefix+'_ligand_distances.csv',index=False)
    # Find candidate heme-coordinating cysteine from 3D by nearest SG to Fe, then read its mapped reference index.
    fe_atoms=[a for r in hemes for a in r if a.name.upper()=='FE']
    axial=[]
    if fe_atoms:
        for refi,r in ref_to_res.items():
            if ref[refi-1]=='C' and r.has_id('SG'):
                axial.append((float(min_dist([r['SG']],fe_atoms)),refi,r.id[1]))
        axial.sort()
    ev=scored[scored.status=='evaluable']
    top=scored.head(args.top_k);te=top[top.status=='evaluable']
    lines=[f'Structure: {args.structure}, chain: {args.chain}',
           f'Mapped identical residues: {n_match}/{n_mapped} ({identity:.2%})',
           f'Ligands found: HEM={len(hemes)}, substrate C8H/C8B={len(substrate)}',
           f'Top{args.top_k} evaluable pairs: {len(te)}/{len(top)}',
           f'Top{args.top_k} min-heavy <=5A: {sum(te.min_heavy_A<=5)}/{len(te)} (only evaluable denominator)',
           f'Top{args.top_k} Cb-Cb <=8A: {sum(te.Cb_Cb_A<=8)}/{len(te)} (only evaluable denominator)',
           f'Background min-heavy <=5A: {sum(ev.min_heavy_A<=5)}/{len(ev)}',
           f'Background Cb-Cb <=8A: {sum(ev.Cb_Cb_A<=8)}/{len(ev)}',
           f'Nearest Cys-SG to heme Fe: {axial[0] if axial else "not found"} (distance A, reference pos, PDB author pos)']
    Path(args.out_prefix+'_summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('\n'.join(lines))

if __name__=='__main__':main()
