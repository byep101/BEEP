"""Run once on a reviewed, full-length curated FASTA, then enforce round-2 QC."""
from pathlib import Path
import argparse, csv, datetime, hashlib, importlib.util, json, math, os, subprocess, sys
import re

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True)
    parser.add_argument('--curation-note',required=True,help='Existing audit note documenting review of this FASTA')
    parser.add_argument('--analysis-id',default='TamI_round2')
    args=parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_]+',args.analysis_id): raise ValueError('analysis-id must be alphanumeric/underscore')
    r=Path.home()/'projects/tami_dca/TamI_round2'
    inp=Path(args.input).resolve(); note=Path(args.curation_note).resolve()
    if not note.is_file() or not note.read_text().strip(): raise ValueError('Curation audit note is missing')
    spec=importlib.util.spec_from_file_location('tami',r/'scripts/tami_minimal_v1.py')
    tami=importlib.util.module_from_spec(spec); spec.loader.exec_module(tami)
    rows=tami.read_fasta(inp)
    if not 2<=len(rows)<=1000: raise ValueError('Input count outside 2..1000; do not truncate')
    if any('-' in seq for sid,seq in rows): raise ValueError('Input must be unaligned full-length proteins')
    ref=tami.read_fasta(r/'00_baseline/TamI_reference.fasta')[0][1]
    if rows[0]!=('TamI_REF_ADC79647.1',ref): raise ValueError('Exact reference must be first')
    output=r/'03_msa'/f'{args.analysis_id}_full_MSA.faa'
    prepare=r/'04_prepare'/f'{args.analysis_id}_prepare_v1'
    if output.exists() or prepare.exists(): raise ValueError('Outputs already exist; no overwrite')
    prefix=r/'tools/prefix/usr'
    env=os.environ.copy(); env['MAFFT_BINARIES']=str(prefix/'lib/mafft/lib/mafft')
    strategy=['--localpair','--maxiterate','1000'] if len(rows)<200 else ['--auto']
    cmd=[str(prefix/'bin/mafft'),*strategy,'--thread','4',str(inp)]
    record={'argv':cmd,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_sha256':tami.digest(inp),'curation_note':str(note),'curation_note_sha256':tami.digest(note),'strategy_changed_from_v1':'--auto' in strategy,'baseline_full_length_restored':True,'functional_scope_changed_by_user':True}
    (r/'03_msa'/f'{args.analysis_id}_alignment_record.json').write_text(json.dumps(record,indent=2)+'\n')
    with output.open('w') as stdout,(r/'logs'/f'{args.analysis_id}_mafft.log').open('w') as stderr:
        subprocess.run(cmd,stdout=stdout,stderr=stderr,env=env,check=True)
    aligned=tami.read_fasta(output)
    assert len(aligned)==len(rows)
    assert len({len(s) for _,s in aligned})==1
    assert {i:s.replace('-','') for i,s in aligned}==dict(rows)
    assert dict(aligned)['TamI_REF_ADC79647.1'].replace('-','')==ref
    command=[sys.executable,str(r/'scripts/tami_minimal_v1.py'),'prepare','--msa',str(output),'--reference',str(r/'00_baseline/TamI_reference.fasta'),'--ref-id','TamI_REF_ADC79647.1','--out',str(prepare)]
    with (r/'logs'/f'{args.analysis_id}_prepare.log').open('w') as log: subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    q=json.loads((prepare/'input_qc.json').read_text())
    model=tami.read_fasta(prepare/'focus_cov70.faa')
    weights=[float(x) for x in (prepare/'weights.txt').read_text().split()]
    mapping=tami.read_table(prepare/'column_map.tsv')
    assert len(model)==len(weights)==q['model_rows']
    assert math.isclose(sum(weights),q['neff'],rel_tol=1e-12)
    assert all(0<w<=1 for w in weights)
    assert {len(s) for _,s in model}=={q['model_length']}
    assert q['baseline_reference_hash_match']
    for pos,residue in {101:'L',244:'L',295:'L',287:'E',290:'R',305:'E',359:'C'}.items():
        row=mapping[pos-1]
        assert row['reference_aa']==residue and row['status']=='retained',(pos,row)
        assert dict(model)['TamI_REF_ADC79647.1'][int(row['model_column'])-1]==residue
    assert min(float(x['column_occupancy_at_mask']) for x in mapping if x['status']=='retained')>=0.7
    focus=tami.read_fasta(prepare/'focus_all.faa')
    assert min(sum(c in tami.AA for c in seq)/413 for _,seq in focus)>=0.7
    motif=[]
    for sid,seq in focus:
        motif.append({'sequence_id':sid,'E287':seq[286],'R290':seq[289],'E305_manual_checkpoint':seq[304],'C359':seq[358],'ExxR_window':seq[286:290],'heme_window':seq[351:361],'three_anchor_match':seq[286]=='E' and seq[289]=='R' and seq[358]=='C'})
    tami.table(r/'05_qc'/f'{args.analysis_id}_motif_alignment_audit.tsv',list(motif[0]),motif)
    motif_fraction=sum(x['three_anchor_match'] for x in motif)/len(motif)
    baseline=json.loads((r/'00_baseline/TamI_v1_input_qc.json').read_text())
    control=json.loads((r/'05_qc/baseline_curated_prepare_control_v2/input_qc.json').read_text())
    compare=[{'dataset':'v1_frozen','N':baseline['N_model'],'L':baseline['L_model'],'Neff':baseline['Neff'],'Neff_over_L':baseline['Neff_over_L_model'],'scope':'historical broad P450 baseline'},
             {'dataset':'v1_curation_only_control','N':control['model_rows'],'L':control['model_length'],'Neff':control['neff'],'Neff_over_L':control['neff_over_model_length'],'scope':'45 explicit other-function annotations removed; inherited alignment'},
             {'dataset':'round2_curated','N':q['model_rows'],'L':q['model_length'],'Neff':q['neff'],'Neff_over_L':q['neff_over_model_length'],'scope':'reviewed provisional homologs; not confirmed TamI orthologs'}]
    tami.table(r/'05_qc'/f'{args.analysis_id}_v1_comparison.tsv',list(compare[0]),compare)
    stop=[]
    if q['neff_over_model_length']<=baseline['Neff_over_L_model']: stop.append('Neff/L did not improve versus frozen v1')
    if motif_fraction<0.7: stop.append('fewer than 70% have all three P450 anchors; inspect alignment/sequence membership')
    summary={'stage':'MSA and Neff completed','QC_pass':not stop,'stop_reasons':stop,'N_model':q['model_rows'],'L_model':q['model_length'],'Neff':q['neff'],'Neff_over_L':q['neff_over_model_length'],'planning_target_0.3_met':q['neff_over_model_length']>=0.3,'three_motif_anchor_fraction':motif_fraction,'motif_review_threshold':'0.7 operational review trigger, not a biological reliability threshold','scope_change_confounds_direct_comparison':True,'orthology_and_BGC_not_proven':True,'MIp_DCA_run':False}
    summary['analysis_id']=args.analysis_id
    summary['input_fasta']=str(inp)
    summary['input_sha256']=tami.digest(inp)
    (r/'05_qc'/f'{args.analysis_id}_validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    files=[prepare/'focus_cov70.faa',prepare/'weights.txt',prepare/'column_map.tsv']
    (r/'05_qc'/f'{args.analysis_id}_input_sha256.txt').write_text(''.join(tami.digest(p)+'  '+str(p.relative_to(r))+'\n' for p in files))
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
