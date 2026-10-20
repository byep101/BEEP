"""Apply the user's approved functional-scope decision as a versioned control."""
from pathlib import Path
import csv, json, importlib.util, subprocess, sys, datetime

r=Path.home()/'projects/tami_dca/TamI_round2'
spec=importlib.util.spec_from_file_location('tami',r/'scripts/tami_minimal_v1.py')
tami=importlib.util.module_from_spec(spec); spec.loader.exec_module(tami)
with (r/'05_qc/baseline_function_scope_audit.tsv').open(encoding='utf-8-sig') as f:
    audit=list(csv.DictReader(f,delimiter='\t'))
excluded={row['sequence_id'] for row in audit if row['scope_category']=='explicit_other_function_annotation'}
rows=tami.read_fasta(r/'00_baseline/TamI_v1_original_full_msa.faa')
kept=[(sid,seq) for sid,seq in rows if sid not in excluded]
dest=r/'02_intake/baseline_curated_v1'
dest.mkdir(exist_ok=False)
tami.write_fasta(dest/'baseline_curated_unaligned.faa',[(sid,seq.replace('-','')) for sid,seq in kept])
tami.write_fasta(dest/'baseline_curated_original_alignment.faa',kept)
for row in audit:
    row['curation_decision']='excluded_explicit_other_function' if row['sequence_id'] in excluded else 'provisional_homolog_not_orthology_validated'
    if row['scope_category']=='reference_TamI': row['curation_decision']='retain_reference'
tami.table(dest/'curation_audit.tsv',list(audit[0]),audit)
decision={'approved_scope':'TamI-focused; freeze v1 and exclude explicit other-function annotations from round2','user_decision_recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'N_original':len(rows),'N_excluded':len(excluded),'N_retained_including_reference':len(kept),'retained_nonreference_orthology_status':'unresolved; P450 annotation alone is not TamI functional evidence','purpose':'curation-only control; NOT expanded round2 dataset','excluded_ids':sorted(excluded)}
(dest/'curation_summary.json').write_text(json.dumps(decision,indent=2)+'\n')
out=r/'05_qc/baseline_curated_prepare_control'
subprocess.run([sys.executable,str(r/'scripts/tami_minimal_v1.py'),'prepare','--msa',str(dest/'baseline_curated_original_alignment.faa'),'--reference',str(r/'00_baseline/TamI_reference.fasta'),'--ref-id','TamI_REF_ADC79647.1','--out',str(out)],check=True)
summary=json.loads((r/'05_qc/baseline_audit_summary.json').read_text())
summary['scope_decision_pending']=False
summary['scope_decision']='TamI-focused curation approved by user'
(r/'05_qc/baseline_audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Excluded explicit other-function sequences:',len(excluded),'provisional retained:',len(kept),flush=True)
