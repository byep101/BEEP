from pathlib import Path
import csv,json

r=Path.home()/'projects/tami_dca/TamI_round2'
from collect_and_filter import read_fasta,table
rows=[]
for label,path in [('v1_frozen',r/'00_baseline/TamI_v1_focus_all.faa'),('v1_curated_control',r/'05_qc/baseline_curated_prepare_control/focus_all.faa')]:
    for sid,seq in read_fasta(path):
        rows.append({'dataset':label,'sequence_id':sid,'E287':seq[286],'R290':seq[289],'E305_manual_checkpoint':seq[304],'C359':seq[358],'anchor_match':seq[286]=='E' and seq[289]=='R' and seq[358]=='C','ExxR_window':seq[286:290],'heme_window':seq[351:361]})
table(r/'05_qc/baseline_motif_alignment_audit.tsv',rows)
summary={label:{'N':sum(x['dataset']==label for x in rows),'all_three_anchor_matches':sum(x['dataset']==label and x['anchor_match'] for x in rows)} for label in ('v1_frozen','v1_curated_control')}
(r/'05_qc/baseline_motif_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
status=json.loads((r/'05_qc/run_status.json').read_text())
status['user_email_received']=True
status['stage']='awaiting_NCBI_web_BLAST_result'
status['search_status']='first_BLASTplus_attempt_timed_out; web_API_RID_DBX23EUJ014_waiting'
status['pending']=[x for x in status['pending'] if 'email' not in x]
(r/'05_qc/run_status.json').write_text(json.dumps(status,indent=2)+'\n')
print(json.dumps(summary,indent=2))
