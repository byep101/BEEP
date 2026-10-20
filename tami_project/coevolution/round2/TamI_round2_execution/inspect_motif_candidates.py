from pathlib import Path
import re,json
from collect_and_filter import read_fasta
r=Path.home()/'projects/tami_dca/TamI_round2'
rows=read_fasta(r/'00_baseline/TamI_v1_focus_all.faa'); ref=dict(rows)['TamI_REF_ADC79647.1']
out=[]
for m in re.finditer('E..R',ref):
    i=m.start()
    out.append({'reference_E_position':i+1,'reference_R_position':i+4,'reference_window':m.group(),'E_matches':sum(s[i]=='E' for _,s in rows),'R_matches':sum(s[i+3]=='R' for _,s in rows),'both_match':sum(s[i]=='E' and s[i+3]=='R' for _,s in rows),'N':len(rows)})
print(json.dumps(out,indent=2))
(r/'05_qc/reference_ExxR_candidate_audit.json').write_text(json.dumps(out,indent=2)+'\n')
