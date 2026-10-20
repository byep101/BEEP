from pathlib import Path
import csv,json,subprocess,sys,datetime
from collect_and_filter import read_fasta,fasta,table
r=Path.home()/'projects/tami_dca/TamI_round2'
source=r/'05_qc/baseline_uniprot_review/uniprot_annotations.tsv'
with source.open() as f: entries={row['Entry']:row for row in csv.DictReader(f,delimiter='\t')}
additional={'E0D207':'GfsF: FD-891/FD-892 biosynthesis and oxidation',
 'P18327':'CYP105B1/SU2: sulfonylurea herbicide metabolism',
 'P33006':'CYP108: alpha-terpineol hydroxylation',
 'P43492':'CYP116/ThcB: thiocarbamate herbicide degradation',
 'Q59831':'CYP105A3: ML-236B carboxylate hydroxylation to pravastatin'}
with (r/'02_intake/baseline_curated_v1/curation_audit.tsv').open(encoding='utf-8-sig') as f: audit=list(csv.DictReader(f,delimiter='\t'))
for row in audit:
    acc=row['sequence_id'].rsplit('.',1)[0]
    entry=entries.get(acc,{})
    row['current_uniprot_function']=entry.get('Function [CC]','')
    row['current_uniprot_source']='https://www.uniprot.org/uniprotkb/'+acc+'/entry' if entry else ''
    row['additional_exclusion_reason']=additional.get(acc,'')
    if acc in additional: row['curation_decision']='excluded_additional_explicit_other_function'
    row['interpretation_note']='Unspecified, probable or in-vitro binding annotation alone is not confirmed TamI function or a definitive alternative physiological function' if entry and acc not in additional else ''
excluded={x['sequence_id'] for x in audit if x['curation_decision'].startswith('excluded')}
original=read_fasta(r/'00_baseline/TamI_v1_original_full_msa.faa')
kept=[(sid,seq) for sid,seq in original if sid not in excluded]
out=r/'02_intake/baseline_curated_v2';out.mkdir(exist_ok=False)
table(out/'curation_audit.tsv',audit)
fasta(out/'baseline_curated_unaligned.faa',[(sid,seq.replace('-','')) for sid,seq in kept])
fasta(out/'baseline_curated_original_alignment.faa',kept)
fasta(out/'excluded_baseline_full_length.faa',[(sid,seq.replace('-','')) for sid,seq in original if sid in excluded])
summary={'N_original':69,'initial_explicit_exclusions':40,'additional_explicit_exclusions':5,'N_excluded':len(excluded),'N_retained_including_reference':len(kept),'retained_nonreference_orthology_status':'unresolved provisional homologs','additional_evidence':additional,'annotation_snapshot':str(source),'previous_version_retained':'baseline_curated_v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(out/'curation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
subprocess.run([sys.executable,str(r/'scripts/tami_minimal_v1.py'),'prepare','--msa',str(out/'baseline_curated_original_alignment.faa'),'--reference',str(r/'00_baseline/TamI_reference.fasta'),'--ref-id','TamI_REF_ADC79647.1','--out',str(r/'05_qc/baseline_curated_prepare_control_v2')],check=True)
print(json.dumps(summary,indent=2))
