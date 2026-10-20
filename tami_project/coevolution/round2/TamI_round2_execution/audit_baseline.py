"""Verify reproducibility and expose the baseline's functional scope."""
from pathlib import Path
import csv, json, hashlib, importlib.util, collections, shutil

r = Path.home()/'projects/tami_dca/TamI_round2'
repo = Path('/mnt/c/Users/31441/Desktop/work/BEEP')
spec = importlib.util.spec_from_file_location('tami',r/'scripts/tami_minimal_v1.py')
tami = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tami)
base = r/'00_baseline'
qc = r/'05_qc'
old = json.loads((base/'TamI_v1_input_qc.json').read_text())
repro = qc/'baseline_prepare_reproduction'
new = json.loads((repro/'input_qc.json').read_text())
checks = {
    'N_equal': old['N_model']==new['model_rows']==69,
    'L_equal': old['L_model']==new['model_length']==384,
    'Neff_equal': abs(old['Neff']-new['neff'])<1e-12,
    'reference_hash_equal': old['reference_sequence_sha256']==new['reference_sequence_sha256'],
    'model_fasta_byte_identical': (base/'TamI_v1_focus_cov70.faa').read_bytes()==(repro/'focus_cov70.faa').read_bytes(),
    'weights_byte_identical': (base/'TamI_v1_weights.txt').read_bytes()==(repro/'weights.txt').read_bytes(),
}
for item in json.loads((base/'source_manifest.json').read_text()):
    checks['unchanged_source:'+item['source']] = hashlib.sha256(Path(item['source']).read_bytes()).hexdigest()==item['sha256']
assert all(checks.values()),checks
(qc/'baseline_reproduction_checks.json').write_text(json.dumps(checks,indent=2)+'\n')

baseline = tami.read_fasta(base/'TamI_v1_focus_all.faa')
original = dict(tami.read_fasta(base/'TamI_v1_original_full_msa.faa'))
unaligned = dict(tami.read_fasta(base/'TamI_v1_original_unaligned.faa'))
ref = original['TamI_REF_ADC79647.1']
columns = [i for i,a in enumerate(ref) if a!='-']
restore = []
for sid,focus in baseline:
    full = original[sid].replace('-','')
    assert ''.join(original[sid][i] for i in columns)==focus
    assert full==unaligned[sid]
    restore.append({'sequence_id':sid,'full_length':len(full),'focus_ungapped_length':len(focus.replace('-','')),'residues_lost_if_ungapping_focus':len(full)-len(focus.replace('-','')),'original_full_sequence_verified':True})
tami.table(qc/'baseline_full_length_audit.tsv',list(restore[0]),restore)
tami.write_fasta(base/'TamI_v1_verified_full_length.faa',[(sid,unaligned[sid]) for sid,_ in baseline])

annotation_source=repo/'tami_project/data/raw/homologs/TamI/TamI_blast_filtered.tsv'
shutil.copyfile(annotation_source,base/'TamI_v1_blast_annotations.tsv')
with annotation_source.open(encoding='utf-8-sig') as f:
    annotation={row['saccver']:row for row in csv.DictReader(f,delimiter='\t')}
# Evidence categories are explicit annotation strings, not inferred orthology.
other_function_terms = ['PikC','Vitamin D','deoxyerythronolide','PksS','Mevastatin','Nocardicin',
    'Mycinamicin','Biotin biosynthesis','Pentalenic','Neocarzinostatin','Cyclomarin','Biflaviolin',
    'TylH1','Steroid','Vancomycin','Pentalenolactone','Erythromycin','Epothilone',
    'Methyl-branched lipid','dihydromenaquinone','Linalool','methoxybenzoate',
    'Pulcherriminic','Azinomycin','DoxA']
audit=[]
for sid,_ in baseline:
    source=annotation.get(sid,{})
    title=source.get('stitle','TamI reference ADC79647.1')
    evidence=[term for term in other_function_terms if term.lower() in title.lower()]
    category = 'reference_TamI' if sid=='TamI_REF_ADC79647.1' else ('explicit_other_function_annotation' if evidence else 'P450_annotation_does_not_establish_TamI_function')
    audit.append({'sequence_id':sid,'scope_category':category,'annotation_evidence_terms':';'.join(evidence),'organism':source.get('sscinames','Streptomyces sp. 307-9'),'taxid':source.get('staxids',''),'identity_to_reference_percent':source.get('pident','100'),'annotation':title,'source':'existing_v1_BLAST_annotation_table','action_for_TamI_focused_set':'reference' if category=='reference_TamI' else ('exclude_after_scope_decision' if evidence else 'requires_orthology_and_BGC_review')})
tami.table(qc/'baseline_function_scope_audit.tsv',list(audit[0]),audit)
summary={'baseline_N':69,'baseline_L':384,'baseline_Neff':old['Neff'],'baseline_Neff_over_L':old['Neff_over_L_model'],'reproduction_checks_pass':all(checks.values()),'annotation_categories':dict(collections.Counter(a['scope_category'] for a in audit)),'full_length_sequences_shortened_by_focus_projection':sum(x['residues_lost_if_ungapping_focus']>0 for x in restore),'full_length_restoration_verified_for_all_69':True,'stop_condition':'manual section 7: clearly many functionally different P450s in inherited baseline','scope_decision_pending':True,'round2_MSA_computed':False,'round2_Neff_computed':False}
(qc/'baseline_audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
