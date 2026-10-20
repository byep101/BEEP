"""Read current primary UniProt annotations for the 28 unresolved old proteins."""
from pathlib import Path
import csv, io, json, urllib.parse, urllib.request, datetime, hashlib
from collect_and_filter import read_fasta
r=Path.home()/'projects/tami_dca/TamI_round2'
rows=read_fasta(r/'02_intake/baseline_curated_v1/baseline_curated_unaligned.faa')
accessions=[sid.rsplit('.',1)[0] for sid,seq in rows if not sid.startswith('TamI_REF')]
params={'query':'('+' OR '.join('accession:'+x for x in accessions)+')','format':'tsv','fields':'accession,id,protein_name,gene_names,organism_name,cc_function,cc_pathway,cc_catalytic_activity,sequence'}
url='https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode(params)
out=r/'05_qc/baseline_uniprot_review';out.mkdir(exist_ok=False)
with urllib.request.urlopen(url,timeout=60) as response: data=response.read();headers=dict(response.headers)
(out/'uniprot_annotations.tsv').write_bytes(data)
records=list(csv.DictReader(io.StringIO(data.decode()),delimiter='\t'))
record={'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'url':url,'requested_accessions':accessions,'received_records':len(records),'sha256':hashlib.sha256(data).hexdigest(),'response_headers':headers}
(out/'retrieval_record.json').write_text(json.dumps(record,indent=2)+'\n')
byacc={sid.rsplit('.',1)[0]:seq for sid,seq in rows if not sid.startswith('TamI_REF')}
for row in records:
    print(json.dumps({k:v for k,v in row.items() if k!='Sequence'},ensure_ascii=False))
    assert row['Sequence']==byacc[row['Entry']],('Sequence mismatch',row['Entry'])
print('Exact full-length sequence matches for',len(records),'entries')
