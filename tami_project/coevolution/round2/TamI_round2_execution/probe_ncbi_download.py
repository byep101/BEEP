from pathlib import Path
import urllib.parse,urllib.request,hashlib,json,datetime
from collect_and_filter import parse_genbank,read_fasta
r=Path.home()/'projects/tami_dca/TamI_round2'
out=r/'05_qc/ncbi_download_probe';out.mkdir(exist_ok=False)
params={'db':'protein','id':'ADC79647.1','rettype':'gp','retmode':'xml','tool':'TamI_round2_audit','email':'byep101@gmail.com'}
url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
req=urllib.request.Request(url,data=urllib.parse.urlencode(params).encode(),headers={'User-Agent':'TamI-round2-research/1.0'})
with urllib.request.urlopen(req,timeout=60) as response: data=response.read()
(out/'ADC79647.1.xml').write_bytes(data)
entry=parse_genbank(data)['ADC79647.1']
assert entry['sequence']==read_fasta(r/'00_baseline/TamI_reference.fasta')[0][1]
record={'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sequence_matches_frozen_reference':True,'reference_sequence_sha256':hashlib.sha256(entry['sequence'].encode()).hexdigest(),'raw_xml_sha256':hashlib.sha256(data).hexdigest(),**{k:v for k,v in entry.items() if k!='sequence'}}
(out/'probe_result.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
