"""Sequential fallback to NCBI's web BLAST interface with an auditable RID."""
from pathlib import Path
import argparse, datetime, hashlib, json, re, urllib.parse, urllib.request, xml.etree.ElementTree as ET

ROOT=Path.home()/'projects/tami_dca/TamI_round2'
OUT=ROOT/'01_search/url_api_v1'
URL='https://blast.ncbi.nlm.nih.gov/Blast.cgi'

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def request(params):
    req=urllib.request.Request(URL,data=urllib.parse.urlencode(params).encode(),headers={'User-Agent':'TamI-round2-research/1.0'})
    with urllib.request.urlopen(req,timeout=90) as response: return response.read()
def save(record): (OUT/'request_record.json').write_text(json.dumps(record,indent=2)+'\n')

def submit(email):
    previous=json.loads((ROOT/'01_search/search_record.json').read_text())
    if 'returncode' not in previous: raise RuntimeError('Original remote request still running; no parallel submission')
    OUT.mkdir(exist_ok=False)
    params={'CMD':'Put','PROGRAM':'blastp','DATABASE':'refseq_protein','QUERY':(ROOT/'00_baseline/TamI_reference.fasta').read_text(),'ENTREZ_QUERY':'Bacteria[Organism]','EXPECT':'1e-20','HITLIST_SIZE':'900','EMAIL':email,'TOOL':'TamI_round2_audit','MATRIX':'BLOSUM62','WORD_SIZE':'3','GAPCOSTS':'11 1','COMPOSITION_BASED_STATISTICS':'2','FILTER':'F'}
    # blastp CLI default gap extension is 1; record all explicit API parameters.
    record={'started_utc':now(),'method':'NCBI web BLAST Common URL API fallback','prior_attempt':str(ROOT/'01_search/search_record.json'),'parameters':{k:v for k,v in params.items() if k not in ('QUERY','EMAIL')},'email_provided':True,'reference_sha256':hashlib.sha256(params['QUERY'].encode()).hexdigest()}
    save(record)
    data=request(params); (OUT/'submit_response.html').write_bytes(data)
    text=data.decode('utf-8','replace')
    rid=re.search(r'RID\s*=\s*([A-Z0-9-]+)',text)
    if not rid: raise RuntimeError('No RID returned; inspect saved response')
    record['RID']=rid.group(1); record['submitted_utc']=now()
    rtoe=re.search(r'RTOE\s*=\s*(\d+)',text)
    record['RTOE_seconds']=int(rtoe.group(1)) if rtoe else None
    save(record); print(json.dumps(record,indent=2))

def poll():
    record=json.loads((OUT/'request_record.json').read_text())
    last=record.get('last_checked_utc',record['submitted_utc'])
    elapsed=(datetime.datetime.now(datetime.timezone.utc)-datetime.datetime.fromisoformat(last)).total_seconds()
    if elapsed<60:
        print('POLL_NOT_DUE; wait at least',round(60-elapsed+1),'seconds'); return
    data=request({'CMD':'Get','RID':record['RID'],'FORMAT_OBJECT':'SearchInfo'})
    (OUT/'latest_search_info.txt').write_bytes(data)
    text=data.decode('utf-8','replace')
    status=re.search(r'Status=(\w+)',text)
    record['last_checked_utc']=now(); record['status']=status.group(1) if status else 'UNKNOWN_RESPONSE'
    save(record); print(json.dumps(record,indent=2))
    print('\n'.join(line.strip() for line in text.splitlines() if any(x in line for x in ('Status=','ThereAreHits=','Error:','Message ID#'))))

def retrieve():
    record=json.loads((OUT/'request_record.json').read_text())
    data=request({'CMD':'Get','RID':record['RID'],'FORMAT_TYPE':'XML2_S','DESCRIPTIONS':'900','ALIGNMENTS':'900'})
    (OUT/'result.xml').write_bytes(data)
    tree=ET.fromstring(data)
    # Remove XML namespace only from parsed tags; retain untouched raw response.
    for element in tree.iter(): element.tag=element.tag.split('}')[-1]
    search=tree.find('.//Search')
    if search is None: raise ValueError('No Search element in BLAST XML2')
    qlen=int(search.findtext('query-len'))
    if qlen!=413: raise ValueError('Returned query length differs from TamI reference')
    report=tree.find('.//Report')
    database=report.findtext('./search-target/Target/db')
    if database!='refseq_protein': raise ValueError('Unexpected returned database: '+str(database))
    output=[]; meta=[]
    for hit in search.findall('./hits/Hit'):
        descr=hit.find('./description/HitDescr')
        acc=descr.findtext('accession')
        identifier=descr.findtext('id','')
        versioned=re.search(r'(?:ref\|)?([A-Z]{1,4}_?\d+\.\d+)',identifier)
        if versioned: acc=versioned.group(1)
        if not re.search(r'\.\d+$',acc): raise ValueError('BLAST did not supply a versioned accession: '+str(acc)+' '+identifier)
        slen=int(hit.findtext('len'))
        meta.append({'accession':acc,'title':descr.findtext('title'),'taxid':descr.findtext('taxid'),'scientific_name':descr.findtext('sciname'),'length':slen,'all_descriptions':[{el.tag:el.text for el in d} for d in hit.findall('./description/HitDescr')]})
        for hsp in hit.findall('./hsps/Hsp'):
            length=int(hsp.findtext('align-len')); ident=int(hsp.findtext('identity'))
            output.append(['ADC79647.1',acc,f'{100*ident/length:.6f}',length,hsp.findtext('query-from'),hsp.findtext('query-to'),hsp.findtext('hit-from'),hsp.findtext('hit-to'),hsp.findtext('evalue'),hsp.findtext('bit-score'),qlen,slen])
    (OUT/'TamI_refseq_bacteria_blast.tsv').write_text(''.join('\t'.join(map(str,row))+'\n' for row in output))
    (OUT/'hit_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    record['retrieved_utc']=now(); record['unique_hits']=len(meta); record['HSP_rows']=len(output)
    record['result_sha256']=hashlib.sha256(data).hexdigest()
    record['database_statistics']={x.tag:x.text for x in search.findall('./stat/Statistics/*')}
    record['server_blast_version']=report.findtext('version')
    record['returned_database']=database
    record['returned_parameters']={x.tag:x.text for x in report.findall('./params/Parameters/*')}
    record['status']='RETRIEVED'; save(record)
    print(json.dumps(record,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('action',choices=['submit','poll','retrieve']); p.add_argument('--email')
    a=p.parse_args()
    if a.action=='submit':
        if not a.email: p.error('--email required')
        submit(a.email)
    elif a.action=='poll': poll()
    else: retrieve()
