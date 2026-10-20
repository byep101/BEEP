"""Reconstructed, audited round-2 intake. Not the missing manual ZIP's script.

NCBI GenBank XML supplies complete proteins, definition, taxonomy and taxid.
The best-scoring HSP per accession is used for identity and query-span coverage.
Retained records remain provisional homologs; this does not establish orthology.
"""
from pathlib import Path
import argparse, csv, hashlib, json, re, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

AA=set('ACDEFGHIKLMNPQRSTVWY')
FIELDS='qseqid saccver pident length qstart qend sstart send evalue bitscore qlen slen'.split()
OTHER_FUNCTION_TERMS=['pikc','vitamin d','deoxyerythronolide','pkss','mevastatin',
 'nocardicin','mycinamicin','biotin biosynthesis','pentalenic','neocarzinostatin',
 'cyclomarin','biflaviolin','tylh1','steroid','vancomycin','pentalenolactone',
 'erythromycin','epothilone','methyl-branched lipid','dihydromenaquinone',
 'linalool','methoxybenzoate','pulcherriminic','azinomycin','doxa','gfsf','fd-891','fd-892',
 'sulfonylurea','terpineol','thiocarbamate','pravastatin']

def read_fasta(path):
    records=[]; name=None; chunks=[]
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        if line.startswith('>'):
            if name is not None: records.append((name,''.join(chunks).upper()))
            name=line[1:].split()[0]; chunks=[]
        elif line.strip(): chunks.append(line.strip())
    if name is not None: records.append((name,''.join(chunks).upper()))
    if len({a for a,b in records})!=len(records): raise ValueError('Duplicate FASTA identifiers')
    return records

def fasta(path,records):
    with Path(path).open('w',encoding='ascii',newline='\n') as out:
        for sid,seq in records:
            out.write('>'+sid+'\n')
            for i in range(0,len(seq),80): out.write(seq[i:i+80]+'\n')

def table(path,rows):
    if not rows: return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def select_hits(path):
    best={}; count=0
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'): continue
        values=line.split('\t')
        if len(values)!=12: raise ValueError('Expected exactly the manual\'s 12 BLAST fields')
        h=dict(zip(FIELDS,values)); count+=1
        for key in ('pident','evalue','bitscore'): h[key]=float(h[key])
        for key in ('length','qstart','qend','sstart','send','qlen','slen'): h[key]=int(h[key])
        if h['qlen']!=413: raise ValueError('Query length is not 413')
        if not (1<=h['qstart']<=h['qend']<=413): raise ValueError('Bad query coordinates')
        if h['saccver'] not in best or h['bitscore']>best[h['saccver']]['bitscore']:
            best[h['saccver']]=h
    audit=[]; selected=[]
    for accession,h in best.items():
        coverage=(h['qend']-h['qstart']+1)/h['qlen']
        ratio=h['slen']/413
        failures=[]
        if not 25<=h['pident']<=95: failures.append('identity_outside_25_95')
        if coverage<0.7: failures.append('query_span_coverage_lt_0.7')
        if not 0.75<=ratio<=1.25: failures.append('length_ratio_outside_0.75_1.25')
        if h['evalue']>1e-20: failures.append('evalue_gt_1e-20')
        row={**h,'query_span_coverage':coverage,'length_ratio':ratio,'numeric_decision':'exclude' if failures else 'pass','reason':';'.join(failures)}
        audit.append(row)
        if not failures: selected.append(accession)
    return count,audit,selected

def parse_genbank(xml):
    root=ET.fromstring(xml)
    if root.tag!='GBSet': raise ValueError('NCBI response is not a GBSet')
    out={}
    for gb in root.findall('GBSeq'):
        acc=gb.findtext('GBSeq_accession-version')
        seq=(gb.findtext('GBSeq_sequence') or '').upper()
        if not acc or not seq: raise ValueError('Missing accession or complete sequence')
        declared=gb.findtext('GBSeq_length')
        if declared and int(declared)!=len(seq): raise ValueError('Downloaded record length mismatch: '+acc)
        quals=gb.findall('.//GBQualifier')
        taxid=[q.findtext('GBQualifier_value') for q in quals if q.findtext('GBQualifier_name')=='db_xref' and (q.findtext('GBQualifier_value') or '').startswith('taxon:')]
        coded=[q.findtext('GBQualifier_value') for q in quals if q.findtext('GBQualifier_name')=='coded_by']
        out[acc]={'sequence':seq,'annotation':gb.findtext('GBSeq_definition',''),'organism':gb.findtext('GBSeq_organism',''),'taxonomy':gb.findtext('GBSeq_taxonomy',''),'taxid':';'.join(taxid),'coded_by':';'.join(coded)}
    return out

def run(args):
    if not re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+',args.email) or 'example.' in args.email.lower():
        raise ValueError('Provide the user\'s real email, not a placeholder')
    out=Path(args.out); out.mkdir(parents=True,exist_ok=False)
    cache=out/'ncbi_raw'; cache.mkdir()
    count,audit,accessions=select_hits(args.blast_tsv)
    table(out/'blast_numeric_audit.tsv',audit)
    if not accessions: raise ValueError('No hits pass the fixed numeric filters')
    raw={}; download_audit=[]
    for start in range(0,len(accessions),50):
        batch=accessions[start:start+50]
        params={'db':'protein','id':','.join(batch),'rettype':'gp','retmode':'xml','tool':'TamI_round2_audit','email':args.email}
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        request=urllib.request.Request(url,data=urllib.parse.urlencode(params).encode(),headers={'User-Agent':'TamI-round2-research/1.0'})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request,timeout=60) as response: data=response.read()
                parsed=parse_genbank(data)
                missing=set(batch)-set(parsed)
                if missing: raise ValueError('Missing requested versions: '+','.join(sorted(missing)))
                break
            except Exception:
                if attempt==2: raise
                time.sleep(2*(attempt+1))
        raw.update(parsed)
        name=f'batch_{start//50+1:03d}.xml'; (cache/name).write_bytes(data)
        download_audit.append({'batch':name,'requested':len(batch),'received':len(parsed),'sha256':hashlib.sha256(data).hexdigest(),'accessions':';'.join(batch)})
        table(out/'download_audit.tsv',download_audit)
        print('Downloaded complete proteins:',len(raw),flush=True)
        time.sleep(0.5)
    baseline=read_fasta(args.baseline)
    reference=read_fasta(args.reference)
    if len(reference)!=1 or len(reference[0][1])!=413: raise ValueError('Bad reference')
    refseq=reference[0][1]
    if hashlib.sha256(refseq.encode()).hexdigest()!='fbc563a7aeebba572bac52607531eb227351b690037bf8108f176a3f8701027a': raise ValueError('Reference hash mismatch')
    combined=[('TamI_REF_ADC79647.1',refseq)]
    excluded_sequences={seq:sid for sid,seq in read_fasta(args.excluded_fasta)} if args.excluded_fasta else {}
    seen={refseq:'TamI_REF_ADC79647.1'}; intake=[]
    for sid,seq in baseline:
        if sid=='TamI_REF_ADC79647.1': continue
        if '-' in seq: raise ValueError('Baseline must be verified full-length unaligned proteins')
        duplicate=seen.get(seq,'')
        intake.append({'sequence_id':sid,'origin':'curated_baseline','decision':'exclude_duplicate' if duplicate else 'retain_provisional','duplicate_of':duplicate,'annotation_review':'see baseline curation audit; orthology unresolved'})
        if not duplicate: combined.append((sid,seq)); seen[seq]=sid
    downloaded=[]
    blast_lengths={row['saccver']:row['slen'] for row in audit}
    for sid in accessions:
        rec=raw[sid]; seq=rec['sequence']; downloaded.append((sid,seq))
        failures=[]; nonstd=sum(a not in AA for a in seq)/len(seq)
        if len(seq)!=blast_lengths[sid]: raise ValueError('BLAST/downloaded full-length mismatch for '+sid)
        if not 0.75<=len(seq)/413<=1.25: failures.append('downloaded_length_ratio')
        if nonstd>0.01: failures.append('nonstandard_fraction_gt_0.01')
        if 'Bacteria' not in [x.strip() for x in rec['taxonomy'].split(';')]: failures.append('bacterial_lineage_not_verified')
        function_terms=[s for s in OTHER_FUNCTION_TERMS if s in rec['annotation'].lower()]
        if function_terms: failures.append('explicit_other_function_annotation')
        if seq in excluded_sequences: failures.append('exact_sequence_matches_excluded_baseline:'+excluded_sequences[seq])
        if re.search(r'(.)\1{9,}',seq): failures.append('homopolymer_at_least_10_requires_review')
        if seq in seen: failures.append('exact_sequence_duplicate')
        heme=[m.start()+1 for m in re.finditer(r'F..G...C.G',seq)]
        exxr=[m.start()+1 for m in re.finditer(r'E..R',seq)]
        decision='exclude' if failures else 'retain_provisional'
        row={'sequence_id':sid,'origin':'RefSeq','decision':decision,'reason':';'.join(failures),'duplicate_of':seen.get(seq,''),'length':len(seq),'nonstandard_fraction':nonstd,'other_function_terms':';'.join(function_terms),'heme_pattern_starts':';'.join(map(str,heme)),'ExxR_starts':';'.join(map(str,exxr)),**{k:v for k,v in rec.items() if k!='sequence'},'orthology_status':'not_established','BGC_status':'not_checked'}
        intake.append(row)
        if not failures: combined.append((sid,seq)); seen[seq]=sid
    if len(combined)>1000: raise ValueError('More than 1000 rows: do not truncate')
    fasta(out/'downloaded_full_length.faa',downloaded)
    fasta(out/'combined_curated_provisional_v1.faa',combined)
    table(out/'intake_audit.tsv',intake)
    summary={'raw_HSP_rows':count,'unique_accessions':len(audit),'numeric_pass':len(accessions),'downloaded':len(raw),'baseline_rows_including_reference':len(baseline),'combined_rows':len(combined),'new_retained':sum(x['origin']=='RefSeq' and x['decision']=='retain_provisional' for x in intake),'curation_status':'provisional; inspect annotations, motifs and taxonomy before alignment','BGC_orthology_not_proven':True,'email_supplied':True}
    (out/'intake_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('blast-tsv','baseline','reference','out','email'): p.add_argument('--'+name,required=True)
    p.add_argument('--excluded-fasta',help='Prevent reintroduction of known other-function baseline proteins under different accessions')
    run(p.parse_args())
