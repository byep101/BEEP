from pathlib import Path
import os, subprocess, json, datetime, shutil

r = Path.home()/'projects/tami_dca/TamI_round2'
p = r/'tools/prefix/usr'
env = os.environ.copy()
env['LD_LIBRARY_PATH'] = str(p/'lib/ncbi-blast+')+':'+str(p/'lib/x86_64-linux-gnu')
env['MAFFT_BINARIES'] = str(p/'lib/mafft/lib/mafft')
env['NCBI_DATA'] = str(p/'share/ncbi/data')
versions = {}
for binary, flag in [('blastp','-version'),('mafft','--version')]:
    result = subprocess.run([str(p/'bin'/binary),flag],env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
    versions[binary] = result.stdout.strip()
(r/'05_qc/tool_versions.json').write_text(json.dumps(versions,indent=2)+'\n')
print(versions,flush=True)
out = r/'01_search/TamI_refseq_bacteria_blast.tsv'
if out.exists():
    raise RuntimeError('Search file already exists; inspect before retrying')
cmd = [str(p/'bin/blastp'),'-query',str(r/'00_baseline/TamI_reference.fasta'),'-db','refseq_protein','-remote','-entrez_query','Bacteria[Organism]','-max_target_seqs','900','-evalue','1e-20','-outfmt','6 qseqid saccver pident length qstart qend sstart send evalue bitscore qlen slen','-out',str(out)]
record = {'argv':cmd,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'database_version':'remote current; server database metadata not exposed by this TSV format','max_target_seqs_is_return_limit':True,'timeout_seconds':1200}
(r/'01_search/search_record.json').write_text(json.dumps(record,indent=2)+'\n')
print('Starting one NCBI remote BLAST request',flush=True)
with (r/'logs/remote_blast.stderr.log').open('w') as err:
    try:
        result = subprocess.run(cmd,env=env,stdout=err,stderr=subprocess.STDOUT,timeout=1200)
        record['returncode'] = result.returncode
    except subprocess.TimeoutExpired:
        record['returncode'] = 'timeout'
record['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
record['output_bytes'] = out.stat().st_size if out.exists() else 0
(r/'01_search/search_record.json').write_text(json.dumps(record,indent=2)+'\n')
print(record,flush=True)
