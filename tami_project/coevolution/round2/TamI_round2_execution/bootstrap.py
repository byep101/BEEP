"""Stage a separate WSL round-2 run without changing round 1."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import datetime

ROOT = Path.home() / 'projects/tami_dca'
R2 = ROOT / 'TamI_round2'
V1 = ROOT / 'TamI_v1'
REPO = Path('/mnt/c/Users/31441/Desktop/work/BEEP')
R2.mkdir(exist_ok=False)
for sub in ('00_baseline', '01_search', '02_intake', '03_msa', '04_prepare', '05_qc', 'logs', 'scripts', 'tools/packages', 'tools/prefix'):
    (R2/sub).mkdir(parents=True, exist_ok=True)
freeze = ['focus_all.faa', 'focus_cov70.faa', 'input_qc.json', 'weights.txt', 'column_map.tsv', 'MIp_pairs.csv', 'dca/DCA_pairs.csv']
manifest = []
for name in freeze:
    src = V1 / name
    dst = R2/'00_baseline'/('TamI_v1_' + src.name)
    shutil.copyfile(src, dst)
    manifest.append({'source': str(src), 'copy': str(dst), 'sha256': hashlib.sha256(src.read_bytes()).hexdigest()})
for src, name in [
    (REPO/'tami_project/data/raw/ncbi/ADC79647.1_TamI_Streptomyces_307-9.fasta', 'TamI_reference.fasta'),
    (REPO/'tami_project/coevolution/round1_pilot/TamI_coev_minimal_kit/inputs/TamI_full_msa.faa', 'TamI_v1_original_full_msa.faa'),
    (REPO/'tami_project/data/processed/msa/TamI/TamI_pilot_input.faa', 'TamI_v1_original_unaligned.faa'),
]:
    shutil.copyfile(src, R2/'00_baseline'/name)
shutil.copyfile(REPO/'tami_project/coevolution/round1_pilot/TamI_coev_minimal_kit/scripts/tami_minimal.py', R2/'scripts/tami_minimal_baseline_actual.py')
shutil.copyfile(ROOT/'TamI_DCA_EVcouplings_Continuation_StarterKit_v1.1/tami_dca_resume_v1_1/scripts/tami_minimal.py', R2/'scripts/tami_minimal_v1.py')
(R2/'00_baseline/source_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
(R2/'00_baseline/SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in sorted((R2/'00_baseline').iterdir()) if p.is_file()))
(R2/'05_qc/run_status.json').write_text(json.dumps({'status':'environment_preparation','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'homolog intake, MSA, weights, Neff only'},indent=2)+'\n')
print('Frozen baseline:', R2, flush=True)
packages = R2/'tools/packages'
with (R2/'logs/tool_download.log').open('w') as log:
    subprocess.run(['apt-get','download','mafft','ncbi-blast+','ncbi-data','libmbedcrypto7','libmbedtls14','libmbedx509-1'],cwd=packages,stdout=log,stderr=subprocess.STDOUT,check=True)
for package in sorted(packages.glob('*.deb')):
    subprocess.run(['dpkg-deb','-x',str(package),str(R2/'tools/prefix')],check=True)
print('Tools extracted into independent run directory', flush=True)
