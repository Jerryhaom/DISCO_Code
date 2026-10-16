"""Run retained analysis sources in a new local output directory."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys
import hashlib

ROOT=Path(__file__).resolve().parent
JOBS={
 'theory-identities':('theory',['verify_identities.py']),
 'theory-geometry':('theory',['check_actual_geometry.R']),
 'ed4-plot':('ed4',['plot_recomputed.R']),
 'ed4f':('ed4',['run_ed4f.R']),
 'ed4h':('ed4',['run_ed4h.R']),
 'fig4c-plot':('fig4c',['plot_verified_raw_branch.R']),
 'fig4c':('fig4c',['run_fig4c.R','compare_full_raw_branch.R','plot_verified_raw_branch.R']),
 's11b':('s11b',['01_recompute_raw_and_pairs.R','02_current_endpoint_targets.R','04_server_core_fingerprint.R','05_verify_id_score_bridge.R','03_power_statistics_and_plot.R']),
 'ed9c':('ed9c',['run_ed9c.R']),
 'ph-check':('ph',['check_frozen_PH.py']),
 'ed1-contract':('source_data',['inspect_ed1_contract.py']),
 'ed1-observe':('source_data',['export_ED1_original.R']),
 'ed1-manifest-check':('source_data',['inspect_ed1_contract.py']),
 'charls':('charls',['build_interval_endpoints.R','build_retrospective_sensitivity.R','fit_interval_sensitivity.R:primary:200:observed_wave','fit_interval_sensitivity.R:repeat:200:observed_wave','fit_interval_sensitivity.R:retrospective_primary:200:retrospective_strict','fit_interval_sensitivity.R:retrospective_repeat:200:retrospective_strict','summarize_and_verify.R','plot_interval_sensitivity.R']),
}
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('mode',choices=JOBS)
 p.add_argument('--output',type=Path,required=True,help='New directory; existing directories are rejected')
 p.add_argument('--config',type=Path,help='Explicit local JSON configuration for controlled-input modes')
 p.add_argument('--rscript',default='Rscript')
 a=p.parse_args();key,jobs=JOBS[a.mode]
 config=json.loads(a.config.read_text()) if a.config else {}
 contract=json.loads((ROOT/'required_inputs.json').read_text())[a.mode]
 missing=[k for k in contract['required_config_keys'] if not config.get(k)]
 if missing:raise SystemExit('Missing explicit configuration: '+', '.join(missing))
 for k in contract['required_config_keys']:
  value=Path(config[k]).expanduser().resolve()
  if not value.exists():raise SystemExit('Configured input does not exist: '+k)
  config[k]=str(value)
 for k in contract['required_config_keys']:
  if config.get(k):
   base=Path(config[k]).resolve()
   for rel in contract.get('required_'+k+'_paths',[]):
    if not (base/rel).exists():raise SystemExit('Required input is missing: '+k+'/'+rel)
 out=a.output.resolve()
 if out==ROOT.parent or ROOT.parent in out.parents:raise SystemExit('Output must be outside the code repository.')
 if a.mode=='ed1-observe':
  source_contract=json.loads((ROOT/'modules/source_data/original_source_contract.json').read_text())
  actual=hashlib.sha256(Path(config['simulation_source']).read_bytes()).hexdigest()
  if actual!=source_contract['whole_source_sha256']:raise SystemExit('Caller Simulation source hash differs from the recorded original source; inspect before running.')
 if out.exists():raise SystemExit('Output already exists; choose a new directory.')
 for k in contract['required_config_keys']:
  src=Path(config[k])
  if src.is_dir() and (out==src or src in out.parents):raise SystemExit('Output must be outside configured input roots.')
 shutil.copytree(ROOT/'modules'/key,out)
 for d in ('input/private','results','logs','figures','tables','compile_cache','private/input','private/models'):(out/d).mkdir(parents=True,exist_ok=True)
 for d in ('input/private','private','private/input','private/models'):os.chmod(out/d,0o700)
 env=os.environ.copy();env.update(DISCO_RELEASE_ROOT=str(ROOT),DISCO_RUN_MODULE=str(out),DISCO_INPUT_ROOT=config.get('input_root',''),DISCO_RESULTS_INPUT=config.get('results_root',''),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 for k,v in {'geometry_input':'DISCO_GEOMETRY_INPUT','organ_input':'DISCO_ORGAN_INPUT','mortality_endpoint':'DISCO_MORTALITY_ENDPOINT','disease_endpoint':'DISCO_DISEASE_ENDPOINT','organ_disease_patch':'DISCO_ORGAN_DISEASE_PATCH','theory_font':'DISCO_THEORY_FONT','simulation_source':'DISCO_SIMULATION_SOURCE','r_library':'DISCO_R_LIBRARY','ed1_manifest':'DISCO_ED1_MANIFEST'}.items():
  if config.get(k):env[v]=config[k]
 if a.mode=='ed1-observe':env['DISCO_EXPORT_SCOPE']='ED1_original'
 log=[]
 for job in jobs:
  name,*extra=job.split(':');source=out/'code'/name
  args=extra or ([str(out/'results')] if name in ('02_current_endpoint_targets.R','04_server_core_fingerprint.R') else [str(out)])
  command=[a.rscript,str(source),*args] if name.endswith('.R') else [sys.executable,str(source)]
  dest=out/'logs'/(job.replace(':','_')+'.log')
  with dest.open('w') as f:r=subprocess.run(command,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT)
  log.append(dict(source=name,args=args,returncode=r.returncode,log=str(dest.relative_to(out))))
  (out/'execution.json').write_text(json.dumps(dict(mode=a.mode,scope=contract['execution_scope'],jobs=log,cohort_refits_executed=a.mode in {'ed4h','fig4c','s11b','charls'},published=False),indent=2)+'\n')
  if r.returncode:raise SystemExit('Source failed: '+str(dest))
 print('Completed '+a.mode+'; execution evidence: '+str(out/'execution.json'))
if __name__=='__main__':main()
