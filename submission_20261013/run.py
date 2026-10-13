#!/usr/bin/env python3
"""Run the current submission modules against an authorized local data tree."""
import argparse,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('module',choices=['dm','matched']);p.add_argument('--data-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--rscript',default='Rscript');p.add_argument('--check-inputs',action='store_true')
a=p.parse_args();root=a.data_root.resolve(strict=True);out=a.output.resolve()
if out==root or root in out.parents and any(t in out.parts for t in ['inputs_private_20261007','private_materialized']):p.error('Output must be separate from source inputs.')
if out.exists() and any(out.iterdir()):p.error('Output must be new or empty.')
for d in ['results','tables','logs','qa','output/pdf','staging','input']:(out/d).mkdir(parents=True,exist_ok=True)
env=dict(os.environ,DISCO_RESULTS_ROOT=str(root),DISCO_RUN_OUTPUT=str(out))
base=Path(__file__).resolve().parent/a.module
for f in base.glob('*.R'):subprocess.run([a.rscript,'-e',"parse(commandArgs(TRUE)[1]);cat('Parsed\\n')",str(f)],check=True,stdout=subprocess.DEVNULL,env=env)
if a.check_inputs:
    # Input paths are enumerated without loading participant records.
    import re,json
    src='\n'.join(f.read_text() for f in sorted(base.glob('*.R')))
    paths=[root/s for s in re.findall(r"file.path\(root,'([^']+)'\)",src) if '.' in Path(s).name]
    if a.module=='dm':
        private=root/'03_multimodal_UKB/inputs_private_20261007/private_materialized/data/hm/UKB'
        paths += [private/s for s in re.findall(r"file.path\(private,'([^']+)'\)",src)]
    paths=sorted(set(paths))
    missing=[str(s.relative_to(root)) for s in paths if not s.is_file()]
    result={'module':a.module,'checked_files':len(paths),'missing':missing,'statistical_models_run':False}
    (out/'results/input_check.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
    raise SystemExit(bool(missing))
for name in (['fit_logDM.R','plot_Fig3c.R'] if a.module=='dm' else ['fit_matched.R','plot_matched.R']):
    with (out/'logs'/f'{name}.log').open('w') as log:subprocess.run([a.rscript,str(base/name)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
print(out)
