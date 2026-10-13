"""Run a dated correction module without writing to the input tree."""
from pathlib import Path
import argparse, os, shutil, subprocess

parser=argparse.ArgumentParser()
parser.add_argument('module',choices=['organ_abc','organ_d','protein_sensitivity','hrs_scores','hrs_audit','math_geometry'])
parser.add_argument('--data-root',type=Path,required=True,help='Approved DISCO input tree, including Results_final snapshots')
parser.add_argument('--output',type=Path,required=True,help='A new output directory, outside the input snapshots')
parser.add_argument('--rscript',default='Rscript')
a=parser.parse_args()
code_root=Path(__file__).resolve().parent
data_root=a.data_root.resolve(strict=True)
out=a.output.resolve()
if out.exists(): raise SystemExit('Output already exists; use a new directory to protect prior results.')
if str(out).startswith(str(data_root)+os.sep) and 'Results_final' not in out.parts:
    raise SystemExit('Choose an output directory outside the raw-input directories.')
out.mkdir(parents=True)
env=os.environ.copy()
env.update(DISCO_CODE_ROOT=str(code_root),DISCO_DATA_ROOT=str(data_root),DISCO_RUN_OUTPUT=str(out),
           OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
jobs={'organ_abc':['organ_abc_compute.R','organ_plot_source.R','organ_joint_PH.R'],
      'organ_d':['organ_d_compute.R','organ_d_plot.R'],
      'protein_sensitivity':['protein_sensitivity_compute.R','protein_sensitivity_plot.R'],
      'hrs_scores':['hrs_median_compute.R'], 'hrs_audit':['hrs_independent_audit.R'],
      'math_geometry':['geometry_corrected.R']}
for source in jobs[a.module]:
    dest=out/source
    shutil.copy2(code_root/'sources'/source,dest)
    extra=['pass_1'] if a.module=='math_geometry' else []
    with (out/(source+'.log')).open('w') as log:
        process=subprocess.run([a.rscript,str(dest),*extra],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT)
    if process.returncode: raise SystemExit(f'{source} failed: see {source}.log')
print('Completed',a.module,'in',out)
