"""Run approved dated analyses from immutable sources into a separate local output.

Participant inputs remain controlled; generated PRIVATE data are never release files.
This launcher does not certify provisional survival-date definitions or PH validity.
"""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess

MODULES = {
 'table1': [('recompute_table1.R', ['@OUT@/numeric'])],
 'ipaq': [('run_IPAQ_sensitivity.R', ['@DATA@', '@OUT@/numeric'])],
 'organ_cause': [('01_compute.R', ['pass_1']), ('02_plot.R', ['pass_1'])],
 'organ_cause_ph': [('01_compute.R', ['pass_1'])],
 'CIF': [('01_compute.R', ['pass_1']), ('02_plot.R', [])],
 'organ_70_30': [('01_prepare.R', []), ('02_train_clocks.py', []),
   ('03_DISCO_scores.R', []), ('04_compare.R', []), ('05_plot_d.R', []), ('06_plot_matched.R', [])],
 'organ_70_30_compare': [('04_compare.R', []), ('05_plot_d.R', []), ('06_plot_matched.R', [])],
 'nmr': [('nmr_dimension_diagnostic.R', []), ('repair_nmr_core.R', ['pass_1']),
   ('recalculate_downstream.R', ['pass_1'])],
}
CONTRACTS = {
 'table1': ['Result4/CHARLS.DISCO.rds','Result4/CLHLS.DISCO.rds',
   'Result4/CLHLS/New/CLHLS.2014_2018.longitudinal.datasets.rds','Result4/RuLAS.DISCO.rds',
   'Result1/NHANES.Clinical.DISCO2.rds','Result1/NHANES.Clinical.DISCO.rds','UKB_Data/NHANES_covariates.rds',
   'UKB_Data/extended_covariates_clinical.rds','UKB_Data/covariates_raw_full.rds','UKB_Data/ukb_race_21000.csv'],
 'ipaq': ['UKB_Data/extended_covariates_clinical.rds','UKB_Data/covariates_raw_full.rds','UKB_Data/ukb_race_21000.csv'],
 'organ_cause': ['Result3/Organ/UKB.odis.rds','Result3/Organ/UKB.prot.imputed.rds',
   'UKB_Data/extended_covariates_clinical.rds','UKB_Data/ukb_race_21000.csv','UKB_Data/UKB.Mortality.csv',
   'Result3/Organ/Plot_Fig4E.R','Result3/Organ/Tissue.order.rds','Result3/Organ/mortality/res.rds',
   'Result3/Organ/mortality/pl.rds',
   'Results_final/05_organ_network/outputs/organ_expanded_covariates/organ_cause_specific_expanded.csv',
   'Results_final/05_organ_network/outputs/organ_expanded_covariates/organ_chronic_disease_expanded.csv'],
 'CIF': ['Result2/UKB.ClinicalProt.DISCO.rds','UKB_Data/dat_disease.rds','UKB_Data/UKB.Mortality.csv',
   'Results_final/03_multimodal_UKB/reproduce_ED7_KM_facet.R'],
 'organ_70_30': ['UKB_Data/20240321_prot_instance0.RData','Result2/Proteins/coding143.tsv',
   'Result3/Organ/UKB.Prot.class.rds',
   'Results_final/05_organ_network/organ_clock_comparison/20261002_local_rebuild/PRIVATE/meta.rds',
   'Results_final/05_organ_network/organ_clock_comparison/20261002_local_rebuild/PRIVATE/common_cohorts.rds',
   'Results_final/05_organ_network/organ_clock_comparison/20261002_local_rebuild/results/feature_sets.json',
   'Results_final/05_organ_network/organ_clock_comparison/20260928/results/organ_mapping.csv',
   'Results_final/05_organ_network/organ_clock_comparison/20261002_rebuilt_figures/02_plot.R',
   'Results_final/05_organ_network/organ_clock_comparison/20261002_local_rebuild/vendor/clock_utils.py',
   'Results_final/05_organ_network/organ_clock_comparison/20261002_official_reproduction/vendor/Oh/UKB_Olink_3k_aging_model_github.ipynb',
   'Results_final/05_organ_network/organ_clock_comparison/20261002_official_reproduction/vendor/Goeminne/scripts/R/Olink_Processing.R'],
 'nmr': ['Result6/Cal_Rugao_omics.R','Result6/Rugao/NMR.csv','Result6/Rugao.Phe.DISCO2.rds',
   'Result6/Rugao.PCs.met.rds','Result6/Rugao.Met.DISCO.rds','Result6/Rugao.Met.DISCO2.rds',
   'Result6/Rugao.2019.mortality.rds','Result6/Plot_surv.R','Result6/Plot_subsystems2.R',
   'Result6/microbe/Species.txt','Result6/Rugao.Meta.DISCO2.rds','Result6/Rugao.selected.microbe.rds',
   'Result6/Rugao.subsystems.rds','Result6/Rugao.subsystems2.rds','Result6/Analysis.2019.x.xlsx',
   'Result6/Rugao/2019Rugao.data.y1.csv'],
}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('module',choices=MODULES)
    p.add_argument('--data-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True,help='New directory outside raw source directories')
    p.add_argument('--questionnaire',type=Path)
    p.add_argument('--rscript',default='Rscript')
    p.add_argument('--python',default='python3')
    p.add_argument('--preflight',action='store_true',help='Inspect dependency existence; no computation or individual output')
    p.add_argument('--frozen-score-snapshot',type=Path,help='Explicit new 70:30 score snapshot for comparison-only; never a substitute for full training')
    p.add_argument('--allow-provisional',action='store_true',help='Explicitly accept unreconciled saved date definitions for local numerical reproduction')
    a=p.parse_args(); code=Path(__file__).resolve().parent
    data=a.data_root.resolve(strict=True); out=a.output.resolve()
    if out.exists(): p.error('Output exists; use a new directory')
    if out==data or (data in out.parents and 'Results_final' not in out.relative_to(data).parts):
        p.error('Outputs cannot be placed in raw input trees')
    if code in out.parents: p.error('Generated participant outputs cannot be inside the upload candidate')
    source_module={'organ_70_30_compare':'organ_70_30','organ_cause_ph':'organ_cause'}.get(a.module,a.module)
    required=[data / x for x in CONTRACTS[source_module]]
    if a.module=='organ_70_30_compare':
        if a.frozen_score_snapshot is None: p.error('Comparison-only requires --frozen-score-snapshot')
        frozen=a.frozen_score_snapshot.resolve(strict=True)
        required += [frozen/x for x in ['PRIVATE/DISCO_original.rds','PRIVATE/clock_predictions.csv',
          'numeric/new_train_selected_proteins.json','qa/clocks_complete.json']]
    elif a.frozen_score_snapshot is not None: p.error('Frozen scores are permitted only in explicitly named comparison-only mode')
    questionnaire=a.questionnaire or (Path(os.environ['DISCO_RULAS_QUESTIONNAIRE']) if os.environ.get('DISCO_RULAS_QUESTIONNAIRE') else None)
    if a.module=='table1':
        if questionnaire is None: p.error('Table1 requires --questionnaire or DISCO_RULAS_QUESTIONNAIRE')
        required.append(questionnaire.resolve())
    missing=[str(x.relative_to(data)) if data in x.parents else x.name for x in required if not x.is_file()]
    source_hashes={str(x.relative_to(code)):hashlib.sha256(x.read_bytes()).hexdigest()
      for x in sorted((code/source_module).rglob('*')) if x.is_file()}
    out.mkdir(parents=True,mode=0o700); os.chmod(out,0o700)
    report={'module':a.module,'missing':missing,'source_sha256':source_hashes,
      'preflight_only':a.preflight,'date_status':'PROVISIONAL: saved dates retained, final raw endpoint policy not approved',
      'algorithm_change':'none; source/data/output path routing only',
      'third_party':'Original Oh/Goeminne inputs are read from controlled local source; not relicensed or redistributed',
      'scientific_final_release':False}
    (out/'entry_metadata.json').write_text(json.dumps(report,indent=2))
    if missing: raise SystemExit('Missing required inputs; see entry_metadata.json')
    if a.preflight: print('PREFLIGHT_DEPENDENCIES_EXIST',a.module); return
    if not a.allow_provisional: raise SystemExit('Saved endpoint definitions remain provisional; rerun with --allow-provisional for numerical reproduction only')
    for x in (code/source_module).iterdir():
        if x.is_file(): shutil.copy2(x,out/x.name)
    if a.module.startswith('organ_70_30'):
        for name in ['PRIVATE','numeric','qa','logs','panels']:
            (out/name).mkdir(mode=0o700,exist_ok=True)
    env=os.environ.copy(); env.update(DISCO_DATA_ROOT=str(data),DISCO_APPROVED_CODE_ROOT=str(code),
      OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    if questionnaire: env['DISCO_RULAS_QUESTIONNAIRE']=str(questionnaire.resolve(strict=True))
    if a.module=='organ_70_30_compare': env['DISCO_FROZEN_SCORE_SNAPSHOT']=str(frozen)
    else: env.pop('DISCO_FROZEN_SCORE_SNAPSHOT',None)
    for file,argv in MODULES[a.module]:
        params=[x.replace('@OUT@',str(out)).replace('@DATA@',str(data)) for x in argv]
        command=([a.python] if file.endswith('.py') else [a.rscript,'--vanilla'])+[str(out/file),*params]
        with (out/(file+'.launcher.log')).open('w') as log:
            r=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT)
        if r.returncode: raise SystemExit(f'{file} failed ({r.returncode}); preserve logs and use a new output for retry')
    report['completed_numeric_steps']=[x[0] for x in MODULES[a.module]]
    (out/'entry_metadata.json').write_text(json.dumps(report,indent=2))
    print('COMPLETED_CONDITIONAL_NUMERICAL_REPRODUCTION',a.module)

if __name__=='__main__': main()
