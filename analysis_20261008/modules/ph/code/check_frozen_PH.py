import os
#!/usr/bin/env python3
"""Audit already-fitted models; no Cox model, score or participant is changed."""
import csv, hashlib, json, math, shutil
from pathlib import Path
import numpy as np

BASE = Path(os.environ['DISCO_RESULTS_INPUT'])
MOD = Path(os.environ['DISCO_RUN_MODULE'])
UP = MOD/'frozen_aggregate_input'
for d in ('input','code','results','logs','figures','tables'): (MOD/d).mkdir(parents=True,exist_ok=True)
sources = [UP/x for x in ('training_PH_both_bases.csv','training_HR_both_scales.csv',
 'training_calibration.csv','exploratory_5yr_piecewise_HR.csv','execution_QA.json',
 'paired_Cindex_differences.csv')]
sources += [UP/'02_fit_clock_endpoint.R',
 UP/'acceptance.json']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
before = {str(p):sha(p) for p in sources}
for p in sources: shutil.copy2(p,MOD/'input'/p.name)
def rows(name):
 with (UP/name).open(newline='') as f: return list(csv.DictReader(f))
ph, hr, cal, pw = [rows(x) for x in ('training_PH_both_bases.csv','training_HR_both_scales.csv',
 'training_calibration.csv','exploratory_5yr_piecewise_HR.csv')]
qa=json.loads((UP/'execution_QA.json').read_text())
assert len(ph)==180 and len(hr)==144 and len(pw)==72
assert all(int(r['N'])==31381 and int(r['deaths'])==3318 for r in ph+hr)
assert all(int(r['N'])==31381 for r in pw)
def count(basis,terms,col='P'):
 a=[r for r in ph if r['basis']==basis and r['term'] in terms]
 return dict(tested=len(a),nominalP05=sum(float(r[col])<.05 for r in a),
  minP=min(float(r[col]) for r in a),maxP=max(float(r[col]) for r in a))
global_err=[]
for r in ph:
 if r['basis']=='age_residual_SD' and r['term']=='GLOBAL':
  z=next(q for q in ph if q['basis']=='raw_score_SD' and q['term']=='GLOBAL' and
         q['organ']==r['organ'] and q['method']==r['method'])
  global_err.append(abs(float(r['chisq'])-float(z['chisq'])))
coeff_err=[];se_err=[];p_relative=[]
for r in hr:
 if r['scale']!='age_residual_SD': continue
 sc={'Da':'D','Oa':'O','Ga':'G'}[r['term']]
 c=next(q for q in cal if q['organ']==r['organ'] and q['score']==sc)
 z=next(q for q in hr if q['scale']=='raw_score_SD' and q['organ']==r['organ'] and
 q['method']==r['method'] and q['model']==r['model'] and q['term']=={'Da':'Dz','Oa':'Oz','Ga':'Gz'}[r['term']])
 ratio=float(c['age_residual_SD'])/float(c['raw_SD'])
 coeff_err.append(abs(float(r['beta'])-float(z['beta'])*ratio))
 se_err.append(abs(float(r['SE'])-float(z['SE'])*ratio))
 p_relative.append(abs(float(r['P'])-float(z['P']))/max(float(r['P']),float(z['P']),1e-300))
assert max(global_err)<1e-8 and max(coeff_err)<1e-10 and max(se_err)<1e-10
pwD=[r for r in pw if r['term']=='Da']
assert len(pwD)==36
ints=[r for r in pwD if r['period']=='0-60 months']
assert len(ints)==18 and all('post-hoc exploratory' in r['status'] for r in pw)
ordered=sorted((float(r['interaction_P']),i) for i,r in enumerate(ints));bh=[0.0]*18
acc=1.0
for rank in range(18,0,-1):
 p,i=ordered[rank-1];acc=min(acc,p*18/rank);bh[i]=acc
detail=[]
for i,r in enumerate(ints):
 late=next(q for q in pwD if q['organ']==r['organ'] and q['method']==r['method'] and q['period']=='>60 months')
 detail.append(dict(organ=r['organ'],method=r['method'],earlyHR=float(r['HR']),earlyLower=float(r['lower95']),
  lateHR=float(late['HR']),lateLower=float(late['lower95']),interactionP=float(r['interaction_P']),
  interactionBH=bh[i],earlyEvents=int(r['events']),lateEvents=int(late['events'])))
jd=[r for r in hr if r['model']=='joint' and r['term']=='Da' and r['scale']=='age_residual_SD']
summary=dict(status='PASS',modelRefits=0,participantsChanged=0,PHTestsAdded=0,
 trainingN=31381,trainingEvents=3318,heldoutN=13450,heldoutEvents=1423,
 residualDISCO=count('age_residual_SD',{'Da'}),rawDISCO=count('raw_score_SD',{'Dz'}),
 residualGLOBAL=count('age_residual_SD',{'GLOBAL'}),rawGLOBAL=count('raw_score_SD',{'GLOBAL'}),
 maxGlobalChisqBasisDifference=max(global_err),maxBetaRescalingError=max(coeff_err),
 maxSERescalingError=max(se_err),maxPRelativeDifference=p_relative and max(p_relative),
 jointResidualDISCOHRrange=[min(float(r['HR']) for r in jd),max(float(r['HR']) for r in jd)],
 all18JointResidualDISCOLowerCIabove1=all(float(r['lower95'])>1 for r in jd),
 posthoc60months=dict(DISCOMethodPairs=18,periodEstimates=36,allHRabove1=all(float(r['HR'])>1 for r in pwD),
  allLowerCIabove1=all(float(r['lower95'])>1 for r in pwD),
  earlyHRrange=[min(float(r['HR']) for r in ints),max(float(r['HR']) for r in ints)],
  lateHRrange=[min(x['lateHR'] for x in detail),max(x['lateHR'] for x in detail)],
  nominalInteractionP05=sum(x['interactionP']<.05 for x in detail),BHInteractionQ05=sum(x['interactionBH']<.05 for x in detail),
  earlyEvents=700,lateEvents=2618,details=detail,
  limitation='Post-hoc training sensitivity, not a new primary estimate, not a PH-pass claim. Baseline N is not late-period risk-set N.'),
 basisInterpretation='Same full constant-coefficient design space. Individual cox.zph alternatives depend on residual/raw basis; GLOBAL is invariant to numerical precision.',
 HRInterpretation='Working-model association summary over the observed follow-up, not a constant or causal effect and not an arithmetic average of period HRs.',
 sparseAlzheimer=json.loads(sources[-1].read_text()),upstreamExecutionQA=qa,sourceSHA256=before)
assert before=={str(p):sha(p) for p in sources}
(MOD/'results/PH_interpretation_checks.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
(MOD/'logs/runtime.txt').write_text('Python NumPy '+np.__version__+'\nNo models refitted. Sources unchanged.\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ('sparseAlzheimer','sourceSHA256','upstreamExecutionQA')},indent=2))
