#!/usr/bin/env python3
"""Recover eight endpoint rows for organ IDs outside the original two-score union.

Matches the approved rebuild_dates.py rules, without inventing no-disease status.
Individual endpoint rows stay in the 0700 task directory; stdout is aggregate only.
"""
# File-system roots are explicit; no participant data are bundled.
import os as _disco_os
def _disco_path(kind, suffix=""):
    key = {"private": "DISCO_PRIVATE_ROOT", "results": "DISCO_RESULTS_ROOT", "legacy": "DISCO_LEGACY_ROOT"}[kind]
    root = _disco_os.environ.get(key)
    if not root:
        raise RuntimeError("Required environment variable: " + key)
    return str(__import__("pathlib").Path(root) / suffix) if suffix else root
import csv, datetime as dt, hashlib, io, json, os, sys, zipfile
from collections import Counter, defaultdict
from pathlib import Path
os.umask(0o077)
BASE=Path(_disco_path("private",""));OUT=Path(sys.argv[1]);CUT=dt.date(2022,12,31)
SPECIAL={'1900-01-01','1901-01-01','1902-02-02','1903-03-03','1909-09-09','2037-07-07'}
MAP={'obesity':['130792-0.0'],'hypertension':['131286-0.0','131294-0.0'],
 'diabetes':['130708-0.0'],'mi':['42000-0.0'],'stroke':['42006-0.0'],
 'copd':['42016-0.0'],'dementia':['42018-0.0'],'cancer':['40005-'+str(i)+'.0' for i in range(22)]}
def date(v):
 v=(v or '').strip()
 if v in ('','NA','nan','NaT'):return None,False
 if v[:10] in SPECIAL:return None,True
 try:return dt.date.fromisoformat(v[:10]),False
 except ValueError:return None,True
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def rows(z,n):return csv.DictReader(io.TextIOWrapper(z.open(n),encoding='utf-8-sig'))
org=BASE/'month_organ_execution_20261007/input_organ.csv'
endpoint=BASE/'month_disease_prs_execution_20261007/private_endpoint_8.csv'
before={str(org):sha(org),str(endpoint):sha(endpoint)}
with org.open() as f:oi={r['eid'] for r in csv.DictReader(f)}
with endpoint.open() as f:ei={r['eid'] for r in csv.DictReader(f) if r['disease']=='obesity'}
ids=oi-ei;assert len(ids)==52
input_presence={}
for mod in ('clinical','proteomics'):
 with (BASE/'month_disease_prs_execution_20261007'/('input_'+mod+'.csv')).open() as f:
  ix={r['eid'] for r in csv.DictReader(f)}
 input_presence[mod]=len(ids & ix)
assert input_presence=={'clinical':0,'proteomics':0}
blood={};raw={};cancer_types={};qa=Counter()
with zipfile.ZipFile(BASE/'UKB网站数据下载.zip') as z:
 tubes=next(n for n in z.namelist() if '_3166_' in n and n.endswith('.csv') and '__MACOSX' not in n)
 main=next(n for n in z.namelist() if '100001_' in n and n.endswith('.csv') and '__MACOSX' not in n)
 for r in rows(z,tubes):
  i=r['eid']
  if i not in ids:continue
  dd=[date(v)[0] for k,v in r.items() if k.startswith('3166-0.')];valid={d for d in dd if d}
  blood[i]=next(iter(valid)) if len(valid)==1 else None
  qa['blood_conflicting_dates']+=len(valid)>1
 for r in rows(z,main):
  i=r['eid']
  if i not in ids:continue
  assert i not in raw
  dd=[date(r.get('40000-'+str(j)+'.0')) for j in (0,1)]
  loss,lu=date(r.get('191-0.0'));death=min((x[0] for x in dd if x[0]),default=None)
  raw[i]=(loss,death,lu or any(x[1] for x in dd))
  cancer_types[i]={c.split('-',1)[1]:r[c] for c in r if c.startswith('40006-') and r[c]}
  visit,_=date(r.get('53-0.0'));b=blood.get(i)
  qa['blood_visit_dates_differ']+=bool(b and visit and b!=visit)
print('targeted main and tube lookup complete; aggregate only',flush=True)
ev8=defaultdict(dict);un8=defaultdict(set);seen8=defaultdict(set)
reverse={c:d for d,cs in MAP.items() for c in cs}
with zipfile.ZipFile(BASE/'新发疾病原始数据.zip') as z:
 for n in z.namelist():
  if not n.endswith('.csv') or '__MACOSX' in n:continue
  reader=rows(z,n);fields=[c for c in reader.fieldnames if c in reverse]
  if not fields:continue
  for r in reader:
   i=r['eid']
   if i not in ids:continue
   for c in fields:
    d=reverse[c];seen8[d].add(i);v=r[c];day,unknown=date(v)
    if unknown:un8[i].add(d)
    if day and (d not in ev8[i] or day<ev8[i][d]):ev8[i][d]=day
    if d=='cancer':
     typ=cancer_types.get(i,{}).get(c.split('-',1)[1],'')
     if typ and not day:un8[i].add(d)
     qa['cancer_date_without_type']+=bool(day and not typ)
fields=['eid','disease','days','time','event','tie','reason','after_death','after_loss','after_cutoff','unknown_date','prevalent','mortality_end','disease_end']
counts=[]
with (OUT/'private_organ_52_endpoint_8.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for d in MAP:
  cnt=Counter()
  for i in sorted(ids):
   b=blood.get(i);loss,death,bad=raw.get(i,(None,None,False));a=ev8.get(i,{}).get(d)
   end=min(x for x in [loss,death,CUT] if x);mort=death if death and death<=CUT else min(x for x in [loss,CUT] if x)
   prev=bool(a and b and a<=b);unk=d in un8.get(i,set());tie=bool(a and death and a==death and a<=end)
   reason=('missing_raw' if i not in raw else 'invalid_blood' if not b else 'invalid_end_date' if bad else
    'missing_disease_source' if i not in seen8[d] else 'prevalent' if prev else 'known_disease_unknown_date' if unk else 'nonpositive_followup' if end<=b else '')
   event=1 if a and a<=end and not tie else 2 if death and death==end else 0
   stop=a if event==1 else end;days=(stop-b).days if b else None
   if not reason:assert days>0
   row=dict(eid=i,disease=d,days=days,time=days/(365.25/12) if days is not None else None,event=event,tie=int(tie),reason=reason,
    after_death=int(bool(a and death and a>death)),after_loss=int(bool(a and loss and a>loss)),after_cutoff=int(bool(a and a>CUT)),
    unknown_date=int(unk),prevalent=int(prev),mortality_end=mort,disease_end=end)
   w.writerow(row);cnt['input_N']+=1;cnt[reason or 'eligible']+=1
   if not reason:cnt['events']+=event==1;cnt['competing_deaths']+=event==2;cnt['ties']+=tie
  counts.append(dict(disease=d,source_rows=len(seen8[d]),**cnt))
after={p:sha(p) for p in before};assert before==after
result=dict(status='EXECUTED',organ_N=len(oi),previous_overlap_N=len(oi&ei),outside_old_endpoint_union_N=len(ids),
 reason_for_omission='Original endpoints were restricted to union of clinical and pure-protein saved-score input populations; these 52 are in neither. Not an absence of disease evidence.',
 original_score_input_presence=input_presence,raw_main_found_N=len(raw),raw_blood_found_N=len(blood),counts=counts,qa=dict(qa),
 private_augmented_endpoint=str(OUT/'private_organ_52_endpoint_8.csv'),private_data_returned=False,
 protocol='Identical eight-disease field map and approved raw_dates death_first policy; administrative date not changed',
 source_hashes=before,sources_unchanged=True,rule_source=_disco_path("private","month_disease_prs_execution_20261007/rebuild_dates.py"))
(OUT/'organ_52_endpoint_aggregate.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
