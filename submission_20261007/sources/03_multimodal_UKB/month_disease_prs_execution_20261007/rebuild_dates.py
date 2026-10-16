#!/usr/bin/env python3
"""Raw dates only; all participant records remain in this 0700 server directory."""
# File-system roots are explicit; no participant data are bundled.
import os as _disco_os
def _disco_path(kind, suffix=""):
    key = {"private": "DISCO_PRIVATE_ROOT", "results": "DISCO_RESULTS_ROOT", "legacy": "DISCO_LEGACY_ROOT"}[kind]
    root = _disco_os.environ.get(key)
    if not root:
        raise RuntimeError("Required environment variable: " + key)
    return str(__import__("pathlib").Path(root) / suffix) if suffix else root
import ast,csv,datetime as dt,hashlib,io,json,os,re,zipfile
from pathlib import Path
from collections import Counter,defaultdict
os.umask(0o077)
ROOT=Path(__file__).resolve().parent
BASE=Path(_disco_path("private","")); CUT=dt.date(2022,12,31)
SPECIAL={'1900-01-01','1901-01-01','1902-02-02','1903-03-03','1909-09-09','2037-07-07'}
MAP={'obesity':['130792-0.0'],'hypertension':['131286-0.0','131294-0.0'],
 'diabetes':['130708-0.0'],'mi':['42000-0.0'],'stroke':['42006-0.0'],
 'copd':['42016-0.0'],'dementia':['42018-0.0'],'cancer':['40005-'+str(i)+'.0' for i in range(22)]}
qa=Counter(); hashes={}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def date(v):
 v=(v or '').strip()
 if v in ('','NA','nan','NaT'):return None,False
 if v[:10] in SPECIAL:return None,True
 try:return dt.date.fromisoformat(v[:10]),False
 except ValueError:return None,True
def rows(z,n):return csv.DictReader(io.TextIOWrapper(z.open(n),encoding='utf-8-sig'))
ids=set()
for mod in ('clinical','proteomics'):
 p=ROOT/('input_'+mod+'.csv');hashes[str(p)]=sha(p)
 with open(p) as f:ids.update(r['eid'] for r in csv.DictReader(f))
source=BASE/'bioage_v2/scripts/rebuild_diseases_v2.py';hashes[str(source)]=sha(source)
tree=ast.parse(source.read_text()); rules=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='ICD_RULES' for t in n.targets))
sr=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SR_CODE' for t in n.targets))
lookup={}
for ch in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
 for num in range(100):lookup[ch+f'{num:02d}']=[d for d,rr in rules.items() if any(ch==c and lo<=num<=hi for c,lo,hi in rr)]
raw={};blood={};events28=defaultdict(dict);unknown28=defaultdict(set);prev28=defaultdict(set);cancer_types={}
def add(i,ds,v):
 day,unknown=date(v)
 for d in ds:
  if unknown or day is None:unknown28[i].add(d)
  elif d not in events28[i] or day<events28[i][d]:events28[i][d]=day
zpath=BASE/'UKB网站数据下载.zip';hashes[str(zpath)]=sha(zpath)
with zipfile.ZipFile(zpath) as z:
 main=next(n for n in z.namelist() if '100001_' in n and n.endswith('.csv') and '__MACOSX' not in n)
 tubes=next(n for n in z.namelist() if '_3166_' in n and n.endswith('.csv') and '__MACOSX' not in n)
 for r in rows(z,tubes):
  i=r['eid']
  if i not in ids:continue
  dd=[date(v)[0] for k,v in r.items() if k.startswith('3166-0.')];valid=set(d for d in dd if d)
  blood[i]=next(iter(valid)) if len(valid)==1 else None
  qa['blood_conflicting_dates']+=len(valid)>1
 reader=rows(z,main);header=reader.fieldnames
 hc=[c for c in header if c.startswith('41270-')];pairs=[(c,c.replace('41270-','41280-')) for c in hc]
 assert len(pairs)==259 and all(d in header for c,d in pairs)
 dc=[c for c in header if c.startswith(('40001-','40002-'))]
 sc=[c for c in header if c.startswith('20002-')]
 cc=[c for c in header if c.startswith('40006-')]
 for r in reader:
  i=r['eid']
  if i not in ids:continue
  assert i not in raw
  dd=[date(r.get('40000-'+str(j)+'.0')) for j in (0,1)]
  loss,lu=date(r.get('191-0.0'));death=min((x[0] for x in dd if x[0]),default=None)
  raw[i]=(loss,death,lu or any(x[1] for x in dd))
  cancer_types[i]={c.split('-',1)[1]:r[c] for c in cc if r[c]}
  for c,d in pairs:
   code=r[c].strip().upper();ds=lookup.get(code[:3],[])
   if ds:add(i,ds,r[d]);qa['HES_classified_pairs']+=bool(ds)
  for c in dc:
   ds=lookup.get(r[c].strip().upper()[:3],[])
   if ds:add(i,ds,r.get('40000-'+c.split('-')[1].split('.')[0]+'.0',''))
  # Preserve self-report definitions, but never invent a date from integer age.
  # At/before blood visit establishes prevalence; later reports have unknown onset.
  for c in sc:
   d=sr.get(r[c].strip())
   if not d:continue
   instance=c.split('-')[1].split('.')[0]; visit,_=date(r.get('53-'+instance+'.0'))
   b=blood.get(i)
   if visit and b and visit<=b:prev28[i].add(d)
   else:unknown28[i].add(d);qa['self_report_no_exact_onset']+=1
  if len(raw)%100000==0:print('raw people',len(raw),flush=True)
print('full raw dates read',len(raw),flush=True)
ev8=defaultdict(dict);un8=defaultdict(set);seen8=defaultdict(set)
zp=BASE/'新发疾病原始数据.zip';hashes[str(zp)]=sha(zp)
reverse={c:d for d,cs in MAP.items() for c in cs}
with open(ROOT/'private_cancer_pairs.csv','w') as f:
 w=csv.writer(f);w.writerow(['eid','instance','raw_date','icd10_type','special_or_invalid'])
 with zipfile.ZipFile(zp) as z:
  for n in z.namelist():
   if not n.endswith('.csv') or '__MACOSX' in n:continue
   reader=rows(z,n);fields=[c for c in reader.fieldnames if c in reverse]
   for r in reader:
    i=r['eid']
    if i not in ids:continue
    for c in fields:
     d=reverse[c];seen8[d].add(i);v=r[c];day,unknown=date(v)
     if unknown:un8[i].add(d)
     if day and (d not in ev8[i] or day<ev8[i][d]):ev8[i][d]=day
     if d=='cancer':
      typ=cancer_types.get(i,{}).get(c.split('-',1)[1],'')
      if v or typ:
       w.writerow([i,c.split('-',1)[1],v,typ,int(unknown)])
       if typ and not day:un8[i].add(d)
       qa['cancer_date_without_type']+=bool(day and not typ)
def output(tag,diseases,events,unknown,prevalent):
 counts=[]
 with open(ROOT/('private_endpoint_'+tag+'.csv'),'w') as f:
  fields=['eid','disease','days','time','event','tie','reason','after_death','after_loss','after_cutoff','unknown_date','prevalent','mortality_end','disease_end']
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for d in diseases:
   cnt=Counter()
   for i in sorted(ids):
    b=blood.get(i);loss,death,bad=raw.get(i,(None,None,False));a=events.get(i,{}).get(d)
    end=min(x for x in [loss,death,CUT] if x);mort=death if death and death<=CUT else min(x for x in [loss,CUT] if x)
    prev=bool((a and b and a<=b) or d in prevalent.get(i,set()))
    unk=d in unknown.get(i,set());tie=bool(a and death and a==death and a<=end)
    reason=('missing_raw' if i not in raw else 'invalid_blood' if not b else 'invalid_end_date' if bad else
      'missing_disease_source' if tag=='8' and i not in seen8[d] else 'prevalent' if prev else 'known_disease_unknown_date' if unk else 'nonpositive_followup' if end<=b else '')
    event=1 if a and a<=end and not tie else 2 if death and death==end else 0
    stop=a if event==1 else end;days=(stop-b).days if b else None
    if not reason:assert days>0
    row=dict(eid=i,disease=d,days=days,time=days/(365.25/12) if days is not None else None,event=event,tie=int(tie),reason=reason,
      after_death=int(bool(a and death and a>death)),after_loss=int(bool(a and loss and a>loss)),after_cutoff=int(bool(a and a>CUT)),
      unknown_date=int(unk),prevalent=int(prev),mortality_end=mort,disease_end=end)
    w.writerow(row);cnt['input_N']+=1;cnt[reason or 'eligible']+=1
    if not reason:cnt['events']+=event==1;cnt['competing_deaths']+=event==2;cnt['ties']+=tie
    for k in ['after_death','after_loss','after_cutoff','unknown_date']:cnt[k]+=row[k]
   counts.append(dict(disease=d,**cnt))
 return counts
q8=output('8',MAP,ev8,un8,{})
json.dump(dict(counts=q8,qa=dict(qa),field_map=MAP,source_hashes=hashes),open(ROOT/'endpoint8_QA.json','w'),indent=2)
print('8 endpoints ready',flush=True)
q28=output('28',rules,events28,unknown28,prev28)
json.dump(dict(counts=q28,qa=dict(qa),rules=rules,self_report=sr,HES_pairs=pairs,source_hashes=hashes,
 unknown_policy='Exclude known disease with unresolved onset; self-report not converted to dates',
 disease_end='min(loss,death,2022-12-31)',mortality_end='death through cutoff overrides loss; otherwise min(loss,cutoff)',
 tie_policy='death-first analysis convention; disease-first sensitivity retained'),open(ROOT/'endpoint28_QA.json','w'),indent=2)
assert all(sha(Path(p))==h for p,h in hashes.items())
print('28 endpoints ready; source hashes unchanged',flush=True)
