#!/usr/bin/env python3
# File-system roots are explicit; no participant data are bundled.
import os as _disco_os
def _disco_path(kind, suffix=""):
    key = {"private": "DISCO_PRIVATE_ROOT", "results": "DISCO_RESULTS_ROOT", "legacy": "DISCO_LEGACY_ROOT"}[kind]
    root = _disco_os.environ.get(key)
    if not root:
        raise RuntimeError("Required environment variable: " + key)
    return str(__import__("pathlib").Path(root) / suffix) if suffix else root
import csv,datetime as dt,hashlib,io,json,os,collections,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
os.umask(0o077)
ZIP=Path(_disco_path("private","UKB\u7f51\u7ad9\u6570\u636e\u4e0b\u8f7d.zip"))
CUTOFF=dt.date(2022,12,31)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
def date(v):
    v=(v or '').strip()
    if v in ('','NA','NaN','nan','NULL','None'): return None
    try: return dt.date.fromisoformat(v[:10])
    except ValueError: raise ValueError('Unparsed nonmissing source date; no participant values emitted')
def read(z,member,fields):
    data={}
    with z.open(member) as f:
        r=csv.DictReader(io.TextIOWrapper(f,encoding='utf-8-sig'))
        assert all(x in r.fieldnames for x in ['eid']+fields)
        for row in r:
            eid=str(int(row['eid']))
            assert eid not in data,'Duplicate raw ID'
            data[eid]={c:date(row[c]) for c in fields}
    return data
before={str(p):sha(p) for p in [ZIP,ROOT/'analysis_input.csv']}
with zipfile.ZipFile(ZIP) as z:
    main=[n for n in z.namelist() if '100001_' in n and n.endswith('.csv') and '__MACOSX' not in n]
    blood=[n for n in z.namelist() if '_3166_' in n and n.endswith('.csv') and '__MACOSX' not in n]
    assert len(main)==len(blood)==1
    raw=read(z,main[0],['53-0.0','191-0.0','40000-0.0','40000-1.0'])
    tubes=read(z,blood[0],['3166-0.'+str(i) for i in range(7)])
counts=collections.Counter(); reasons=collections.Counter()
fields=['eid','visit','blood','loss','death0','death1','death','end','days','time','status','reason','blood_visit_different','tube_not_performed']
with open(ROOT/'analysis_input.csv') as fi,open(ROOT/'private_endpoint.csv','w') as fo:
    w=csv.DictWriter(fo,fieldnames=fields);w.writeheader()
    seen=set()
    for r in csv.DictReader(fi):
        i=r['eid'];assert i not in seen;seen.add(i)
        a=raw.get(i,{});t=tubes.get(i,{})
        valid=set(v for v in t.values() if v and v!=dt.date(1900,1,1))
        b=next(iter(valid)) if len(valid)==1 else None
        visit=a.get('53-0.0');loss=a.get('191-0.0');ds=[a.get('40000-'+str(j)+'.0') for j in (0,1)]
        death=min((d for d in ds if d),default=None)
        reason=('missing_raw_ID' if i not in raw else 'missing_blood_ID' if i not in tubes else
                'inconsistent_tube_dates' if len(valid)>1 else 'no_valid_blood_date' if not b else '')
        end=None;status=None;days=None
        if not reason:
            if death and death<=b: reason='baseline_or_prebaseline_death'
            elif death and death<=CUTOFF: end=death;status=1
            else: end=min(loss,CUTOFF) if loss else CUTOFF;status=0
            if end:
                days=(end-b).days
                if days<=0: reason='nonpositive_followup'
        different=bool(b and visit and b!=visit)
        notperformed=any(v==dt.date(1900,1,1) for v in t.values())
        counts['input_n']+=1;counts['blood_visit_different']+=different;counts['has_not_performed_tube']+=notperformed
        counts['death_dates_disagree']+=bool(all(ds) and ds[0]!=ds[1])
        if not reason:
            counts['eligible_n']+=1;counts['events']+=status
            counts['death_overrides_earlier_loss']+=bool(status and loss and loss<end)
            counts['future_death_loss_censored']+=bool(death and death>CUTOFF and loss and loss<=CUTOFF)
            counts['loss_censored']+=bool(not status and loss and loss<=CUTOFF)
        reasons[reason or 'eligible']+=1
        w.writerow(dict(eid=i,visit=visit or '',blood=b or '',loss=loss or '',death0=ds[0] or '',death1=ds[1] or '',
            death=death or '',end=end or '',days=days if days is not None else '',time=days/(365.25/12) if days is not None else '',
            status=status if status is not None else '',reason=reason,blood_visit_different=int(different),tube_not_performed=int(notperformed)))
after={p:sha(p) for p in before};assert before==after
json.dump(dict(counts=counts,exclusive_reasons=reasons,source_hashes_before=before,source_hashes_after=after,
 source_unchanged=True,main_member=main[0],blood_member=blood[0],cutoff=str(CUTOFF),
 time_unit='fractional months: raw days/(365.25/12)',private_endpoint=str(ROOT/'private_endpoint.csv'),
 endpoint_sha256=sha(ROOT/'private_endpoint.csv')),open(ROOT/'endpoint_aggregate.json','w'),indent=2)
print(json.dumps(dict(counts=counts,exclusive_reasons=reasons)))
