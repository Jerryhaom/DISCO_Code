#!/usr/bin/env python3
"""Independent aggregate verification; does not refit or inspect individual data."""
import csv, hashlib, json, math, sys
from pathlib import Path
from pypdf import PdfReader

root=Path(__file__).resolve().parents[1]
read=lambda path:list(csv.DictReader(path.open(newline='')))
num=lambda r,k:float(r[k])
key=lambda r:(r['modality'],r['disease'],r['metric'],r['model'])
model_sources={key(r):r for r in read(root/'inputs/three_modalities_raw_unit_models.csv')}
adapt={'ClinicalProt':'clinicalprot','Prot':'proteomics','Clinical':'clinical'}
report={'independent_QA':True,'scope':'Frozen aggregate results and displayed numeric mappings, not a raw-person endpoint re-audit','checks':{}}
max_error=0
for model,prefix in [('extended_BMI','Fig3c_raw_unit_expanded_BMI_NEAR_SOURCE'),('crude','Fig3c_raw_unit_crude_NEAR_SOURCE')]:
    rows=read(root/f'qa/{prefix}_plotted.csv')
    assert len(rows)==48
    assert len({(r['modality'],r['disease'],r['metric']) for r in rows})==48
    for r in rows:
        s=model_sources[(adapt[r['modality']],r['disease'],r['metric'],model)]
        c=model_sources[(adapt[r['modality']],r['disease'],r['metric'],'crude')]
        for field in ['N','events','beta','se','HR','low','high','C','C_se']:
            err=abs(num(r,field)-num(s,field));max_error=max(err,max_error);assert err<1e-12
        assert abs(num(r,'C_crude')-num(c,'C'))<1e-12
        assert num(r,'N_crude')==num(c,'N') and num(r,'events_crude')==num(c,'events')
        assert abs(num(r,'percent_crude')-100*num(c,'events')/num(c,'N'))<1e-11
        assert abs(num(r,'HR')-math.exp(num(r,'beta')))<1e-12
        for field,sign in [('low',-1),('high',1)]:
            assert abs(num(r,field)-math.exp(num(r,'beta')+sign*1.959963984540054*num(r,'se')))<1e-11
        assert num(r,'low')<=num(r,'HR')<=num(r,'high')
    reader=PdfReader(root/f'output/pdf/{prefix}.pdf')
    text='\n'.join(p.extract_text() or '' for p in reader.pages)
    assert len(reader.pages)==1
    for r in rows:
        assert f"{num(r,'HR'):.2f}" in text
        assert f"{num(r,'C_crude'):.3f}" in text
    report['checks'][model]={'rows':48,'maximum_mapping_error':max_error,'raw_unit_scale':'saved DISCO and raw DM, not Table1 logDM','all_HR_and_C_rounded_text_present':True}

cif_source={(r['disease'],r['policy'],r['group'],r['time_months']):r for r in read(root/'inputs/ClinicalProt_CIF_aggregate.csv')}
cif=read(root/'qa/ClinicalProt_CIF_plotted.csv')
names={'Obesity':'obesity','Hypertension':'hypertension','Type 2 diabetes':'diabetes','Cancer':'cancer','MI':'mi','Stroke':'stroke','COPD':'copd','Dementia':'dementia'}
assert len(cif)==480
by={}
for r in cif:
    s=cif_source[(names[r['disease']],r['policy'],r['group'],r['time_months'])]
    for field in ['CIF','CIF_low','CIF_high','death_CIF','KM','N','at_risk']:
        assert abs(num(r,field)-num(s,field))<1e-14
    assert abs(num(r,'years')-num(r,'time_months')/12)<1e-14
    assert 0<=num(r,'CIF_low')<=num(r,'CIF')<=num(r,'CIF_high')<=1
    assert num(r,'CIF')<=num(r,'KM')+1e-12
    assert num(r,'CIF')+num(r,'death_CIF')<=1+1e-12
    assert num(r,'at_risk')>=num(r,'N')*.1
    by.setdefault((r['disease'],r['group']),[]).append(r)
for group in by.values():
    group.sort(key=lambda r:num(r,'years'))
    assert len(group)==15 and num(group[0],'years')==0 and num(group[-1],'years')==14
    for a,b in zip(group,group[1:]):
        assert num(b,'CIF')>=num(a,'CIF')-1e-12
        assert num(b,'KM')>=num(a,'KM')-1e-12
        assert num(b,'at_risk')<=num(a,'at_risk')
qa=read(root/'inputs/ClinicalProt_CIF_QA.csv')
assert len(qa)==80 and max(num(r,'max_AJ_error') for r in qa)<1e-12
for r in qa:
    if r['policy']=='death_first' and r['group']=='All':
        for metric in ('DISCO','DM'):
            s=model_sources[('clinicalprot',r['disease'],metric,'crude')]
            assert num(r,'N')==num(s,'N') and num(r,'events')==num(s,'events')
reader=PdfReader(root/'output/pdf/ClinicalProt_CIF_sensitivity.pdf')
assert len(reader.pages)==2
for page,metric in zip(reader.pages,('DISCO','DM')):
    text=page.extract_text() or ''
    assert f'{metric} Q1' in text and f'{metric} Q4' in text
    assert all(title in text for title in names)
report['checks']['CIF']={'plotted_rows':480,'strata':32,'pages':2,'grid_years':list(range(15)),
    'minimum_at_risk_fraction':min(num(r,'at_risk')/num(r,'N') for r in cif),
    'max_precomputed_independent_AJ_error':max(num(r,'max_AJ_error') for r in qa),
    'all_counts_match_crude_ClinicalProt_models':True,
    'annual_interpolation_not_event_time_curve':True,'pointwise_not_simultaneous_CI':True}
report['result']='PASS'
json.dump(report,sys.stdout,indent=2,ensure_ascii=False)
print()
