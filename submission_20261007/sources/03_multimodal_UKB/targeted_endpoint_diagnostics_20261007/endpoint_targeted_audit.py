#!/usr/bin/env python3
"""Server-only read-only participant audit; public stdout contains aggregates only."""
# File-system roots are explicit; no participant data are bundled.
import os as _disco_os
def _disco_path(kind, suffix=""):
    key = {"private": "DISCO_PRIVATE_ROOT", "results": "DISCO_RESULTS_ROOT", "legacy": "DISCO_LEGACY_ROOT"}[kind]
    root = _disco_os.environ.get(key)
    if not root:
        raise RuntimeError("Required environment variable: " + key)
    return str(__import__("pathlib").Path(root) / suffix) if suffix else root
import collections
import csv
import datetime as dt
import hashlib
import io
import json
import os
import re
import sys
import zipfile
from pathlib import Path

BASE = Path(_disco_path("private",""))
MORT = BASE / 'month_endpoint_execution_20261007_0WHBMrmO'
DIS = BASE / 'month_disease_prs_execution_20261007'
OUT = Path(sys.argv[1])
os.umask(0o077)

def date(s):
    return dt.date.fromisoformat(s[:10]) if s and s not in ('NA', 'NaN', 'nan') else None

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def bucket(x):
    if x < -365: return 'earlier_more_than_365d'
    if x < -30: return 'earlier_31_to_365d'
    if x < -1: return 'earlier_2_to_30d'
    if x == -1: return 'earlier_1d'
    if x == 0: return 'same_day'
    if x == 1: return 'later_1d'
    if x <= 30: return 'later_2_to_30d'
    if x <= 365: return 'later_31_to_365d'
    return 'later_more_than_365d'

counts = collections.Counter()
offsets = []
abnormal = {}
loss_death = {}
before = {}
for p in [MORT / 'private_endpoint.csv', DIS / 'input_proteomics.csv']:
    before[str(p)] = sha(p)
with (MORT / 'private_endpoint.csv').open() as f:
    for row in csv.DictReader(f):
        counts['mortality_input_n'] += 1
        if row['reason']: continue
        counts['eligible_n'] += 1
        b, v = date(row['blood']), date(row['visit'])
        death, loss = date(row['death']), date(row['loss'])
        delta = (b-v).days
        counts['baseline_offset_' + bucket(delta)] += 1
        if delta:
            abnormal[row['eid']] = row
            offsets.append(delta)
            counts['anomalous_date_before_first_2006_assessment'] += b < dt.date(2006,3,13)
            counts['anomalous_date_after_last_2010_assessment'] += b > dt.date(2010,10,1)
            counts['anomalous_deaths'] += int(row['status'])
            counts['anomalous_loss_censored'] += int(row['reason'] == '' and row['status'] == '0' and bool(loss))
        if death and loss and loss < death:
            loss_death[row['eid']] = row
            counts['registered_death_after_loss_n'] += 1

prot = set()
with (DIS / 'input_proteomics.csv').open() as f:
    for row in csv.DictReader(f): prot.add(row['eid'])
counts['proteomics_input_n'] = len(prot)
counts['table1_anomalies_in_proteomics_input'] = len(set(abnormal) & prot)
counts['loss_death_in_proteomics_input'] = len(set(loss_death) & prot)

source_zip = BASE / 'UKB网站数据下载.zip'
metadata = []
header_presence = {}
tube_checks = collections.Counter()
with zipfile.ZipFile(source_zip) as z:
    selected = [n for n in z.namelist() if ('100001_' in n or '_3166_' in n) and n.endswith('.csv') and '__MACOSX' not in n]
    for n in selected:
        info = z.getinfo(n)
        with z.open(n) as stream:
            rd = csv.DictReader(io.TextIOWrapper(stream, encoding='utf-8-sig'))
            fields = rd.fieldnames
            metadata.append(dict(member=n, uncompressed_bytes=info.file_size, CRC=info.CRC,
                                 zip_member_time=list(info.date_time), field_count=len(fields)))
            if '100001_' in n:
                wanted = ['53-0.0','54-0.0','190-0.0','191-0.0','21842-0.0','40000-0.0','40000-1.0',
                          '40020-0.0','40021-0.0']
                header_presence.update({k: k in fields for k in wanted})
                header_presence['sample_time_fields'] = [c for c in fields if c.startswith(('21842-', '3166-'))]
                # Field IDs listed here are a header inventory only; no field is
                # asserted to establish provider/geographical coverage.
                # Header only: no additional scan of the 6 GB main table.
            else:
                blood_fields = [c for c in fields if c.startswith('3166-')]
                inst = collections.Counter(c.split('-')[1].split('.')[0] for c in blood_fields)
                header_presence['3166_instance_array_counts'] = dict(inst)
                for r in rd:
                    if r['eid'] not in abnormal: continue
                    b = date(abnormal[r['eid']]['blood'])
                    v = date(abnormal[r['eid']]['visit'])
                    dd = [date(r[c]) for c in blood_fields if c.startswith('3166-0.') and r[c] and r[c][:10] != '1900-01-01']
                    tube_checks['anomaly_matching_raw3166_n'] += 1
                    tube_checks['all_valid_baseline_tubes_same_selected_date'] += bool(dd and all(x == b for x in dd))
                    tube_checks['any_baseline_tube_equals_visit'] += v in dd
                    tube_checks['one_or_more_not_performed_tubes'] += any(r[c][:10] == '1900-01-01' for c in blood_fields if c.startswith('3166-0.'))
                    tube_checks['repeat_instance_nonmissing_values'] += sum(bool(r[c]) and r[c][:10] != '1900-01-01' for c in blood_fields if not c.startswith('3166-0.'))
    ancillary_members = [n for n in z.namelist() if '__MACOSX' not in n and re.search(r'metadata|readme|manifest|\.html?$|\.json$|\.xml$|\.log$', n, re.I)]

# Existing private disease endpoint: collect only six loss/death subjects and ovarian
# source diagnostics. Individual records stay in RAM and are never returned.
window_counts = collections.Counter()
ovarian_counts = collections.Counter()
with (DIS / 'private_endpoint_28.csv').open() as f:
    for row in csv.DictReader(f):
        if row['eid'] in loss_death:
            window_counts['loss_death_pairs_total'] += 1
            window_counts['eligible_pairs'] += not bool(row['reason'])
            window_counts['eligible_censor_at_loss_pairs'] += not row['reason'] and row['event'] == '0'
            window_counts['eligible_competing_death_pairs'] += not row['reason'] and row['event'] == '2'
            window_counts['diagnosis_after_loss_pairs'] += int(row.get('after_loss','0'))
            window_counts['mortality_end_differs_disease_end_pairs'] += row['mortality_end'] != row['disease_end']
        if row['disease'] == 'ovarian_ca' and row['eid'] in prot:
            ovarian_counts['proteomics_endpoint_rows'] += 1
            ovarian_counts['eligible_before_sex_covariate_prs_filter'] += not bool(row['reason'])

# Metadata-file inventory, no arbitrary log contents or patient rows emitted.
inventory = []
pruned = {'node_modules', '.git', '__pycache__', 'DNAmAge', 'UKB_hzx', 'bioage_v2',
          'month_disease_prs_execution_20261007', 'month_endpoint_execution_20261007_0WHBMrmO'}
for root, dirs, files in os.walk(BASE):
    dirs[:] = [d for d in dirs if d not in pruned and not d.startswith('targeted_endpoint_diagnostics_')]
    for name in files:
        if re.search(r'metadata|manifest|provider|censor|coverage|release|readme|download', name, re.I) and not name.endswith(('.csv','.rds','.RData','.npy','.npz','.parquet')):
            p = Path(root) / name
            try: size=p.stat().st_size
            except OSError: continue
            inventory.append(dict(path=str(p),bytes=size))

after = {p:sha(Path(p)) for p in before}
assert before == after
result = dict(scope='targeted read-only aggregate; no endpoint or model replaced',
              counts=dict(counts), anomaly_offset_days=dict(min=min(offsets),max=max(offsets),median=sorted(offsets)[len(offsets)//2]),
              raw_tube_checks=dict(tube_checks), raw_export_metadata=metadata,
              field_presence=header_presence, ancillary_zip_members=ancillary_members,
              loss_death_disease_window=dict(window_counts), ovarian_input_counts=dict(ovarian_counts),
              metadata_candidates=inventory, source_hashes_before=before, source_hashes_after=after,
              original_sources_unchanged=True,
              limitations=['ZIP timestamps/export filename are export provenance, not provider completeness evidence',
                           'Metadata inventory excludes DNAmAge, UKB_hzx and bioage_v2 trees; absence is not exhaustive server proof',
                           'Unique tube dates do not prove anomalous clock dates biologically correct',
                           'Field54 assessment centre is not treated as actual provider/geography',
                           'No disease date after loss is asserted observed merely because a later death was registered'])
(OUT / 'endpoint_targeted_aggregate.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps(result,ensure_ascii=False,indent=2))
