"""Inspect public ED1 source excerpts or a supplied aggregate/synthetic export manifest."""
from pathlib import Path
import hashlib
import json
import os
import csv

ROOT=Path(__file__).resolve().parents[1]
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
contract=json.loads((ROOT/'original_source_contract.json').read_text())
source_records=[]
for block in contract['blocks']:
 path=ROOT/'source'/Path(block['member']).name
 assert digest(path)==block['sha256']
 source_records.append(dict(source='source/'+path.name,first_line=block['first_line'],last_line=block['last_line'],sha256=digest(path)))
assert 'set.seed(123)' in (ROOT/'source/ED1_function_681_788.R').read_text()
result=dict(status='PASS_PUBLIC_ED1_SOURCE_CONTRACT',source_records=source_records,seed=123,observed_generator_executed=False,participant_data_read=False,whole_original_source_required_from_caller=True)
manifest_path=os.environ.get('DISCO_ED1_MANIFEST')
if manifest_path:
 manifest=Path(manifest_path).resolve();rows=json.loads(manifest.read_text())
 data_root=manifest.parents[1]
 assert rows and all(row['privacy']=='aggregate_public_candidate' for row in rows)
 checked=[]
 for row in rows:
  logical=row['output'];prefix='caller:DISCO_RUN_MODULE/'
  assert logical.startswith(prefix)
  rel=Path(logical[len(prefix):]);assert not rel.is_absolute() and '..' not in rel.parts
  path=data_root/rel;assert digest(path)==row['output_sha256']
  with path.open(newline='') as f:
   data=list(csv.DictReader(f))
  assert len(data)==row['rows']
  assert not set(row['columns'])&{'ID','id','id1','eid','seqn','participant_id','death_date','diagnosis_date'}
  checked.append(dict(panel=row['panel'],component=row['component'],rows=len(data),sha256=digest(path)))
 result.update(status='PASS_CALLER_AGGREGATE_SYNTHETIC_MANIFEST',manifest_sha256=digest(manifest),checked_components=checked,scope='Verify supplied ED1 synthetic/aggregate files and logical manifest; not rerun original generator or any cohort model')
(ROOT/'results/ED1_public_contract_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('source_records','checked_components')},indent=2))
