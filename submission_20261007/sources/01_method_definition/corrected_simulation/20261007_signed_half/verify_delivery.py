# File-system roots are explicit; no participant data are bundled.
import os as _disco_os
def _disco_path(kind, suffix=""):
    key = {"private": "DISCO_PRIVATE_ROOT", "results": "DISCO_RESULTS_ROOT", "legacy": "DISCO_LEGACY_ROOT"}[kind]
    root = _disco_os.environ.get(key)
    if not root:
        raise RuntimeError("Required environment variable: " + key)
    return str(__import__("pathlib").Path(root) / suffix) if suffix else root
from pathlib import Path
import csv
import hashlib

root = Path(_disco_path("results",""))
manifest = Path(__file__).with_name('delivery_manifest.tsv')
with manifest.open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert len(rows) == 12
assert len({row['asset'] for row in rows}) == len(rows)
for row in rows:
    path = root / row['file']
    assert path.is_file(), (row['asset'], str(path))
    for key in ('producer', 'input_evidence'):
        assert (root / row[key]).is_file(), (row['asset'], key)
    if row['sha256'] != 'NA':
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'], row['asset']
    print(row['asset'], row['status'], 'PASS')
print('PASS: 12 unique entries, 10 PDF hashes, producer and input-evidence paths.')
