from project_paths import CODE_ROOT, DATA_ROOT, OUTPUT_ROOT
import csv
from pathlib import Path

out = OUTPUT_ROOT / '09_figures_tables/outputs'
out.mkdir(parents=True, exist_ok=True)
with open(DATA_ROOT / 'manifests/figure_table_manifest.csv', newline='') as f:
    rows = list(csv.DictReader(f))
for row in rows:
    row['final_exists'] = (OUTPUT_ROOT / row['final_file']).exists()
    row['ready'] = row['status'] == 'completed' and row['final_exists']
with open(out/'figure_table_audit.csv','w',newline='') as f:
    writer=csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader(); writer.writerows(rows)
(out/'audit_status.txt').write_text('status=PLANNING_ONLY\nfinal_outputs=0\nlegacy_files_modified=false\n')
print('PLANNING_ONLY')
