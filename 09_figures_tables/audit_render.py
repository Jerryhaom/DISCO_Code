from project_paths import CODE_ROOT, DATA_ROOT, OUTPUT_ROOT
"""Read-only PDF evidence capture. Original PDFs are never modified."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import subprocess
from PIL import Image, ImageOps, ImageDraw

ROOT = OUTPUT_ROOT
OUT = ROOT / '09_figures_tables/outputs/content_audit'
SOURCE = DATA_ROOT / 'Manuscript/PNAS/Submission_Rev1/Figures.pdf'
OUT.mkdir(parents=True, exist_ok=True)

def run(args):
    return subprocess.run(args, capture_output=True, text=True)

def render(path, key, page=1, width=950):
    prefix = OUT / key
    info = run(['pdfinfo', str(path)])
    rec = {'file': str(path), 'page': page, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
           'pdfinfo': info.stdout, 'error': info.stderr, 'rendered': False}
    tx = run(['pdftotext', '-f', str(page), '-l', str(page), '-layout', str(path), '-'])
    prefix.with_suffix('.txt').write_text(tx.stdout)
    if info.returncode == 0:
        proc = run(['pdftoppm', '-f', str(page), '-l', str(page), '-scale-to', str(width),
                    '-png', '-singlefile', str(path), str(prefix)])
        rec['rendered'] = proc.returncode == 0
        rec['error'] += proc.stderr
    return key, rec

def sheet(keys, name, cols=2, cell=(730, 540)):
    canvas = Image.new('RGB', (cols * cell[0], ((len(keys)+cols-1)//cols)*cell[1]), 'white')
    draw = ImageDraw.Draw(canvas)
    for n, key in enumerate(keys):
        x, y = (n % cols)*cell[0], (n//cols)*cell[1]
        draw.text((x+8, y+5), key, fill='black')
        path = OUT / (key+'.png')
        if path.exists():
            im = Image.open(path).convert('RGB')
            im.thumbnail((cell[0]-16,cell[1]-30))
            canvas.paste(im, (x+(cell[0]-im.width)//2, y+26))
        else:
            draw.text((x+20,y+100), 'INVALID / NOT RENDERED', fill='red')
    canvas.save(OUT/(name+'.png'))

jobs = [(SOURCE, f'submission_{p:02}', p, 1450) for p in range(1,16)]
groups = {}
for module in sorted(ROOT.glob('0[1-8]_*')):
    pdfs = sorted(p for p in module.rglob('*.pdf') if p.name != 'Rplots.pdf')
    keys=[]
    for i, path in enumerate(pdfs):
        key=f'{module.name[:2]}_{i:02}_{path.stem}'
        keys.append(key)
        jobs.append((path,key,1,1050))
    groups[module.name]=keys

records={}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for key, rec in pool.map(lambda args: render(*args),jobs):
        records[key]=rec
for module,keys in groups.items():
    for start in range(0,len(keys),6):
        sheet(keys[start:start+6],f'sheet_{module[:2]}_{start//6}',cols=2)
(OUT/'inventory.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
(OUT/'groups.json').write_text(json.dumps(groups,indent=2))
print(json.dumps({'pdfs':len(jobs),'failed':[k for k,r in records.items() if not r['rendered']], 'out':str(OUT)}))
