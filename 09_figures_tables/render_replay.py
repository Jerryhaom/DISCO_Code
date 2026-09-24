from project_paths import CODE_ROOT, DATA_ROOT, OUTPUT_ROOT
"""Render each controlled R PDF and record page count, dimensions, and hash."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess, hashlib, json
from PIL import Image, ImageDraw

base = OUTPUT_ROOT / '09_figures_tables/outputs/original_block_replay'
out = base / 'qa'
out.mkdir(parents=True, exist_ok=True)
cached = {r['file']:r for r in json.loads((out/'inventory.json').read_text())} if (out/'inventory.json').exists() else {}
def render(path):
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    old=cached.get(str(path))
    if old and old['sha256']==digest and old['render_ok'] and (out/(path.stem+'.png')).exists():
        return old
    result = subprocess.run(['pdfinfo', str(path)],capture_output=True,text=True)
    proc = subprocess.run(['pdftoppm','-scale-to','1500','-singlefile','-png',str(path),str(out/path.stem)],capture_output=True,text=True)
    txt = subprocess.run(['pdftotext','-layout',str(path),'-'],capture_output=True,text=True)
    (out/(path.stem+'.txt')).write_text(txt.stdout)
    return {'file':str(path),'sha256':digest,
            'info':result.stdout,'render_ok':proc.returncode==0,'error':result.stderr+proc.stderr}
paths=sorted(base.glob('*.pdf'))
with ThreadPoolExecutor(max_workers=4) as pool:
    records=list(pool.map(render,paths))
(out/'inventory.json').write_text(json.dumps(records,indent=2))
for start in range(0,len(paths),6):
    subset=paths[start:start+6]
    canvas=Image.new('RGB',(1500,580*((len(subset)+1)//2)),'white')
    d=ImageDraw.Draw(canvas)
    for i,path in enumerate(subset):
        x,y=(i%2)*750,(i//2)*580
        d.text((x+8,y+6),path.stem,fill='black')
        png=out/(path.stem+'.png')
        if png.exists():
            im=Image.open(png).convert('RGB'); im.thumbnail((735,550))
            canvas.paste(im,(x+(750-im.width)//2,y+27))
    canvas.save(out/f'sheet_{start//6}.png')
print(json.dumps({'count':len(records),'failed':[r['file'] for r in records if not r['render_ok']]}))
