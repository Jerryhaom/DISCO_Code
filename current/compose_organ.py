"""Vector composition using the verified original a-c/d panel geometry."""
from pathlib import Path
from io import BytesIO
import argparse,json
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
parser=argparse.ArgumentParser()
parser.add_argument('--abc',type=Path,required=True)
parser.add_argument('--d',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
a=parser.parse_args()
if a.output.exists(): raise SystemExit('Refusing to replace an existing figure.')
paths={n:a.abc/'panels'/name for n,name in dict(a='a_correlations.pdf',b='b_DISCO_HR.pdf',c='c_ProtAge_HR.pdf').items()}
paths['d']=a.d/'panels/d_Cindex.pdf'
W=691.2;ah=158.4;bh=295.2;dh=331.2;gap=9;head=22;H=head+ah+bh+gap+dh
w=PdfWriter();page=w.add_blank_page(W,H)
for n,(x,y) in {'d':(0,0),'b':(0,dh+gap),'c':(W/2,dh+gap),'a':(0,dh+gap+bh)}.items():
 p=PdfReader(paths[n]);assert len(p.pages)==1
 page.merge_transformed_page(p.pages[0],Transformation().translate(x,y))
b=BytesIO();c=canvas.Canvas(b,pagesize=(W,H));c.setFont('Helvetica',11)
c.drawString(13.95,H-14,'Organ protein age and DISCO');c.setFont('Helvetica-Bold',12)
for tag,x,y in [('a',10,H-head-12),('b',10,dh+gap+bh-12),('c',W/2+10,dh+gap+bh-12),('d',10,dh-16)]:c.drawString(x,y,tag)
c.showPage();c.save();page.merge_page(PdfReader(BytesIO(b.getvalue())).pages[0])
w.add_metadata({'/Title':'Organ protein age and DISCO','/Subject':'Locally rebuilt clocks; common test n=4294, deaths=460. d Original: DISCO + age; ProtAge alone; no sex. Accel: training-calibrated age residuals.'})
a.output.parent.mkdir(parents=True,exist_ok=True);w.write(a.output)
t=PdfReader(a.output).pages[0].extract_text()
assert t.count('Original')==9 and t.count('Accel')==9 and t.count('DISCO + Age')==9
print('Composed four verified vector panels:',a.output)
