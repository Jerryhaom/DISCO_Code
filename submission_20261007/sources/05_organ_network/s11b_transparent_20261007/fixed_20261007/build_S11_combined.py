# File-system roots are explicit; no participant data are bundled.
import os as _disco_os
def _disco_path(kind, suffix=""):
    key = {"private": "DISCO_PRIVATE_ROOT", "results": "DISCO_RESULTS_ROOT", "legacy": "DISCO_LEGACY_ROOT"}[kind]
    root = _disco_os.environ.get(key)
    if not root:
        raise RuntimeError("Required environment variable: " + key)
    return str(__import__("pathlib").Path(root) / suffix) if suffix else root
from pathlib import Path
import hashlib
import json
import sys
import zipfile

sys.path.insert(0,_disco_path("results","09_figures_tables/checklist_actions_20261006/mathematics/tmp/pdf_edit_dependencies"))
import fitz

base=Path(__file__).resolve().parent
source=Path(_disco_path("results","09_figures_tables/materials_sync_20261007/documents/SI_synchronized_CANDIDATE_TRACKED.docx"))
panel=base/'S11b_weighted_power_CORRECTED.pdf'
sha=lambda b:hashlib.sha256(b).hexdigest()
source_hash=sha(source.read_bytes())
with zipfile.ZipFile(source) as z:
    raw=z.read('word/media/supplementary_figure_S11.png')
# The source drawing is an entire raster page, including the old caption and
# watermark. Retain the plotted a/c/d/e pixels, not its obsolete b caption.
native=fitz.open();p=native.new_page(width=1653,height=2520)
p.insert_image(p.rect,stream=raw)
crop=fitz.Rect(320,345,1304,1700)
clear=fitz.Rect(365,986,819,1305)
replacement=fitz.Rect(365,986,815,1286)
before=p.get_pixmap(matrix=fitz.Matrix(1,1),clip=crop,alpha=False)
p.draw_rect(clear,color=(1,1,1),fill=(1,1,1),overlay=True)
with fitz.open(panel) as b:
    p.show_pdf_page(replacement,b,0,keep_proportion=True,overlay=True)
after=p.get_pixmap(matrix=fitz.Matrix(1,1),clip=crop,alpha=False)
assert before.width==after.width and before.height==after.height
# Compare rendered raster pixels outside the replaced b rectangle exactly.
left=int(clear.x0-crop.x0)-2;right=int(clear.x1-crop.x0)+2
top=int(clear.y0-crop.y0)-2;bottom=int(clear.y1-crop.y0)+2
old=before.samples;new=after.samples;stride=before.stride
outside_diff=0
for y in range(before.height):
    a,b=old[y*stride:(y+1)*stride],new[y*stride:(y+1)*stride]
    if top<=y<=bottom:
        a=a[:left*3]+a[right*3:];b=b[:left*3]+b[right*3:]
    outside_diff+=sum(x!=y for x,y in zip(a,b))
assert outside_diff==0,outside_diff
combined=fitz.open();q=combined.new_page(width=504,height=504*crop.height/crop.width)
q.show_pdf_page(q.rect,native,0,clip=crop,keep_proportion=True)
output=base/'S11_COMBINED_corrected_b_CANDIDATE.pdf'
combined.save(output,garbage=4,deflate=True)
q.get_pixmap(matrix=fitz.Matrix(3,3),alpha=False).save(base/'S11_COMBINED_corrected_b_DOCX.png')
q.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).save(base/'S11_COMBINED_corrected_b_QA.png')
assert sha(source.read_bytes())==source_hash
record={'source_docx':str(source),'source_docx_sha256':source_hash,
 'source_media':'word/media/supplementary_figure_S11.png','source_media_sha256':sha(raw),
 'source_pixels':[1653,2520],'figure_crop_native_pixels':list(crop),
 'panel_b_clear_native_pixels':list(clear),'panel_b_replacement_native_pixels':list(replacement),
 'panel_b_source':str(panel),'panel_b_sha256':sha(panel.read_bytes()),
 'other_panel_rendered_channel_differences':outside_diff,
 'output_pdf':str(output),'output_pdf_sha256':sha(output.read_bytes()),
 'docx_png':str(base/'S11_COMBINED_corrected_b_DOCX.png'),
 'a_c_d_e_identity':'Outside b rectangle rendered pixels exactly preserved; no independent statistical approval supplied.',
 'composition':'b remains vector; a/c/d/e retain source raster pixels. Obsolete embedded caption/header/QA watermark are cropped out.',
 'caption_dependency':'Use S11_caption.txt in editable Word text; a/c/d/e depend on previously audited saved aggregates. Original b final upstream run remains unrecovered.',
 'source_preserved':True,'status':'TECHNICAL_REPLACEMENT_CANDIDATE_NOT_ALL_PANEL_SCIENTIFIC_APPROVAL'}
(base/'S11_combined_provenance.json').write_text(json.dumps(record,indent=2))
print('PASS: b vector replaced; all outside-b rendered pixels exactly preserved; DOCX source untouched')
