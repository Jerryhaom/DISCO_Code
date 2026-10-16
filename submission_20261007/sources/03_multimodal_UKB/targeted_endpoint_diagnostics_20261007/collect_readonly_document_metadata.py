#!/usr/bin/env python3
"""Read-only snapshot of candidates and already rendered PDFs; no document edits."""
# File-system roots are explicit; no participant data are bundled.
import os as _disco_os
def _disco_path(kind, suffix=""):
    key = {"private": "DISCO_PRIVATE_ROOT", "results": "DISCO_RESULTS_ROOT", "legacy": "DISCO_LEGACY_ROOT"}[kind]
    root = _disco_os.environ.get(key)
    if not root:
        raise RuntimeError("Required environment variable: " + key)
    return str(__import__("pathlib").Path(root) / suffix) if suffix else root
import hashlib, json, re, sys, zipfile
from datetime import datetime
from pathlib import Path
from lxml import etree
from pypdf import PdfReader

root = Path(_disco_path("results","09_figures_tables/materials_sync_20261007"))
ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
result = {'created': datetime.now().isoformat(timespec='seconds'), 'read_only': True,
          'important_limit': 'Existing v2 renders were individually viewed, but newer DOCX snapshots are not visually certified by those earlier renders.', 'files': []}
for label in ('MAIN', 'SI', 'RESPONSE'):
    for relative in (f'documents/{label}_synchronized_CANDIDATE_TRACKED.docx',
                     f'qa_accepted/{label}_synchronized_CANDIDATE_QA_ONLY.docx',
                     f'qa/render_{"main" if label == "MAIN" else "si" if label == "SI" else "response"}_v2/{label}_synchronized_CANDIDATE_QA_ONLY.pdf'):
        path = root / relative
        item = {'path': str(path), 'bytes': path.stat().st_size,
                'mtime_local': datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec='seconds'),
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        if path.suffix == '.docx':
            with zipfile.ZipFile(path) as z:
                node = etree.fromstring(z.read('word/document.xml'))
                paragraphs = []
                for p in node.xpath('.//w:p', namespaces=ns):
                    texts = p.xpath('.//w:t[not(ancestor::w:del)]/text()', namespaces=ns)
                    if texts: paragraphs.append(''.join(texts))
                effective = '\n'.join(paragraphs)
                item['effective_keyword_counts'] = {term: len(re.findall(re.escape(term), effective))
                    for term in ('4,294', '18,648', '17 sets', '615', '7,030', '7030', '31,381', '344', 'CIF', 'References')}
                item['short_keyword_contexts'] = [effective[max(0, m.start()-65):m.end()+100]
                    for m in re.finditer(r'4,294|18,648|17 sets|615|7,030|7030|31,381|344', effective)]
                item['native_field_count'] = len(node.xpath('.//w:fldChar|.//w:fldSimple', namespaces=ns))
        else:
            pdf = PdfReader(path)
            item['pages'] = len(pdf.pages)
            item['individually_viewed_pages'] = list(range(1, len(pdf.pages) + 1))
            item['view_mode'] = 'page PNG, original image detail; not inferred from XML'
        result['files'].append(item)
json.dump(result, sys.stdout, ensure_ascii=False, indent=2)
print()
