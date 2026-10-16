"""Insert the verified signed-half plots into the existing technical PDF candidate.

This is a new conditional candidate, not a theoretical approval or a submission
overwrite. Preserve the original PDF and keep new plots as vector PDF objects.
"""
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

ROOT = Path(_disco_path("results",""))
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / '09_figures_tables/checklist_actions_20261006/mathematics/tmp/pdf_edit_dependencies'))
import pymupdf as fitz

SOURCE = ROOT / '09_figures_tables/checklist_actions_20261006/mathematics/output/pdf/Mathematical_Supplement_TECHNICAL_CORRECTIONS_CANDIDATE.pdf'
SIM = HERE / 'pass_2'
OUT = HERE / 'Mathematical_Supplement_SIGNED_HALF_TECHNICAL_CANDIDATE.pdf'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
before = sha(SOURCE)
doc = fitz.open(SOURCE)
assert len(doc) == 16
edits = []
fonts = fitz.Archive('/System/Library/Fonts/Supplemental')
css = """
@font-face {font-family: TimesQA; src: url('Times New Roman.ttf');}
@font-face {font-family: TimesQA; font-weight: bold; src: url('Times New Roman Bold.ttf');}
body {font-family: TimesQA; font-size: 10.56pt; line-height: 1.5; margin: 0; color: #000;}
p {margin: 0; text-align: justify;} sub,sup {font-size: 72%;}
"""

def erase(page, box, graphics=False):
    rect = fitz.Rect(box)
    old = page.get_text('text', clip=rect)
    page.add_redact_annot(rect, fill=(1, 1, 1))
    page.apply_redactions(images=1 if graphics else 0, graphics=2 if graphics else 0, text=0)
    return rect, old

def caption(pno, box, text):
    page = doc[pno-1]
    rect, old = erase(page, box)
    spare, scale = page.insert_htmlbox(rect, '<p>'+text+'</p>', css=css, archive=fonts, scale_low=1)
    assert spare >= 0 and scale == 1, (pno, spare, scale)
    edits.append(dict(page=pno,rect=list(rect),old_text=old,new_html=text,scale=scale))

def plot(pno, box, name):
    page = doc[pno-1]
    rect, old = erase(page, box, graphics=True)
    with fitz.open(SIM / name) as plotdoc:
        assert len(plotdoc) == 1
        page.show_pdf_page(rect, plotdoc, 0, keep_proportion=True)
    edits.append(dict(page=pno,rect=list(rect),plot=name,sha256=sha(SIM/name),vector=True))

caption(3,(89,513,513,772),
  'Let &Sigma;<sub>ref</sub> and &Sigma;<sub>new</sub> be the reference and augmented covariance matrices. '
  'The nonnegative theoretical Gaussian entropy distance is '
  'd<sub>E</sub>(x) = (1/2)|logdet(&Sigma;<sub>new</sub>) - logdet(&Sigma;<sub>ref</sub>)|. '
  'It differs from the signed Gaussian entropy change &Delta;H = '
  '(1/2)[logdet(&Sigma;<sub>new</sub>) - logdet(&Sigma;<sub>ref</sub>)]. '
  'The technically corrected S1-S5 branch and the original S6-S7 calculations use &Delta;H, '
  'without an absolute value. Multiplication by 1/2 preserves ranks and separate min-max '
  'normalization but does not resolve the difference between signed and absolute quantities. '
  'Outside Gaussian distributions, log-determinant change is not an exact general entropy measure. '
  'The final entropy definition and theoretical equivalence require author review.<br/><br/>'
  '<b>(4). KL distance</b><br/>'
  'The simulation computes directed KL divergence from the distribution estimated '
  'after adding an individual to the original reference distribution. For Gaussian '
  'distributions P ~ N(&mu;<sub>P</sub>, &Sigma;<sub>P</sub>) and '
  'Q ~ N(&mu;<sub>Q</sub>, &Sigma;<sub>Q</sub>), KL is defined as follows:')
caption(8,(89,613,514,731),
  '<b>5.1. 2-dimensional dataset:</b> 1,000 observations were generated from a bivariate '
  'normal distribution with zero means and covariance [[1, &rho;], [&rho;, 1]], for '
  '&rho; = 0.1 and 0.7. Each metric was evaluated on a 100 &times; 100 grid spanning '
  '1.5 times the observed range. Scores were separately normalized to [0, 1], '
  'with log(1 + KL) applied before normalization. Entropy denotes the signed Gaussian '
  'entropy change &Delta;H, not the nonnegative theoretical d<sub>E</sub>.')
caption(9,(89,336,528,405),
  '<b>Figure S1.</b> Contour distributions of DISCO, Mahalanobis, Entropy and directed '
  'KL divergence in two dimensions for &rho; = 0.1 (a) and 0.7 (b). Each metric is '
  'separately normalized to [0, 1]. Entropy denotes signed Gaussian entropy change; '
  'its factor 1/2 cancels under this normalization.')
plot(10,(89,305,514,625),'Cmp.DISCO_with_KL_high.pdf')
caption(10,(89,639,514,733),
  '<b>Figure S2.</b> Contours of DISCO, Mahalanobis, Entropy and KL in the PC1-PC2 '
  'plane of the 10-dimensional datasets (&rho; = 0.1 and 0.7). Other PCs are zero. '
  'The corrected inverse PCA restores feature scales and means. Colors are separately '
  'normalized metric values, not common raw units. Entropy is the signed Gaussian '
  'entropy change &Delta;H.')
plot(11,(89,74,514,492),'DM3KL.pdf')
caption(11,(89,501,514,609),
  '<b>Figure S3.</b> Pairwise Spearman correlations for 1,000 generated 10-dimensional '
  'observations at each Toeplitz parameter: 0.1, 0.3, 0.5 and 0.7, top to bottom. '
  'Entropy is &Delta;H, half the signed log-determinant change. Removing the historical '
  'tenfold row recycling and multiplying entropy by 1/2 preserve the displayed rank '
  'correlations. These observations also belong to the reference set.')
plot(12,(89,74,514,492),'DM2KL.pdf')
caption(12,(89,501,514,581),
  '<b>Figure S4.</b> DISCO, Mahalanobis and signed Gaussian entropy change (&Delta;H) '
  'versus directed KL divergence for 1,000 observations per setting. Rows use '
  '&rho; = 0.1, 0.3, 0.5 and 0.7. Spearman correlations are displayed. '
  'KL includes the augmented mean and covariance terms.')
plot(13,(89,168,514,475),'Noise.pdf')
caption(13,(89,482,514,544),
  '<b>Figure S5.</b> DISCO, Mahalanobis, signed Gaussian entropy change and KL under '
  'Gaussian noise. Noise levels are standard deviations (0-2). Correlations describe '
  'repeated perturbations of 1,000 subjects; independent-row P values are not shown.')
doc.save(OUT,garbage=4,deflate=True)
doc.close()
assert sha(SOURCE) == before
# Confirm every page outside the explicitly changed set is pixel-identical.
changed = {3,8,9,10,11,12,13}
qa = []
with fitz.open(SOURCE) as old, fitz.open(OUT) as new:
    for i in range(16):
        if i+1 not in changed:
            a = old[i].get_pixmap(matrix=fitz.Matrix(1,1),alpha=False)
            b = new[i].get_pixmap(matrix=fitz.Matrix(1,1),alpha=False)
            assert a.samples == b.samples, ('unintended page change',i+1)
            qa.append(dict(page=i+1,unchanged_pixels=True))
(HERE/'math_pdf_replacement_log.json').write_text(json.dumps(dict(source=str(SOURCE),source_sha256=before,
  output=str(OUT),output_sha256=sha(OUT),edits=edits,unchanged_page_checks=qa,
  status='TECHNICAL_CONDITIONAL_CANDIDATE_NOT_THEORETICALLY_APPROVED'),ensure_ascii=False,indent=2))
print(OUT)
