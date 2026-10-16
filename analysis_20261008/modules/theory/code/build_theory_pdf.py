import os
"""Align the preface and theory pages 1–8; retain actual simulation pages 9–16."""
from pathlib import Path
import hashlib,json,html
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from pypdf import PdfReader,PdfWriter

M=Path(__file__).resolve().parents[1]
SRC=M/'input/SELECTED_MathematicalSupplement.pdf'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
EXPECTED='ed55e3d733e9785f0d521c1dabc559fa71a0d834b272951deb5fdbd93f5ce00a'
assert sha(SRC)==EXPECTED
FONT=os.environ['DISCO_THEORY_FONT']
pdfmetrics.registerFont(TTFont('Theory',FONT))
STYLE=ParagraphStyle('body',fontName='Theory',fontSize=10.5,leading=14.6,spaceAfter=8)
HEAD=ParagraphStyle('heading',fontName='Theory',fontSize=13,leading=17,spaceAfter=12)
SUB=ParagraphStyle('subheading',fontName='Theory',fontSize=11.5,leading=16,spaceAfter=9)
FORM=ParagraphStyle('formula',fontName='Theory',fontSize=10.25,leading=15.0,leftIndent=8,spaceAfter=10)
WIDTH,HEIGHT=595.2,841.92
pages=[]
def page(title,blocks):pages.append((title,blocks))
def t(s):return ('body',html.escape(s))
def f(s):return ('formula',s)
def h(s):return ('subheading',html.escape(s))

page('Supplementary mathematics: purpose and interpretation',[
 h('Why compare distributional summaries?'),
 t('Biomarkers are continuous, dependent measurements. Shannon entropy applies to a discrete distribution; differential entropy describes a continuous distribution. Both concern a specified distribution, not an intrinsic entropy value assigned to a single cross-sectional observation. A complete joint biological distribution is generally not observed in cohort biomarker data.'),
 t('Age-related changes can involve means, variances and dependence. A reference-based score can summarize how an individual differs from a declared reference or perturbs its fitted second-moment structure. This does not, by itself, measure thermodynamic entropy or reconstruct an individual’s biological entropy.'),
 h('Working hypothesis, not an identity'),
 t('We hypothesize that multivariate deviations from an explicitly defined reference may capture dysregulation associated with aging. Whether those deviations reflect biological entropy requires additional assumptions and evidence. The cohort-specific reference groups need not be young or disease-free.'),
 t('KL and Rényi divergences compare probability distributions. They are not entropy changes, and a positive divergence does not imply that the new distribution has higher entropy. Covariance log-determinants give the differential entropy of a fitted nonsingular Gaussian proxy; they are not the exact entropy of an unspecified biological distribution.'),
 h('Structure of this supplement'),
 t('Page 2 defines the four summaries. Pages 3–8 give their exact finite-reference relationships, conditional approximations and limitations. Pages 9–16 retain the corrected simulation figures. These specified Gaussian illustrations examine geometry and apparent prediction; they do not establish universal equivalence, biological entropy or causality.')])

page('Definitions and comparison roadmap',[
 t('We compare four second-moment summaries of individual deviation and reference perturbation: Mahalanobis distance, DISCO, fitted-Gaussian entropy change and directed Gaussian KL divergence. Their relationships with entropy and divergence require the definitions and conditions below.'),
 h('Mahalanobis distance'),
 t('Mahalanobis distance measures covariance-scaled deviation from the reference mean and requires an invertible covariance. It can be calculated without Gaussian data; Gaussian likelihood and entropy interpretations require Gaussian assumptions. It is not universally robust to outliers or singular covariance.'),
 f('d<sub>M</sub><super>2</super> = (x − μ)<super>T</super>S<super>−1</super>(x − μ).'),
 h('DISCO'),
 t('DISCO measures Pearson-correlation perturbation when x is appended to a reference of n observations. It can be calculated without a Gaussian model, but requires positive reference marginal variances and defined weights. Main age-weighted and historical unweighted simulation conventions are kept distinct.'),
 f('DISCO(x) = log[n<super>2</super> Σ<sub>i,j</sub>w<sub>ij</sub>(R<sub>ref,ij</sub> − R<sub>new,ij</sub>)<super>2</super>].'),
 t('The main weights are nonnegative, diagonal-zero and normalized over the full matrix. Conditional comparisons with Gaussian entropy and divergence do not make this logged, age-weighted score a general entropy estimator.'),
 h('Signed fitted-Gaussian entropy change'),
 f('ΔH = [log det S<sub>new</sub> − log det S<sub>ref</sub>]/2.'),
 t('The displayed corrected S1–S7 simulations use signed ΔH, not its absolute value. The historical ED1 producer uses twice this quantity, the full log-determinant difference. The nonnegative absolute difference has a different ranking and is not substituted silently.'),
 h('Directed Gaussian KL divergence'),
 t('The simulations compare the augmented fitted Gaussian to the reference fitted Gaussian, new || reference, including the change in mean. Exact finite-n identities and ridge conventions are stated on pages 3–8. Local Rényi comparisons use their separately declared direction and assumptions.'),
 t('Observed rank agreement in the specified simulations is not general functional equivalence or an advantage for every distribution. Shared use of second moments alone does not establish equality among the four metrics.')])

page('1. Signed Gaussian entropy and finite-reference updates',[
 t('DISCO quantifies Pearson-correlation perturbation. The entropy hypotheses in the introduction are not identities between this score and biological entropy. Mahalanobis distance itself can be computed without a Gaussian distribution; its Gaussian likelihood and entropy interpretations require the corresponding distributional assumptions.'),
 t('Let a fixed reference contain n ≥ 2 observations in p dimensions, with mean μ and positive-definite covariance S calculated with denominator n − 1. Append x. Put v = x − μ and q = vᵀS⁻¹v = dM². The augmented mean is μ₊ = μ + v/(n + 1). Centering about this new mean gives the scatter update below.'),
 f('(n)S<sub>+</sub> = (n − 1)S + n vv<super>T</super>/(n + 1),'),
 f('S<sub>+</sub> = aS + bvv<super>T</super>,   a = (n − 1)/n,   b = 1/(n + 1),<br/>g = b/a = n/[(n − 1)(n + 1)].'),
 t('The matrix determinant lemma gives an exact finite-n result. It does not require common ordinary eigenvectors of S and S₊.'),
 f('det(S<sub>+</sub>)/det(S) = a<super>p</super>(1 + gq).'),
 t('For a fitted nonsingular Gaussian, H = [p log(2πe) + log det S]/2. Accordingly the signed entropy benchmark actually used in Figs S1–S7 is'),
 f('ΔH = [log det S<sub>+</sub> − log det S]/2<br/>      = [p log a + log(1 + gq)]/2.'),
 t('For fixed n and p this is strictly increasing in q, because its derivative is g/[2(1 + gq)] > 0. It is an invertible monotone transform of squared Mahalanobis distance, not the same numerical quantity.'),
 t('The nonnegative difference dE = |ΔH| is a different quantity. It changes monotonic direction at q₀ = (a⁻ᵖ − 1)/g. Thus an absolute entropy difference need not preserve the ranking by dM. It is not used in the displayed signed-entropy benchmark.'),
 t('If covariance is instead computed with denominator n, the correct update has a = n/(n + 1), b = n/(n + 1)² and g = 1/(n + 1). This convention must not be mixed with the unbiased covariance used by R cov().'),
 t('Outside Gaussian distributions, covariance log-determinants describe the entropy of a fitted Gaussian proxy, not the exact entropy of the unknown data distribution.')])

page('2. Directed Gaussian KL and the mean term',[
 t('The simulation computes directed augmented-to-reference KL, DKL(N₊ || Nref), not a symmetrized divergence. With δμ = μ₊ − μ, the general formula is'),
 f('D<sub>KL</sub> = [tr(S<super>−1</super>S<sub>+</sub>) − p + δμ<super>T</super>S<super>−1</super>δμ<br/>                 + log det S − log det S<sub>+</sub>]/2.'),
 t('Let B = S⁻¹/² S₊ S⁻¹/², and let λᵢ be the eigenvalues of B, equivalently the generalized eigenvalues in S₊u = λSu. Since tr(S⁻¹S₊) = tr B and det B = det S₊/det S,'),
 f('D<sub>KL</sub> = (1/2)Σ<sub>i</sub>[λ<sub>i</sub> − 1 − log λ<sub>i</sub>]<br/>                 + (1/2)δμ<super>T</super>S<super>−1</super>δμ.'),
 t('This is valid for arbitrary positive-definite covariances. Dividing the ordinary eigenvalues of S₊ by those of S is valid only with a common eigenbasis and correct matching. Augmentation need not preserve that basis.'),
 t('For the unregularized rank-one unbiased update in Section 1, B has eigenvalue a repeated p − 1 times and eigenvalue a + bq once. Also δμᵀS⁻¹δμ = q/(n + 1)². Substitution gives'),
 f('D<sub>KL</sub> = [p(a − 1 − log a) + bq − log(1 + gq)<br/>                 + q/(n + 1)<super>2</super>]/2.'),
 t('The last term is the change in mean. In the bounded-q, large-n setting it can have the same inverse-square order as the covariance-only KL term; it cannot generally be dropped as negligible.'),
 t('If two Gaussian distributions have exactly the same covariance S, their KL divergence is half the squared Mahalanobis distance between their means. It is not necessarily half the individual observation’s dM². In the augmented sample, the mean displacement is v/(n + 1).'),
 f('Same covariance: D<sub>KL</sub>(P || Q) = (μ<sub>P</sub> − μ<sub>Q</sub>)<super>T</super><br/>S<super>−1</super>(μ<sub>P</sub> − μ<sub>Q</sub>)/2.'),
 t('Even the complete Gaussian KL need not be globally increasing in q under the unbiased update. Its derivative is negative for 0 ≤ q < 2(n + 1)/[n(n + 2)] and positive above that threshold. The rank correlations shown in the specified simulations therefore do not establish general functional equivalence.'),
 t('These identities use unregularized fitted covariances. Historical KL adds a diagonal ridge, as stated in Section 4.1, which need not preserve the exact rank-one form. KL measures relative distributional difference, not entropy change alone; positive KL does not imply larger entropy.')])

page('3. Local KL and Rényi relationships',[
 h('3.1. A controlled KL approximation'),
 t('Put ΔS = S₊ − S and A = S⁻¹/²ΔS S⁻¹/². For symmetric A with ||A||₂ ≤ η < 1, retaining the mean term gives'),
 f('D<sub>KL</sub>(new || ref) = δμ<super>T</super>S<super>−1</super>δμ/2 + ||A||<sub>F</sub><super>2</super>/4 + R,<br/>|R| ≤ η||A||<sub>F</sub><super>2</super>/[6(1 − η)].'),
 t('Expand e − log(1 + e) for each eigenvalue e of A. The terms from order three have absolute sum at most |e|³/[3(1 − η)], giving the stated bound after multiplication by 1/2. The squared norm equals tr[(S⁻¹ΔS)²] exactly. This is a reference-whitened relative perturbation, not an arbitrary raw covariance or Pearson difference. The mean term vanishes only if δμ = 0.'),
 h('3.2. Gaussian Rényi divergence and its domain'),
 t('For P = N(μP, ΣP), Q = N(μQ, ΣQ), α > 0 and α ≠ 1, define Cα = (1 − α)ΣP + αΣQ. For positive-definite Cα the Gaussian integral gives'),
 f('D<sub>α</sub>(P || Q) = [α log det Σ<sub>Q</sub> + (1 − α)log det Σ<sub>P</sub><br/>                          − log det C<sub>α</sub>]/[2(α − 1)]<br/>                    + (α/2)(μ<sub>P</sub> − μ<sub>Q</sub>)<super>T</super>C<sub>α</sub><super>−1</super>(μ<sub>P</sub> − μ<sub>Q</sub>).'),
 t('To derive this expression, complete the square in ∫p(x)^αq(x)^(1−α)dx with precision M = αΣP⁻¹ + (1−α)ΣQ⁻¹ and linear term αΣP⁻¹μP + (1−α)ΣQ⁻¹μQ. M and Cα have equivalent positive-definiteness conditions. For 0 < α < 1 the domain is automatic; for α > 1 it is not. Failure of this condition means a divergent integral, not a finite score.'),
 h('3.3. A restricted correlation-quadratic correspondence'),
 t('For zero-mean correlation models P = N(0, I), Q = N(0, I + E), with symmetric E, zero diagonal and ||E||₂ < 1/2,'),
 f('D<sub>2</sub>(ref || new) = [2log det(I + E) − log det(I + 2E)]/2<br/>                            = ||E||<sub>F</sub><super>2</super>/2 + O(Σ<sub>i</sub>|e<sub>i</sub>|<super>3</super>).'),
 t('This follows from the two log-determinant series. It connects an unlogged equal-weight correlation quadratic to a specific local Gaussian divergence. It does not equate age-weighted log DISCO to general KL or Rényi, and does not cover augmented means or marginal variances. The direction ref || new here also differs from the simulation KL direction.')])

page('4. Pearson perturbation and reference-size scaling',[
 t('Suppose all reference marginal standard deviations sᵢ are positive. Set zᵢ = (xᵢ − μᵢ)/sᵢ and let R be the reference Pearson matrix. From the exact covariance update,'),
 f('R<sub>+,ij</sub> = (R<sub>ij</sub> + g z<sub>i</sub>z<sub>j</sub>)/<br/>                  sqrt[(1 + g z<sub>i</sub><super>2</super>)(1 + g z<sub>j</sub><super>2</super>)].'),
 t('Thus Pearson perturbation changes with the marginal scales as well as covariance. It is not obtained by simply subtracting two covariances or by a PCA rotation.'),
 f('T<sub>ij</sub> = z<sub>i</sub>z<sub>j</sub> − R<sub>ij</sub>(z<sub>i</sub><super>2</super> + z<sub>j</sub><super>2</super>)/2,<br/>ΔR<sub>ij</sub> = gT<sub>ij</sub> + O(g<super>2</super>).'),
 t('For fixed p, bounded standardized x, and fixed nonnegative reference weights with wᵢᵢ = 0 and Σᵢⱼwᵢⱼ = 1, define'),
 f('D<sub>raw</sub> = Σ<sub>i,j</sub>w<sub>ij</sub>(ΔR<sub>ij</sub>)<super>2</super>.<br/>D<sub>raw</sub> = g<super>2</super>Σ<sub>i,j</sub>w<sub>ij</sub>T<sub>ij</sub><super>2</super> + O(g<super>3</super>).'),
 t('Because ng → 1, n²Draw approaches the weighted sum of Tᵢⱼ². When the reference changes randomly, a limiting claim also requires convergent means, finite positive marginal variances, correlations and weights, and control of the added observation’s standardized deviation. It is not exact finite-sample invariance or a uniform guarantee for observations moving into the tails as n grows.'),
 t('An explicit elementwise bound is available. If max|zᵢ| ≤ M and 0 ≤ g ≤ g₀, then'),
 f('|ΔR<sub>ij</sub> − gT<sub>ij</sub>| ≤ Cg<super>2</super>,<br/>C = (5/2)M<super>4</super> + (3/2)g<sub>0</sub>M<super>6</super>.'),
 t('For f(t) = (r + tc)/sqrt[(1 + ta)(1 + tb)], use |r| ≤ 1, |c| ≤ M² and a,b ≤ M². Direct differentiation gives |f″(t)| ≤ 5M⁴ + 3g₀M⁶. Taylor’s theorem gives the bound; |T| ≤ 2M² then gives'),
 f('|D<sub>raw</sub> − g<super>2</super>Σw<sub>ij</sub>T<sub>ij</sub><super>2</super>| ≤ 4g<super>3</super>M<super>2</super>C + g<super>4</super>C<super>2</super>.'),
 t('The main score is log(n²Draw), not log(1 + n²Draw). A finite logarithmic limit requires a positive limiting quadratic. If Draw = 0, the original score is −∞; a positive floor would define a different algorithm.')])

page('4.1. Weights, non-equivalence and defined domains',[
 t('For the main age-weighted algorithm, let aᵢ be the biomarker–age Pearson correlation in the declared weight-estimation population. The weights are'),
 f('v<sub>ij</sub> = |a<sub>i</sub>a<sub>j</sub>| for i ≠ j;  v<sub>ii</sub> = 0,<br/>w<sub>ij</sub> = v<sub>ij</sub>/Σ<sub>u,v</sub>v<sub>uv</sub>.'),
 t('The denominator must be positive. Normalization and summation are over the full off-diagonal matrix. Using only one triangle with the same weights halves Draw and shifts the log score by −log 2. In held-out analyses, weights are learned in training and frozen; historical whole-cohort weights are not training-only.'),
 t('Normalized equal-pair weights are 1/[p(p − 1)]. The displayed simulations retain their original unweighted convention log[(n²/2)Σᵢⱼ(ΔRᵢⱼ)²]. Its raw quadratic differs from normalized equal-pair Draw by p(p − 1)/2. For the restricted Rényi example in Section 3.3 this same factor appears. Arbitrary age weights do not share that isotropic quadratic.'),
 h('A direct counterexample to general equivalence'),
 t('Append x = μ. Then q = 0, S₊ = aS and R₊ = R, so Draw = 0. Nevertheless,'),
 f('ΔH = p log a/2 &lt; 0,<br/>d<sub>E</sub> = |ΔH| &gt; 0,<br/>D<sub>KL</sub>(new || ref) = p(a − 1 − log a)/2 &gt; 0.'),
 t('For n = 20, p = 3, these equal −0.07693994, 0.07693994 and 0.00193994, respectively. Pearson perturbation is zero while the Gaussian quantities are not. This disproves a general identity or proportionality, even before considering non-Gaussian biology.'),
 h('Singularity and implementation conventions'),
 t('Sample covariance rank is at most min(p, n − 1). Nonsingular Gaussian entropy and inverse-covariance formulas require positive definiteness; they cannot be used on singular or indefinite covariance by taking an absolute determinant. Pearson DISCO need not have an invertible covariance, but requires positive marginal variances and a defined weight normalizer.'),
 t('Historical simulation implementations differ: DM adds 10⁻⁶ to the diagonal when det S < 10⁻¹⁰; KL adds 10⁻¹⁰ to both covariances; entropy uses the unregularized determinant. The displayed results retain those conventions. No silent floor, absolute determinant, unified ridge or log1p replacement has been introduced. Singular Gaussian formulas require a separately defined support or regularization analysis rather than an invented finite value.')])

page('5. Gaussian simulations and scope of interpretation',[
 h('5.1. Two-dimensional datasets'),
 t('The existing simulations generate 1,000 observations from a zero-mean bivariate normal distribution with covariance [[1, ρ], [ρ, 1]], for ρ = 0.1 and 0.7. A 100 × 100 grid spans 1.5 times the observed range. Each metric is separately normalized to [0, 1]; log(1 + KL) is applied before that normalization. This display transformation of KL is distinct from the original logarithmic DISCO definition.'),
 t('Entropy in Figs S1–S7 denotes signed ΔH, one half of the log-determinant difference. It does not denote the absolute difference dE. The following pages retain the actual corrected plots, values and original outcome generators.'),
 h('Independent checks of the mathematical identities'),
 t('Direct covariance recomputation in 120 matrix cases, covering n = 7, 30, 1,000, p = 2, 5 and both covariance denominators, agrees with the finite-n update and determinant identities. Maximum signed-entropy error was 8.72 × 10⁻¹⁵ and maximum Pearson-update error 3.33 × 10⁻¹⁶. Generalized-eigenvalue and full-matrix KL agreed within 6.40 × 10⁻¹⁴.'),
 t('The Gaussian Rényi covariance formula was independently checked against the completed-square integral formula in 40 cases at orders 0.5 and 2 (maximum difference 3.55 × 10⁻¹⁵). Decreasing perturbations in the identity-reference example recovered its leading equal-weight quadratic. These checks verify the stated formulas, not a universal entropy interpretation of DISCO.'),
 h('Checks against the actual displayed objects'),
 t('The eight S3/S4 settings each contain 1,000 observations at p = 10. All reference covariances were positive definite (minimum eigenvalue at least 0.14045; condition number at most 33.32). The exact signed formula reproduced stored entropy to 1.69 × 10⁻¹⁴; the exact Pearson update reproduced DISCO to 2.97 × 10⁻¹². A further 160 direct augmented-correlation calculations agreed to 2.56 × 10⁻¹⁵.'),
 t('Each setting had 555 negative signed entropy changes. Mahalanobis distance had Spearman correlation 1 with ΔH but −0.02241 with |ΔH|. Absolute entropy is therefore a diagnostic comparison here, not a replacement for the plotted benchmark.'),
 h('Geometric illustration versus independent validation'),
 t('S3/S4 score observations that also belong to the reference. S5 repeatedly perturbs the same 1,000 subjects over 11 noise levels. They describe geometry and responsiveness in the specified Gaussian design; repeated rows are not independent participants. Mortality simulations in S6/S7 use one dataset per cell and report apparent C-index ± 1 SE, not paired superiority tests. None of these designs establishes direct biological entropy, a causal aging mechanism, or superiority for every dataset and distribution.')])

(M/'results/theory_pages_spec.json').write_text(json.dumps([dict(page=i+1,title=title,blocks=blocks) for i,(title,blocks) in enumerate(pages)],ensure_ascii=False,indent=2)+'\n')
tmp=M/'results/theory_pages_1_8.pdf'
c=canvas.Canvas(str(tmp),pagesize=(WIDTH,HEIGHT));boxes=[]
for pageNo,(title,blocks) in enumerate(pages,1):
 y=HEIGHT-52
 for typ,content in [('heading',html.escape(title))]+blocks:
  style={'body':STYLE,'heading':HEAD,'subheading':SUB,'formula':FORM}[typ]
  para=Paragraph(content,style);w,ht=para.wrap(WIDTH-104,HEIGHT)
  if y-ht<49:raise RuntimeError(f'Page {pageNo} overflow at {content[:70]} y={y} h={ht}')
  para.drawOn(c,52,y-ht);boxes.append(dict(page=pageNo,type=typ,x=52,y=y-ht,width=w,height=ht,content=content))
  y-=ht+style.spaceAfter
 c.setFont('Theory',9);c.drawCentredString(WIDTH/2,30,str(pageNo));c.showPage()
c.save()
writer=PdfWriter();original=PdfReader(SRC);replacement=PdfReader(tmp)
assert len(original.pages)==16 and len(replacement.pages)==8
for i,p in enumerate(original.pages): writer.add_page(replacement.pages[i] if i<=7 else p)
out=M/'results/Mathematical_Supplement_EXACT_CONDITIONS_CANDIDATE.pdf'
with out.open('wb') as stream:writer.write(stream)
assert sha(SRC)==EXPECTED
prov=dict(source=str(SRC),sourceSHA256=EXPECTED,output=str(out),outputSHA256=sha(out),pages=16,
 onlyReplacedPages=list(range(1,9)),unchangedPages=list(range(9,17)),
 allSimulationPagesRetained=True,simulationValuesChanged=False,primaryAlgorithmChanged=False,
 bodyFont='DejaVu Sans 10.5 pt',equationsVectorText=True,layoutBounds=boxes,
 theoreticalScope='Exact signed entropy identity; non-equivalence counterexample; generalized Gaussian KL, mean term, controlled Taylor remainder, restricted Renyi, exact Pearson and finite-n domains.',
 sourcePrefaceAligned='Pages 1–2 qualify biological-entropy hypotheses and distinguish non-Gaussian computability from Gaussian interpretation.')
(M/'results/math_pdf_provenance.json').write_text(json.dumps(prov,ensure_ascii=False,indent=2)+'\n')
print(out);print(prov['outputSHA256']);print('PASS: eight preface/theory pages built without overflow; all eight simulation pages imported unchanged')
