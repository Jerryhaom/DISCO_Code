# 10-biomarker clinical calculator example

This companion module uses **DISCO4HD >= 0.3.0**. It contains fictional profiles only, and does not access controlled cohort records.

```sh
Rscript run.R 10_biomarker_calculator/run_example.R
```

Outputs under `results/10_biomarker_calculator/outputs/`: available references, example DISCO / raw DM scores for all five references, a panel-matched CHARLS comparison using UKB parameters, and session information. `DISCO_OUTPUT_ROOT` can relocate outputs.

```r
library(DISCO4HD)
disco10_references()
ref <- disco10_reference("NHANES")
x <- data.frame(CRP=1, GLU=5, CHOL=4.5, HDL=1.3, CREA=80,
                UA=300, WBC=6, UREA=5, RBC=4.8, ALB=42)
cal_disco10(x, reference=ref)
cal_disco10(x, reference="UKB", comparison="NHANES")
```

`reference` chooses the actual covariance/correlation reference and fixed age weights. UKB, NHANES and CLHLS use ten markers; CHARLS uses eight (no RBC/ALB); RuLAS uses nine (no UREA). Optional `comparison` describes a UKB-reference score relative to a cohort's stored distribution and is supported only for the full, unchanged UKB reference. CHARLS/RuLAS comparisons recompute the corresponding reduced-panel UKB score. Do not interpret these as the own-reference calibration of those cohorts.

Required units: CRP mg/L; GLU/CHOL/HDL/UREA mmol/L; CREA/UA µmol/L; WBC 10⁹/L; RBC 10¹²/L; ALB g/L. CRP is logged internally as ln(CRP + 1). **DM is not logged**. Read the [full package guide](https://github.com/Jerryhaom/DISCO4HD/blob/main/docs/10-biomarker-guide.md) for equations, reference ages, sample sizes and limitations.

The small clinical panel is a practical calculator implementation. It does not replace the larger clinical/multiomics panels, historical scores or statistical models used throughout the article, and no formal biomarker optimality claim is made.
# New synthetic NHANES-calibrated demonstration

With the audited DISCO4HD >= 0.3.1 installed, run
`Rscript run.R 10_biomarker_calculator/run_synthetic_nhanes.R`.
This writes 400 fictional profiles, scores, generation provenance and a three-panel
PNG. No participant records are read. The full Markdown guide and runnable
visualization are in the companion package's `docs/synthetic-nhanes-example.md`.

Important: the web calculator's comparison-cohort selector changes the z-score
comparison distribution, **not** its fixed UKB scoring reference. In the R API,
`reference="NHANES"` really changes the reference. To reproduce the web's NHANES
comparison, use `reference="UKB", comparison="NHANES"`.
