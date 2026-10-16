#!/usr/bin/env python3
from project_paths import CODE_ROOT, DATA_ROOT, OUTPUT_ROOT
"""Organize figure reproduction outputs by analysis module."""
import os, shutil, json

BASE = str(OUTPUT_ROOT / "09_figures_tables")
OUT_BASE = str(OUTPUT_ROOT)
BASELINE = f"{BASE}/outputs/strict_runs/20260917T102326_baseline_all/pass_1"
CLIN = f"{BASE}/outputs/strict_extensions_clinical/20260917T104101_clinical_original_blocks_verified/pass_1"
EXTERN = f"{BASE}/outputs/strict_extensions_external/20260917T105013_external_original_blocks/pass_1"
ORGAN = f"{BASE}/outputs/strict_extensions_organ"
FRAIL = f"{BASE}/outputs/strict_extensions_frailty/20260917T105428_verified_original_context/pass_1"
RULAS = f"{BASE}/outputs/strict_extensions_rulas/20260917T024148_original_blocks/pass_1"
SIMGEO = f"{BASE}/outputs/strict_extensions_simulation_geometry/20260917T104833_original_S1_S5/pass_1"
SIMSURV = f"{BASE}/outputs/strict_extensions_simulation_survival/20260917T025312_original_four_metrics/pass_1"
RECONSIM = f"{BASE}/outputs/evidence_reconstruction_sim"
RECONROC = f"{BASE}/outputs/evidence_reconstruction_roc"
FIG3C = f"{BASE}/outputs/fig3c_reestimated"

# organ subruns
ORGAN_FIG4B = f"{ORGAN}/20260917T104203_verified_context/pass_1"
ORGAN_FIG4F = f"{ORGAN}/20260917T104111_corrected_source_blocks/pass_1"
ORGAN_ED9E = f"{ORGAN}/20260917T104203_verified_context/pass_1"
ORGAN_ED10A = f"{ORGAN}/20260917T104203_verified_context/pass_1"
ORGAN_ED9C = f"{ORGAN}/20260917T104111_corrected_source_blocks/pass_1"
ORGAN_ED9D = f"{ORGAN}/20260917T104111_corrected_source_blocks/pass_1"
ORGAN_ED10B = f"{ORGAN}/20260917T104430_sensitivity_source_blocks/pass_1"
ORGAN_ED10D = f"{ORGAN}/20260917T104430_sensitivity_source_blocks/pass_1"

# module -> figure name -> list of source paths (files or dirs)
MAP = {
 "01_method_definition": {
   "Fig1b": [f"{RECONSIM}/Fig1b_RECONSTRUCTED.pdf"],
   "ED1": [f"{BASELINE}/ED1__Ex.Fig1.pdf"],
   "ED2": [
      f"{BASELINE}/ED2b__UKB.Entropy.HD.pdf",
      f"{BASELINE}/ED2c__UKB.Entropy.HD.pdf",
      f"{BASELINE}/ED2d__Charls.Entropy.HD.pdf",
      f"{BASELINE}/ED2e__Rulas.Entropy.HD.pdf",
      f"{BASELINE}/ED2f__Clhls.Entropy.HD.pdf",
      f"{CLIN}/ED2a/NHANES.Entropy.HD.pdf",
   ],
   "ED4": [
      f"{BASELINE}/ED4a__NHANES.DISCO.Distribution.pdf",
      f"{BASELINE}/ED4b__DISCO.robust.pdf",
      f"{BASELINE}/ED4c__DISCO.refsel.pdf",
      f"{BASELINE}/ED4d__Bench.Distance.pdf",
      f"{BASELINE}/ED4e__Sim.clinical.pdf",
      f"{BASELINE}/ED4e__Sim.clinicalprot.pdf",
      f"{BASELINE}/ED4f__Leave_one_out2.pdf",
   ],
   "ED5": [f"{RECONSIM}/ED5_RECONSTRUCTED.pdf"],
   "S1_S5_geometry": [
      f"{SIMGEO}/Cmp.DISCO_with_KL_high.pdf",
      f"{SIMGEO}/Cmp.DISCO_with_KL_low.pdf",
      f"{SIMGEO}/DM2KL.pdf", f"{SIMGEO}/DM3KL.pdf", f"{SIMGEO}/Noise.pdf",
   ],
   "S6_S7_survival": [
      f"{SIMSURV}/Standard_cindex.pdf", f"{SIMSURV}/Highnoise_cindex.pdf",
   ],
 },
 "02_primary_clinical": {
   "Fig1c": [f"{CLIN}/Fig1c"],
   "Fig1d": [f"{CLIN}/Fig1d"],
   "Fig3a": [f"{BASELINE}/Fig3a__NHANES.cmp.v2.pdf"],
   "Fig3b": [f"{CLIN}/Fig3b_KM", f"{BASELINE}/Fig3b_age__Age.distribution.pdf"],
   "ED3a": [
      f"{BASELINE}/ED3a__UKB.DISCO.Age.pdf",
      f"{BASELINE}/ED3a__UKB.DISCO.HD.pdf",
      f"{BASELINE}/ED3a__UKB.DISCO2.all.pdf",
   ],
 },
 "03_multimodal_UKB": {
   "Fig2a": [f"{BASELINE}/Fig2a__UKB.Prot.DISCO.all.pdf"],
   "Fig2b": [f"{BASELINE}/Fig2b__UKB.ClinicalProt.DISCO.all.pdf"],
   "Fig3c": [f"{FIG3C}/Fig3c_REESTIMATED.pdf", f"{FIG3C}/Fig3c_48_estimates_vs_submission.csv", f"{FIG3C}/README.md"],
 },
 "04_external_validation": {
   "ED3b": [f"{EXTERN}/ED3b_CHARLS"],
   "ED3c": [f"{EXTERN}/ED3c_RuLAS"],
   "ED3d": [f"{EXTERN}/ED3d_CLHLS"],
   "ED6": [f"{BASELINE}/ED6a__HRS.cmp.v2.pdf", f"{CLIN}/ED6b"],
 },
 "05_organ_network": {
   "Fig4a": [f"{BASELINE}/Fig4a_legacy__Organ.outcomes2.pdf", f"{BASELINE}/Fig4a_legacy_plot_data.rds"],
   "Fig4b": [f"{ORGAN_FIG4B}/Fig4b"],
   "Fig4c": [f"{BASELINE}/Fig4c_CANDIDATE__Prot.de.pdf", f"{BASELINE}/Fig4c_CANDIDATE_plot_data.rds"],
   "Fig4d": [f"{BASELINE}/Fig4d__Organ.pcorr.pdf", f"{BASELINE}/Fig4d_original_pcorr_matrix.rds"],
   "Fig4e": [f"{BASELINE}/Fig4e_legacy__Partial.net.pdf", f"{BASELINE}/Fig4e_legacy_plot_data.rds"],
   "Fig4f": [f"{ORGAN_FIG4F}/Fig4f"],
   "ED9a": [f"{BASELINE}/ED9a_before__Prot.tissue.pdf", f"{BASELINE}/ED9a_selected__Prot.tissue.sel.pdf"],
   "ED9b": [f"{BASELINE}/ED9b__Tissue.Overlapped.prot.pdf"],
   "ED9c": [f"{ORGAN_ED9C}/ED9c"],
   "ED9d": [f"{ORGAN_ED9D}/ED9d"],
   "ED9e": [f"{ORGAN_ED9E}/ED9e"],
   "ED10a": [f"{ORGAN_ED10A}/ED10a"],
   "ED10b": [f"{ORGAN_ED10B}/ED10b"],
   "ED10c_ED11": [
      f"{BASELINE}/ED10c_ED11cde__Organ.cindex.degree.pdf",
      f"{BASELINE}/ED10c_ED11cde__Organ.cindex.sensitivity.pdf",
      f"{BASELINE}/ED10c_ED11cde__Organ.pcorr.cmp.pdf",
      f"{BASELINE}/ED10c_ED11cde__Organ.pcorr.sensitivity.pdf",
      f"{BASELINE}/ED11a_overlap__Degree_cindex.disease.pdf",
   ],
   "ED10d": [f"{ORGAN_ED10D}/ED10d"],
 },
 "06_frailty_function": {
   "ED8a": [f"{FRAIL}/ED8a"],
   "ED8b": [f"{FRAIL}/ED8b"],
   "ED8c": [f"{BASELINE}/ED8c__CHARLS.DISCO.Frailty.pdf", f"{BASELINE}/ED8c__CHARLS.DISCO.Frailty2.pdf"],
   "ED8d": [f"{RECONROC}/ED8d_RECONSTRUCTED.pdf"],
 },
 "07_RuLAS_omics": {
   "Fig2c": [f"{RULAS}/Fig2c__RuLAS.Clinical.DISCO.surv.pdf"],
   "Fig2d": [f"{RULAS}/Fig2d__RuLAS.Metabolomics.DISCO.surv.pdf"],
   "Fig2e": [f"{RULAS}/Fig2e__RuLAS.Microbe.DISCO.surv.pdf"],
   "Fig2f": [f"{RULAS}/Fig2f__RuLAS.Aggregate.DISCO.surv.pdf"],
   "ED9f": [f"{RULAS}/ED9f__Subsystems.cindex.pdf"],
 },
}

# code (replay scripts) per module
CODE = {
 "01_method_definition": [
   f"{BASE}/strict/extensions/simulation_geometry/replay_geometry.R",
   f"{BASE}/strict/extensions/simulation_survival/replay_survival.R",
   f"{BASE}/strict/run_isolated.py",
   f"{BASE}/strict/build_all_figures_inventory.py",
 ],
 "02_primary_clinical": [f"{BASE}/strict/extensions/clinical/replay_clinical.R"],
 "03_multimodal_UKB": [
   f"{BASE}/strict/extensions/fig3c_reestimated/recompute.R",
   f"{BASE}/strict/extensions/fig3c_reestimated/plot_forest.R",
 ],
 "04_external_validation": [],
 "05_organ_network": [f"{BASE}/strict/extensions/organ/replay_organ.R"],
 "06_frailty_function": [f"{BASE}/strict/extensions/frailty/replay_frailty.R"],
 "07_RuLAS_omics": [f"{BASE}/strict/extensions/rulas/replay_rulas.R"],
}

report = {}

def copy_path(src, dst):
    if os.path.isdir(src):
        shutil.copytree(src, dst, dirs_exist_ok=True)
    elif os.path.exists(src):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
    else:
        report.setdefault("missing", []).append(src)

for module, figs in MAP.items():
    rep = os.path.join(OUT_BASE, module, "reproduction")
    os.makedirs(os.path.join(rep, "code"), exist_ok=True)
    for fig, srcs in figs.items():
        dst_dir = os.path.join(rep, fig)
        os.makedirs(dst_dir, exist_ok=True)
        for s in srcs:
            name = os.path.basename(s.rstrip("/"))
            copy_path(s, os.path.join(dst_dir, name))
        report.setdefault(module, {})[fig] = "ok"
    for c in CODE.get(module, []):
        copy_path(c, os.path.join(rep, "code", os.path.basename(c)))

print(json.dumps(report, indent=2, ensure_ascii=False))
