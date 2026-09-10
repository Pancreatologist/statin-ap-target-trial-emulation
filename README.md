# Statin initiation and acute pancreatitis risk — target trial emulation

R analysis scripts for an emulated target trial comparing **statin
initiation with no statin initiation** on the **12-month risk of acute
pancreatitis**, within the **All of Us Research Program** cohort restricted to
individuals with hypertriglyceridaemia (HTG). The design uses a 180-day
landmark to remove immortal-time bias, entropy balancing for confounding,
MICE for missing data, and a weighted discrete-time survival model with
all-cause death as a competing risk.

> **Data access.** This repository contains **code only**. The All of Us data
> used (including the frozen cohort, multiple-imputation frames, propensity
> weights, and bootstrap replicates as `.rds`/`.csv` intermediates) are
> **controlled-access** and cannot be distributed under the All of Us data use
> policies. No participant-level data, derived analytical data files, or study
> result numbers are included or embedded in these scripts. Results of running
> these scripts are therefore not reproducible from this repository alone.

## Study design (summary)

| Element | Specification |
|---|---|
| Target trial | Emulation of a target trial, landmark design |
| Data source | All of Us Research Program, CDR v8, 2008–2023 |
| Cohort | Individuals with hypertriglyceridaemia (TG ≥ 500 mg/dL or an HTG diagnosis code); first qualifying evidence is the index date |
| Time zero (T0) | Index date + 180-day grace period, identical for both arms (eliminates immortal-time bias) |
| Exposure | Statin initiation within the grace period vs non-initiation |
| Primary outcome | Acute pancreatitis within 12 months |
| Competing risk | All-cause mortality |
| Missing data | Multiple imputation by chained equations (5 imputations) |
| Confounding | Entropy balancing (ATT, focal = statin), ~50 baseline covariates across 6 domains |
| Outcome model | Weighted discrete-time survival (pooled logistic, arm × natural spline(month, df = 2) interaction; no PH assumption) |
| Uncertainty | 5 imputations × 200 bootstrap replicates = 1000 replicates, weights re-estimated within each replicate |
| Point estimate / CI | Bootstrap median / percentile 2.5th–97.5th |

Upstream exclusions at/before index: prior statin use, prior acute or chronic
pancreatitis, pancreatic cancer, and insufficient (< 365 days) observation.

## Repository layout

Scripts are grouped by analysis stage (naming reflects the original project):

| Stage | Script(s) |
|---|---|
| Cohort freeze | `freeze_cohort.R`, `refreeze_cohort.R` |
| Imputation / PS / balance | `statin_htg_step1to3_v4.R`, `step1_weightit_balance.R`, `statin_htg_covariates_table1.R` |
| Main bootstrap | `step2_bootstrap_v2.R` |
| Tables and figures | `generate_tables_fresh.R`, `regenerate_figures_final_v2.R` |
| Secondary outcome (new-onset DM) | `ap_dm_secondary_final.R` |
| Subgroup analyses | `subgroup_analysis.R`, `fig_subgroup_forest.R` |
| Sensitivity analyses | `revision_T2b_grace_weightit.R` (grace-period), `revision_T2f_icd_only.R` (ICD-only HTG), `fracture_nco_boot.R` (negative control), `s3_fresh_weightit.R` (active comparator) |
| E-value figure | `fig_evalue.R` |
| Concept definitions | `codelists/supp_tg_stageA_concepts.sql` (OMOP concept resolution) |

## Running

Scripts expect to be run from the repository root, with the All of Us
extraction artefacts (frozen cohort, imputation frames, weights, bootstrap
results) present in the working directory. Because these inputs are
controlled-access data, they are **not** provided here; replace the
`<placeholder>` tokens (and any commented anchor checks) with values from
your own approved extraction before running.

## Methodological notes

- **No research conclusion numbers are embedded.** Hard-coded study results
  (sample sizes, event counts, effect estimates) have been replaced with
  placeholders or read dynamically from the summaries the scripts produce.
- Point estimates use the **median** of bootstrap replicates, not the mean.
- Each bootstrap replicate **re-estimates the weights and refits the model**
  so that weighting uncertainty is propagated.
- Entropy balancing is implemented with `WeightIt` (ebal) to ensure
  convergence on the high-dimensional, skewed covariate space.

## License / caveat

These scripts are provided for **methodological transparency** in connection
with a peer-reviewed submission. They are not a clinical tool and do not
confer any diagnostic or treatment recommendation. See the accompanying
manuscript for the disclaimer and interpretation caveats.

## Reproducing / contact

The node environments, R package versions, and exact seeds used for the
submission run are not reproduced here. For methodological questions please
refer to the corresponding manuscript.