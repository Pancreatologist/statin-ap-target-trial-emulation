# ============================================================================
# Generate Tables 1-3 (fresh weightit)
# Table 1: covariate balance (all-50 SMD=0.00)
# Table 2: primary outcomes (all fresh)
# Table 3: sensitivity analyses (two columns: fresh / old-noted)
# ============================================================================
library(mice); library(dplyr); library(tidyr); library(data.table)

# ---- Load data ----
# reader-facing variable labels (investigator 2026-09-05: three table locations;
# label column only, numbers untouched). Label rationale: smoke_ever=1 is
# survey-reported ever-smoking; current status unobservable (curr_smoke all NA),
# so the label states the survey fact only, consistent with the smoking Limitation.
var_labels <- c(smoke_ever = "Ever smoked (survey-reported)")
pretty_var <- function(v) ifelse(v %in% names(var_labels), unname(var_labels[v]), v)

imp      <- readRDS("imp_statin_htg_v4.rds")
frozen   <- readRDS("frozen_cohort_ids.rds")
frozen_ids <- frozen$person_id
ebal_res <- readRDS("weightit_ebal_fresh.rds")
boot_fresh <- readRDS("bootstrap_fresh_weightit.rds")
grace_sm     <- read.csv("grace_fresh_weightit_summary.csv")   # AUTHORITATIVE re-run (weights re-estimated per rep); old grace_weightit_results.rds = LEGACY
fracture_nco <- read.csv("fracture_nco_fresh_summary.csv")     # AUTHORITATIVE re-run (frozen cohort); old fracture_nco_ci.rds = LEGACY

n_st  <- sum(frozen$arm == "STATIN")
n_nu  <- sum(frozen$arm == "NON_USER")
n_total <- n_st + n_nu

# ---- Bootstrap summaries ----
bs <- data.frame(
  point = apply(boot_fresh, 2, median, na.rm = TRUE),
  lo    = apply(boot_fresh, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
  hi    = apply(boot_fresh, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
)
rownames(bs) <- colnames(boot_fresh)

ap_rr   <- bs["rr_ap", "point"];   ap_lo  <- bs["rr_ap", "lo"];   ap_hi  <- bs["rr_ap", "hi"]
d_rr    <- bs["rr_death", "point"]; d_lo  <- bs["rr_death", "lo"];  d_hi  <- bs["rr_death", "hi"]
rr_cs   <- bs["rr_apcs", "point"]; cs_lo  <- bs["rr_apcs", "lo"];  cs_hi  <- bs["rr_apcs", "hi"]

apST12  <- bs["apST_12", "point"]; apNU12 <- bs["apNU_12", "point"]
dST12   <- bs["dST_12", "point"];  dNU12  <- bs["dNU_12", "point"]

# AP RD (12-mo, /100k): point = median of bootstrap replicate RDs (rd_12 == apST_12 - apNU_12), CI = 2.5/97.5 percentiles
ap_rd    <- median(boot_fresh[, "rd_12"], na.rm = TRUE)
ap_rd_lo <- quantile(boot_fresh[, "rd_12"], 0.025, na.rm = TRUE)
ap_rd_hi <- quantile(boot_fresh[, "rd_12"], 0.975, na.rm = TRUE)
# Death RD (12-mo, /100k): same rule from death CIF columns
d_rd    <- median(boot_fresh[, "dST_12"] - boot_fresh[, "dNU_12"], na.rm = TRUE)
d_rd_lo <- quantile(boot_fresh[, "dST_12"] - boot_fresh[, "dNU_12"], 0.025, na.rm = TRUE)
d_rd_hi <- quantile(boot_fresh[, "dST_12"] - boot_fresh[, "dNU_12"], 0.975, na.rm = TRUE)

# Grace (AUTHORITATIVE re-run summary)
g90  <- grace_sm[grace_sm$grace_days == 90, ]
g365 <- grace_sm[grace_sm$grace_days == 365, ]
g90_rr  <- g90$rr_ap_median;  g90_lo  <- g90$rr_ap_ci_lo;  g90_hi  <- g90$rr_ap_ci_hi
g365_rr <- g365$rr_ap_median; g365_lo <- g365$rr_ap_ci_lo; g365_hi <- g365$rr_ap_ci_hi

# Fracture NC (AUTHORITATIVE re-run, frozen cohort)
frac_rr  <- fracture_nco$rr_median
frac_lo  <- fracture_nco$rr_ci_lo
frac_hi  <- fracture_nco$rr_ci_hi

# ICD-only HTG subpopulation (AUTHORITATIVE re-run 2026-09-06, code pathway
# from t6_htg_pathway.csv; frozen subpopulation; 5x200; weights re-estimated per rep)
icd_sm <- read.csv("icd_only_fresh_summary.csv")
icd_rr  <- icd_sm$rr_median; icd_lo <- icd_sm$rr_ci_lo; icd_hi <- icd_sm$rr_ci_hi


# Median follow-up
d1 <- as.data.table(complete(imp, 1))
d1 <- d1[person_id %in% frozen_ids, ]
d1[, time_days := as.numeric(time_days)]
med_fup <- median(d1$time_days, na.rm = TRUE)
med_fup_yr <- med_fup / 365.25
med_fup_mo <- med_fup / 30.44  # B3-6 fix: months = days/30.44 (old footnote used med_fup_yr*30.44 = days/12)

# 12-month events (imp1)
d1[, idx := .I]
d1[, ap_month := ifelse(ap_event == 1 & time_days <= 365, floor(time_days / 30.44), NA_real_)]
d1[, death_month := ifelse(death_event == 1 & time_days <= 365, floor(time_days / 30.44), NA_real_)]
d1 <- d1[time_days > 0, ]
ap12_st <- sum(d1$arm == "STATIN" & !is.na(d1$ap_month) & d1$ap_month < 12, na.rm = TRUE)
ap12_nu <- sum(d1$arm == "NON_USER" & !is.na(d1$ap_month) & d1$ap_month < 12, na.rm = TRUE)
d12_st  <- sum(d1$arm == "STATIN" & !is.na(d1$death_month) & d1$death_month < 12, na.rm = TRUE)
d12_nu  <- sum(d1$arm == "NON_USER" & !is.na(d1$death_month) & d1$death_month < 12, na.rm = TRUE)

# Person-years
py_st <- sum(d1$arm == "STATIN", na.rm = TRUE) / 12
py_nu <- sum(d1$arm == "NON_USER", na.rm = TRUE) / 12
ir_st <- ap12_st / py_st * 1000
ir_nu <- ap12_nu / py_nu * 1000

# ============================================================================
# TABLE 1: Covariate Balance
# ============================================================================
ps_vars <- c("age","male","white","hispanic","smoke_ever","alcohol",
             "ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
             "n_inpatient_ed","n_outpatient","n_visits_total",
             "htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
             "dementia","as_spond","ra","sle","pulmonary","cancer","neuro",
             "aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
             "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
             "antithromb","nsaid")

covar_groups <- list(
  Demographics = c("age","male","white","hispanic","smoke_ever","alcohol"),
  Lipids = c("ldl","hdl","tc","tg"),
  `Vitals/Renal` = c("bmi","calcium","alt","egfr"),
  `Healthcare Utilisation` = c("n_inpatient_ed","n_outpatient","n_visits_total"),
  Comorbidities = c("htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
                      "dementia","as_spond","ra","sle","pulmonary","cancer","neuro"),
  Medications = c("aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
                   "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
                   "antithromb","nsaid")
)

calc_smd <- function(d, w = NULL, vars = ps_vars) {
  trt <- d$arm == "STATIN"; ctrl <- d$arm == "NON_USER"
  r <- data.frame(variable = vars, stringsAsFactors = FALSE)
  r$mean_statin <- sapply(vars, function(v) mean(as.numeric(d[[v]])[trt], na.rm = TRUE))
  r$mean_nonuser <- sapply(vars, function(v) mean(as.numeric(d[[v]])[ctrl], na.rm = TRUE))
  sd_pooled <- sapply(vars, function(v) {
    sqrt((var(as.numeric(d[[v]])[trt], na.rm = TRUE) + var(as.numeric(d[[v]])[ctrl], na.rm = TRUE)) / 2)
  })
  r$smd_before <- abs(r$mean_statin - r$mean_nonuser) / sd_pooled
  if (!is.null(w)) {
    w_c <- w[ctrl]
    r$mean_nonuser_w <- sapply(vars, function(v) {
      weighted.mean(as.numeric(d[[v]])[ctrl], w_c, na.rm = TRUE)
    })
    r$smd_after <- abs(r$mean_statin - r$mean_nonuser_w) / sd_pooled
  }
  r
}

d1_tab <- as.data.table(complete(imp, 1))
d1_tab <- d1_tab[person_id %in% frozen_ids, ]
d1_tab$arm <- factor(d1_tab$arm, levels = c("NON_USER", "STATIN"))
d1_tab$w <- as.numeric(ebal_res$weightit_models[[1]]$weights)
smd_df <- calc_smd(d1_tab, d1_tab$w)
smd_df$group <- "Other"
for (g in names(covar_groups)) smd_df$group[smd_df$variable %in% covar_groups[[g]]] <- g

is_binary <- function(x) {
  u <- unique(na.omit(x))
  length(u) <= 2 && all(u %in% c(0, 1))
}

table1_rows <- list()
for (g in names(covar_groups)) {
  table1_rows[[paste0("g__", g)]] <- data.frame(
    Variable = g, Group = g,
    `Statin (raw)` = "", `Non-user (raw)` = "", `Non-user (weighted)` = "",
    `SMD before` = NA, `SMD after` = NA,
    is_header = TRUE, stringsAsFactors = FALSE
  )
  for (v in covar_groups[[g]]) {
    row <- smd_df[smd_df$variable == v, ]
    x <- as.numeric(d1_tab[[v]])
    x_st <- x[d1_tab$arm == "STATIN"]
    x_nu <- x[d1_tab$arm == "NON_USER"]
    w_nu <- d1_tab$w[d1_tab$arm == "NON_USER"]
    if (is_binary(x)) {
      st_str <- sprintf("%.1f%%", 100 * mean(x_st, na.rm = TRUE))
      nu_str <- sprintf("%.1f%%", 100 * mean(x_nu, na.rm = TRUE))
      nu_w_str <- sprintf("%.1f%%", 100 * weighted.mean(x_nu, w_nu, na.rm = TRUE))
    } else {
      st_str <- sprintf("%.2f (%.2f)", mean(x_st, na.rm = TRUE), sd(x_st, na.rm = TRUE))
      nu_str <- sprintf("%.2f (%.2f)", mean(x_nu, na.rm = TRUE), sd(x_nu, na.rm = TRUE))
      nu_w_str <- sprintf("%.2f", weighted.mean(x_nu, w_nu, na.rm = TRUE))
    }
    table1_rows[[v]] <- data.frame(
      Variable = pretty_var(v), Group = g,
      `Statin (raw)` = st_str, `Non-user (raw)` = nu_str,
      `Non-user (weighted)` = nu_w_str,
      `SMD before` = row$smd_before, `SMD after` = row$smd_after,
      is_header = FALSE, stringsAsFactors = FALSE
    )
  }
}
table1_df <- bind_rows(table1_rows)
n_gt_01 <- sum(abs(smd_df$smd_after) > 0.1, na.rm = TRUE)
max_smd_after <- max(abs(smd_df$smd_after), na.rm = TRUE)
cat(sprintf("Table 1: %d/50 SMD>0.1, max=%.4f\n", n_gt_01, max_smd_after))
write.csv(table1_df, "table1_covariate_balance.csv", row.names = FALSE)
cat("Table 1 saved.\n")

# ============================================================================
# TABLE 2: Primary Outcomes
# ============================================================================
cat("\n=== Table 2 ===\n")

table2 <- data.frame(
  Outcome = c("Acute pancreatitis", "All-cause mortality"),
  `N statin` = c(n_st, n_st),
  `AP events (12-mo) statin` = c(ap12_st, d12_st),
  `CIF 12-mo statin (/100k)` = c(sprintf("%.1f", apST12), sprintf("%.1f", dST12)),
  `N non-user` = c(n_nu, n_nu),
  `AP events (12-mo) non-user` = c(ap12_nu, d12_nu),
  `CIF 12-mo non-user (/100k)` = c(sprintf("%.1f", apNU12), sprintf("%.1f", dNU12)),
  `RD 12-mo (/100k)` = c(
    sprintf("%.1f (%.1f to %.1f)", ap_rd, ap_rd_lo, ap_rd_hi),
    sprintf("%.1f (%.1f to %.1f)", d_rd, d_rd_lo, d_rd_hi)
  ),
  `RR (95% CI)` = c(
    sprintf("%.2f (%.2f to %.2f)", ap_rr, ap_lo, ap_hi),
    sprintf("%.2f (%.2f to %.2f)", d_rr, d_lo, d_hi)
  ),
  stringsAsFactors = FALSE
)

# Median follow-up footnote (B3-6: months fixed to days/30.44; RD/CIF medians
# are computed independently across bootstrap replicates, so RD != diff of
# median CIFs -- stated explicitly per reviewer request)
fup_note <- sprintf(
  paste0("Median follow-up: %.1f months (%.1f years). 12-month window events: ",
         "statin AP=%d, non-user AP=%d; death: statin=%d, non-user=%d. ",
         "RD = median of bootstrap replicate-specific 12-month risk differences ",
         "(per 100,000); 95%% CI from 2.5th/97.5th percentiles. RD and each CIF are ",
         "independently summarised as bootstrap medians, so the RD does not equal ",
         "the difference between the two median CIFs."),
  med_fup_mo, med_fup_yr, ap12_st, ap12_nu, d12_st, d12_nu
)
table2_foot <- data.frame(
  Outcome = fup_note,
  `N statin` = "", `AP events (12-mo) statin` = "",
  `CIF 12-mo statin (/100k)` = "",
  `N non-user` = "",
  `AP events (12-mo) non-user` = "",
  `CIF 12-mo non-user (/100k)` = "",
  `RD 12-mo (/100k)` = "",
  `RR (95% CI)` = "",
  stringsAsFactors = FALSE
)
table2 <- rbind(table2, table2_foot)

write.csv(table2, "table2_primary_outcomes.csv", row.names = FALSE)
print(table2)
cat("\nTable 2 saved.\n")

# ============================================================================
# TABLE 3: Sensitivity Analyses
# Two columns: Fresh weightit | Old weightthem (noted)
# ============================================================================
cat("\n=== Table 3 ===\n")

# ICD-only HTG row RESTORED (2026-09-06): authoritative re-run completed on the
# code-pathway subpopulation (t6_htg_pathway.csv).
# Fibrate row DELETED (author ruling 2026-09-06): analysis cancelled, replaced
# by the descriptive T5 counts table. MVA row: not estimable (2 events).
table3 <- data.frame(
  Analysis = c(
    "Primary (180-day grace, competing-risk)",
    "  Cause-specific (sensitivity)",
    "All-cause mortality (competing event)",
    "Negative control: fracture/fall",
    "Negative control: motor vehicle accident",
    "Grace period: 90 days",
    "Grace period: 365 days",
    "ICD-only HTG source population",
    "No-interaction model (sensitivity)"
  ),
  `Fresh WeightIt` = c(
    sprintf("%.2f (%.2f to %.2f)", ap_rr, ap_lo, ap_hi),
    sprintf("%.2f (%.2f to %.2f)", rr_cs, cs_lo, cs_hi),
    sprintf("%.2f (%.2f to %.2f)", d_rr, d_lo, d_hi),
    sprintf("%.2f (%.2f to %.2f)", frac_rr, frac_lo, frac_hi),
    "Not estimable (2 events)",
    sprintf("%.2f (%.2f to %.2f)", g90_rr, g90_lo, g90_hi),
    sprintf("%.2f (%.2f to %.2f)", g365_rr, g365_lo, g365_hi),
    sprintf("%.2f (%.2f to %.2f)", icd_rr, icd_lo, icd_hi),
    "—"
  ),
  `Notes` = c(
    "Fresh weightit ebal; bootstrap 5x200; weights re-estimated per rep",
    "Fresh weightit ebal; same replicate weights",
    "Fresh weightit ebal; bootstrap 5x200; weights re-estimated per rep",
    "AUTHORITATIVE re-run on frozen cohort; 5x200; weights re-estimated per rep",
    "Only 2 events; CI not estimable",
    "AUTHORITATIVE re-run; weights re-estimated per rep; competing-risk AJ CIF",
    "AUTHORITATIVE re-run; weights re-estimated per rep; competing-risk AJ CIF",
    sprintf("Code-pathway subpopulation N=%s (STATIN %s / NON_USER %s); 12-mo events %s/%s; 5x200; weights re-estimated per rep", icd_sm$n_total, icd_sm$n_statin, icd_sm$n_nonuser, icd_sm$ap_events_12mo_statin, icd_sm$ap_events_12mo_nonuser),
    "See Methods sensitivity"
  ),
  stringsAsFactors = FALSE
)

write.csv(table3, "table3_sensitivity.csv", row.names = FALSE)
print(table3)
cat("\nTable 3 saved.\n")
cat("\n=== All tables done ===\n")
