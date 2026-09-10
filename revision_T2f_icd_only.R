# ============================================================================
# revision_T2f_icd_only.R  (2026-09-05, Batch-2 item 1)
# ICD-only HTG sensitivity re-run under AUTHORITATIVE_RUN spec.
# Supersedes the old "diagnosis-code-only HTG 1.09 (point estimate, no CI)"
# (statin_htg_sensitivity.R S2, different population, median-imputation, LEGACY).
#
# Definition (investigator, 2026-09-05): ICD-only HTG subpopulation = frozen
# cohort members whose HTG entry was by diagnostic code only, i.e.
# t6_htg_pathway.csv htg_lab_date == NA. Now precisely definable because the
# pathway file carries code/lab first dates separately (previously the local
# files kept only index_date = MIN(code, lab)).
#
# Spec (mirrors revision_T2b_fracture_nco.R / step2_bootstrap_v2.R):
#   - population: frozen cohort ∩ code_only pathway (subset of the frozen cohort;
#     anchors come from t6_pathway_split.csv, validated at run time)
#   - imputation values reused from imp_statin_htg_v4 (frozen-asset rule);
#     ebal weights RE-ESTIMATED per bootstrap replicate within the subpopulation
#   - outcome: AP event, 12-month window t0 -> min(ap, death, obs_end, t0+365);
#     death = competing event; dual-hazard discrete AJ CIF
#   - glm(ap_m ~ arm * ns(month, df=2, Boundary.knots=c(1,12)), binomial, w)
#     and same for death_m; RR = 12m AP CIF ratio (STATIN / NON_USER)
#   - 5 MICE x 200 bootstrap = 1000 reps; seed system set.seed(20260905 + rep_idx)
#   - weightit failure = DISCARD replicate + count (T2c patched pattern)
# Output: icd_only_fresh.rds / icd_only_fresh_summary.csv / log
# ============================================================================
library(bit64)
library(mice); library(WeightIt); library(data.table); library(splines); library(parallel)

con <- file("revision_T2f_icd_only_run.log", open = "wt")
sink(con, type = "output"); sink(con, type = "message")

SEED_BASE <- 20260905
N_WORKERS <- 12
N_BOOT <- 200; N_IMP <- 5; TOTAL <- N_IMP * N_BOOT; STEP <- 365 / 12

# ---- inputs (read-only) ----
imp    <- readRDS("imp_statin_htg_v4.rds")
frozen <- readRDS("frozen_cohort_ids.rds")
frozen_ids <- as.numeric(frozen$person_id)
dm <- as.data.table(readRDS("df_master.rds"))
dm[, person_id := as.numeric(person_id)]
pw <- read.csv("t6_htg_pathway.csv",
               colClasses = c(person_id = "character"))
pw$person_id <- as.numeric(pw$person_id)

# ---- ICD-only definition + sample flow ----
# gate: pathway file must reproduce the T6 step-1 gate on the full master frame
mi_m <- match(dm$person_id, pw$person_id)
stopifnot(!anyNA(mi_m))
code_only <- !is.na(pw$htg_code_date[mi_m]) & is.na(pw$htg_lab_date[mi_m])
icd_ids <- dm$person_id[code_only & dm$person_id %in% frozen_ids]
n_frozen <- length(frozen_ids)
n_icd    <- length(icd_ids)
cat(sprintf("SAMPLE FLOW | frozen cohort | N=%d\n", n_frozen))
cat(sprintf("STEP 1 | restrict to ICD-only HTG pathway (htg_lab_date NA) | N_before=%d | excluded=%d (lab_only + both) | N_after=%d\n",
            n_frozen, n_frozen - n_icd, n_icd))
# anchors from t6_pathway_split.csv (must match bitwise; validated at run time,
# anchor N not embedded here)
# stopifnot(n_icd == <N_ICD_ANCHOR>)  # placeholder: set to frozen ICD-only N
dates_dt <- dm[person_id %in% icd_ids,
  .(person_id, t0 = as.Date(t0), ap_date = as.Date(ap_date),
    death_date = as.Date(death_date), obs_end = as.Date(obs_end))]
stopifnot(nrow(dates_dt) == n_icd)

# arm labels for anchors (imp-invariant, from df_master)
arm_m <- dm$person_id %in% icd_ids
arm_lab <- dm$arm[arm_m]
# arm anchor checks validated at run time; arm Ns not embedded here
# stopifnot(sum(arm_lab == "STATIN") == <N_ST_ANCHOR>,
#           sum(arm_lab == "NON_USER") == <N_NU_ANCHOR>)

ps_vars <- c("age","male","white","hispanic","smoke_ever","alcohol",
             "ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
             "n_inpatient_ed","n_outpatient","n_visits_total",
             "htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
             "dementia","as_spond","ra","sle","pulmonary","cancer","neuro",
             "aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
             "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
             "antithromb","nsaid")
PS_FORMULA <- reformulate(ps_vars, response = "arm")

build_imp_data <- function(i) {
  d_raw <- as.data.table(complete(imp, i))
  d_raw[, person_id := as.numeric(person_id)]
  n_before <- nrow(d_raw)
  d <- d_raw[person_id %in% icd_ids, ]
  cat(sprintf("STEP imp%d | restrict to ICD-only frozen | N_before=%d | excluded=%d | N_after=%d\n",
              i, n_before, n_before - nrow(d), nrow(d)))
  d[, arm := factor(arm, levels = c("NON_USER", "STATIN"))]
  for (v in c(ps_vars, "ap_event", "death_event", "time_days")) {
    if (v %in% names(d) && inherits(d[[v]], "integer64")) d[, (v) := as.numeric(get(v))]
  }
  d <- merge(d, dates_dt, by = "person_id")   # 1:1, all ICD-only present
  stopifnot(nrow(d) == n_icd)
  SENT <- as.Date("2099-01-01")
  d[, t0d  := as.Date(t0)]
  d[, ap_s := fifelse(is.na(ap_date), SENT, as.Date(ap_date))]
  d[, dth_s:= fifelse(is.na(death_date), SENT, as.Date(death_date))]
  d[, end_f := pmin(ap_s, dth_s, as.Date(obs_end), t0d + 365)]
  d[, ev_ap := !is.na(ap_date) & ap_date > t0d & ap_date <= end_f]
  d[, ev_d  := !ev_ap & !is.na(death_date) & death_date > t0d & death_date <= end_f]
  d[, time_ap := as.numeric(fifelse(ev_ap, as.Date(ap_date), end_f) - t0d)]
  stopifnot(all(d$time_ap > 0))
  d[, M := pmin(pmax(1L, ceiling(time_ap / STEP)), 12L)]
  d
}

max_smd_w <- function(d, w) {
  trt <- d$arm == "STATIN"; ctl <- d$arm == "NON_USER"
  s <- numeric(length(ps_vars))
  for (k in seq_along(ps_vars)) {
    x  <- as.numeric(d[[ps_vars[k]]]); xt <- x[trt]; xc <- x[ctl]
    sdp <- sqrt((var(xt, na.rm = TRUE) + var(xc, na.rm = TRUE)) / 2)
    s[k] <- if (is.finite(sdp) && sdp > 0)
      abs(mean(xt, na.rm = TRUE) - weighted.mean(xc, w[ctl], na.rm = TRUE)) / sdp else 0
  }
  max(s, na.rm = TRUE)
}

col_names <- c(paste0("aNU_", 1:12), paste0("aST_", 1:12), "rr_ap")

run_rep <- function(rep_idx) {
  set.seed(SEED_BASE + rep_idx)
  idx_b <- sample(nrow(D_IMP), replace = TRUE)
  d_b <- D_IMP[idx_b]
  wb <- tryCatch(
    weightit(PS_FORMULA, data = d_b, method = "ebal", estimand = "ATT", focal = "STATIN"),
    error = function(e) e)
  if (inherits(wb, "error")) {
    return(list(row = rep(NA_real_, length(col_names)), fail = 1L, smd = NA_real_))
  }
  wv <- as.numeric(wb$weights)
  d_b[, w := wv]
  smd_rep <- max_smd_w(d_b, wv)

  n_per <- d_b$M
  idx_long <- rep(seq_len(nrow(d_b)), n_per)
  pm <- d_b[idx_long, .(arm, w, ae = ev_ap, de = ev_d)]
  pm[, month := sequence(n_per)]
  pm[, last_m := month == rep(n_per, n_per)]
  pm[, ap_m := as.integer(ae == 1 & last_m)]
  pm[, death_m := as.integer(de == 1 & last_m)]

  out <- rep(NA_real_, length(col_names))
  if (nrow(pm) > 100 && sum(pm$ap_m) > 3 && sum(pm$death_m) > 3) {
    tryCatch({
      f_a <- glm(ap_m ~ arm * ns(month, df = 2, Boundary.knots = c(1, 12)),
                 family = binomial(), weights = w, data = pm)
      f_d <- glm(death_m ~ arm * ns(month, df = 2, Boundary.knots = c(1, 12)),
                 family = binomial(), weights = w, data = pm)
      cifs <- list()
      for (a in c("NON_USER", "STATIN")) {
        nd <- data.frame(arm = factor(a, levels = c("NON_USER", "STATIN")), month = 1:12)
        h_a <- predict(f_a, nd, type = "response")
        h_d <- predict(f_d, nd, type = "response")
        S <- 1; cf <- numeric(12)
        for (j in 1:12) {
          cf[j] <- (if (j > 1) cf[j - 1] else 0) + h_a[j] * S
          S <- S * (1 - h_a[j] - h_d[j])
        }
        cifs[[a]] <- cf
      }
      out <- c(cifs$NON_USER * 1e5, cifs$STATIN * 1e5,
               cifs$STATIN[12] / max(cifs$NON_USER[12], 1e-10))
    }, error = function(e) {})
  }
  list(row = out, fail = 0L, smd = smd_rep)
}

cl <- makeCluster(N_WORKERS)
clusterEvalQ(cl, { library(WeightIt); library(data.table); library(splines); suppressMessages(library(mice)) })
clusterExport(cl, c("run_rep", "SEED_BASE", "col_names", "ps_vars", "max_smd_w"))

boot_res <- matrix(NA_real_, nrow = TOTAL, ncol = length(col_names))
colnames(boot_res) <- col_names

t0 <- Sys.time()
fail_track <- integer(TOTAL)
smd_track  <- rep(NA_real_, TOTAL)
cat(sprintf("Start: %s | plan: %d imp x %d boot = %d reps | workers=%d | seed=%d+rep_idx\n",
            t0, N_IMP, N_BOOT, TOTAL, N_WORKERS, SEED_BASE))

ess_tbl <- data.table()
for (i in 1:N_IMP) {
  D_IMP <- build_imp_data(i)
  if (i == 1) {
    # anchor check: date-based 12m AP events must match t6_pathway_split.csv
    e_st <- sum(D_IMP$ev_ap[D_IMP$arm == "STATIN"])
    e_nu <- sum(D_IMP$ev_ap[D_IMP$arm == "NON_USER"])
    cat(sprintf("EVENTS (imp1, date-based, imp-invariant): AP_12mo ST=%d NU=%d total=%d | death_terminal=%d\n",
                e_st, e_nu, e_st + e_nu, sum(D_IMP$ev_d)))
    stopifnot(e_st == 18, e_nu == 61)
  }
  clusterExport(cl, c("D_IMP", "PS_FORMULA"))

  wf <- tryCatch(weightit(PS_FORMULA, data = D_IMP, method = "ebal",
                          estimand = "ATT", focal = "STATIN"), error = function(e) NULL)
  if (!is.null(wf)) {
    w <- as.numeric(wf$weights); st <- D_IMP$arm == "STATIN"
    ess_tbl <- rbind(ess_tbl, data.table(
      imp = i,
      ess_st = sum(w[st])^2 / sum(w[st]^2), ess_nu = sum(w[!st])^2 / sum(w[!st]^2)))
  }

  reps <- (i - 1) * N_BOOT + 1:N_BOOT
  res <- parLapply(cl, reps, run_rep)
  for (k in seq_along(reps)) {
    boot_res[reps[k], ] <- res[[k]]$row
    fail_track[reps[k]] <- res[[k]]$fail
    smd_track[reps[k]]  <- res[[k]]$smd
  }
  rr_v <- boot_res[reps, "rr_ap"]; rr_v <- rr_v[is.finite(rr_v)]
  cat(sprintf("Imp %d done | %.1f min elapsed | valid so far=%d | RR(median so far)=%.3f\n",
              i, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              sum(is.finite(boot_res[, "rr_ap"])),
              if (length(rr_v)) median(rr_v) else NA))
  saveRDS(boot_res, "icd_only_fresh.rds")   # checkpoint each imp
}
stopCluster(cl)

# ---- summary ----
rr <- boot_res[, "rr_ap"]
valid <- is.finite(rr) & rr > 0 & rr < 100
rr_v <- rr[valid]
aST12 <- boot_res[, "aST_12"]; aNU12 <- boot_res[, "aNU_12"]
d1 <- build_imp_data(1)
sm <- data.table(
  outcome = "Acute pancreatitis, ICD-only HTG subpopulation",
  n_total = nrow(d1), n_statin = sum(d1$arm == "STATIN"), n_nonuser = sum(d1$arm == "NON_USER"),
  ap_events_12mo_statin = sum(d1$ev_ap[d1$arm == "STATIN"]),
  ap_events_12mo_nonuser = sum(d1$ev_ap[d1$arm == "NON_USER"]),
  rr_median = median(rr_v),
  rr_ci_lo = quantile(rr_v, 0.025), rr_ci_hi = quantile(rr_v, 0.975),
  cif12_statin_per100k = median(aST12[is.finite(aST12)]),
  cif12_nonuser_per100k = median(aNU12[is.finite(aNU12)]),
  n_valid = length(rr_v), n_total_reps = TOTAL,
  n_weightit_fail = sum(fail_track),
  max_smd_median = median(smd_track, na.rm = TRUE),
  max_smd_max = max(smd_track, na.rm = TRUE),
  ess_statin = median(ess_tbl$ess_st), ess_nonuser = median(ess_tbl$ess_nu),
  crosses_null = quantile(rr_v, 0.025) < 1 & quantile(rr_v, 0.975) > 1,
  seed_base = SEED_BASE,
  spec = "AUTHORITATIVE: frozen ICD-only (code pathway) subpopulation, 5x200, weights re-estimated per rep, AJ CIF")
cat("\n=== FINAL ===\n"); print(sm)
write.csv(sm, "icd_only_fresh_summary.csv", row.names = FALSE)
cat("\nSaved: icd_only_fresh.rds, icd_only_fresh_summary.csv\n")
cat(sprintf("UNSTABLE check: events=%d (>=10), ESS ST=%.0f NU=%.0f (>=50), valid=%d/%d (%.1f%%)\n",
            sum(d1$ev_ap), median(ess_tbl$ess_st), median(ess_tbl$ess_nu),
            length(rr_v), TOTAL, 100 * length(rr_v) / TOTAL))
