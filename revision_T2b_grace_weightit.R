# ============================================================================
# revision_T2b_grace_weightit.R
# Grace-period sensitivity re-run (GRACE = 90 / 365 days) under AUTHORITATIVE_RUN
# spec. Supersedes grace_weightit_bootstrap.R (2026-06-28) -> that file and its
# output grace_weightit_results.rds remain untouched and are LEGACY.
#
# Machinery defects fixed vs old script (DESIGN unchanged, see below):
#   1. weights: ebal weights RE-ESTIMATED inside every bootstrap replicate
#      (old: full-data weights held fixed across replicates -> weight-estimation
#      uncertainty not propagated, CIs too narrow)
#   2. competing risk: death modelled via dual-hazard discrete Aalen-Johansen
#      CIF (old: single-hazard net CIF, death merely censored)
#   3. arm x time interaction restored in both hazard glms:
#      glm(y ~ arm * ns(month, df = 2, Boundary.knots = c(1, 12)))
#      (old: additive arm + spline only)
#   4. person-month structure: follow-up truncated at first event
#      (old: AP cases kept contributing at-risk rows after their event month)
#   5. seed system: set.seed(20260625 + rep_idx), rep_idx = 1..1000 imp-major
#      order, parallel-safe (old: 20260628 + imp index, sequential global stream)
#
# Design preserved from original grace script (= what reviewers saw):
#   - base population: full MICE imputation frame (size not embedded), NOT
#     frozen-restricted -- the frozen cohort embeds the 180-day grace
#     exclusions; each GR rebuilds its own landmark cohort with
#     variant-consistent drops (this is the point of the sensitivity)
#   - dates from df_final_statin_htg_v4.rds (the imputation source frame),
#     attached via match() to preserve mids row order
#   - arm_GR: first statin exposure in [index_date, index_date + GR) -> STATIN
#   - t0 = index_date + GR; drop AP/death before t0; follow-up = t0 to
#     min(AP, death, obs_end, t0 + 365); keep time > 0; 12-month window
#   - AP = event, death = competing event (first-event-wins ordering)
#   - RR = 12-month AP CIF ratio (STATIN / NON_USER)
#
# GATE (before any bootstrap): GR=180 rebuild on imp1 must reproduce the frozen
#   cohort exactly: same id set, arm == frozen arm, 12-month event
#   flags == imp-frame flags with 0 disagreement. Machinery itself validated in
#   revision_T2b_grace_preflight.R: rep1/rep2 rr_ap bitwise identical to
#   bootstrap_fresh_weightit.rds (abs diff 0.000e+00); date-derived events,
#   capped times and month indices all agree with imp-frame values exactly.
#
# Parallel: 12 workers over replicates (per-rep seeded index draw makes
#   parallel == serial). Machine: 16 logical cores.
# Output (NEW files): grace_fresh_weightit.rds / grace_fresh_weightit_summary.csv
# ============================================================================
library(bit64)   # MUST load before any as.numeric() on integer64 ids
library(mice); library(WeightIt); library(data.table); library(splines); library(parallel)

con <- file("revision_T2b_grace_weightit_run.log", open = "wt")
sink(con, type = "output"); sink(con, type = "message")

t_start <- Sys.time()
SEED_BASE <- 20260625
N_WORKERS <- 12
N_BOOT <- 200; N_IMP <- 5; TOTAL <- N_IMP * N_BOOT; STEP <- 365 / 12
GRACE_SET <- c(90L, 365L)

imp    <- readRDS("imp_statin_htg_v4.rds")
frozen <- readRDS("frozen_cohort_ids.rds")
frozen_ids <- as.numeric(frozen$person_id)
df_v4  <- readRDS("df_final_statin_htg_v4.rds")

ps_vars <- c("age","male","white","hispanic","smoke_ever","alcohol",
             "ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
             "n_inpatient_ed","n_outpatient","n_visits_total",
             "htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
             "dementia","as_spond","ra","sle","pulmonary","cancer","neuro",
             "aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
             "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
             "antithromb","nsaid")
PS_FORMULA <- reformulate(ps_vars, response = "arm")

# ---- date source: df_v4 (imputation source frame) ----
dv <- as.data.table(df_v4)
dv[, person_id := as.numeric(person_id)]
dates_dt <- dv[, .(person_id,
                   index_date   = as.Date(index_date),
                   first_statin = as.Date(first_statin),
                   ap_date      = as.Date(ap_date),
                   death_date   = as.Date(death_date),
                   obs_end      = as.Date(obs_end))]
stopifnot(!anyNA(dates_dt$obs_end), !any(duplicated(dates_dt$person_id)))
imp_ids <- as.numeric(imp$data$person_id)
stopifnot(all(imp_ids %in% dates_dt$person_id))   # imp frame fully covered 1:1
fr_arm_dt <- as.data.table(frozen)[, .(person_id = as.numeric(person_id),
                                       arm_frozen = as.character(arm))]

# ---- per-GR landmark rebuild (rules validated in preflight at GR=180) ----
build_grace <- function(i, GR) {
  d <- as.data.table(complete(imp, i))
  d[, person_id := as.numeric(person_id)]
  # match-attach preserves mids row order (merge() would reorder; see preflight)
  mi <- match(d$person_id, dates_dt$person_id)
  stopifnot(!anyNA(mi))
  d[, index_date   := dates_dt$index_date[mi]]
  d[, first_statin := dates_dt$first_statin[mi]]
  d[, ap_date      := dates_dt$ap_date[mi]]
  d[, death_date   := dates_dt$death_date[mi]]
  d[, obs_end      := dates_dt$obs_end[mi]]
  for (v in c(ps_vars, "ap_event", "death_event", "time_days")) {
    if (v %in% names(d) && inherits(d[[v]], "integer64")) d[, (v) := as.numeric(get(v))]
  }
  d[, t0 := index_date + GR]
  d[, arm := factor(
       ifelse(!is.na(first_statin) & first_statin >= index_date & first_statin < t0,
              "STATIN", "NON_USER"), levels = c("NON_USER", "STATIN"))]
  n0 <- nrow(d)
  d <- d[is.na(ap_date) | ap_date >= t0]
  cat(sprintf("STEP GR%d imp%d | drop AP before t0    | %d -> %d (-%d)\n",
              GR, i, n0, nrow(d), n0 - nrow(d)))
  n0 <- nrow(d)
  d <- d[is.na(death_date) | death_date >= t0]
  cat(sprintf("STEP GR%d imp%d | drop death before t0 | %d -> %d (-%d)\n",
              GR, i, n0, nrow(d), n0 - nrow(d)))
  # 12-month window from t0; first-event-wins (AP = event, death = competing)
  SENT <- as.Date("2099-01-01")
  d[, ap_s  := fifelse(is.na(ap_date),     SENT, ap_date)]
  d[, dth_s := fifelse(is.na(death_date),  SENT, death_date)]
  d[, end_f := pmin(ap_s, dth_s, obs_end, t0 + 365)]
  d[, ev_ap := as.integer(!is.na(ap_date) & ap_date > t0 & ap_date <= end_f)]
  d[, ev_d  := as.integer(!ev_ap & !is.na(death_date) & death_date > t0 & death_date <= end_f)]
  d[, time  := as.numeric(end_f - t0)]
  n0 <- nrow(d)
  d <- d[time > 0]
  cat(sprintf("STEP GR%d imp%d | drop time<=0        | %d -> %d (-%d)\n",
              GR, i, n0, nrow(d), n0 - nrow(d)))
  stopifnot(all(d$time > 0))
  d[, M := pmin(pmax(1L, ceiling(time / STEP)), 12L)]
  d
}

# ---- GATE: GR=180 rebuild (imp1) must reproduce the frozen cohort exactly ----
cat("\n=== GATE: GR=180 rebuild vs frozen cohort (imp1) ===\n")
g180 <- build_grace(1, 180L)
cat(sprintf("GR180 N=%d | frozen N=%d | same id set: %s\n",
            nrow(g180), length(frozen_ids), setequal(g180$person_id, frozen_ids)))
stopifnot(nrow(g180) == length(frozen_ids), setequal(g180$person_id, frozen_ids))
mi_f <- match(g180$person_id, fr_arm_dt$person_id)
cat(sprintf("arm agreement vs frozen: %.4f | STATIN n=%d\n",
            mean(g180$arm == fr_arm_dt$arm_frozen[mi_f]), sum(g180$arm == "STATIN")))
stopifnot(all(g180$arm == fr_arm_dt$arm_frozen[mi_f]))
g180[, ev_ap_imp := as.integer(ap_event == 1 & time_days <= 365)]
g180[, ev_d_imp  := as.integer(death_event == 1 & time_days <= 365)]
cat(sprintf("AP 12mo: date-derived=%d imp-frame=%d disagree=%d\n",
            sum(g180$ev_ap), sum(g180$ev_ap_imp), sum(g180$ev_ap != g180$ev_ap_imp)))
cat(sprintf("Death 12mo: date-derived=%d imp-frame=%d disagree=%d\n",
            sum(g180$ev_d), sum(g180$ev_d_imp), sum(g180$ev_d != g180$ev_d_imp)))
stopifnot(sum(g180$ev_ap != g180$ev_ap_imp) == 0,
          sum(g180$ev_d  != g180$ev_d_imp)  == 0)
cat("GATE PASSED\n")
rm(g180); gc()

# ---- one replicate: seeded index draw -> weightit -> person-month -> AJ CIF ----
# PATCH 2026-09-05 (T2c): weightit failure = DISCARD replicate + count.
# The original w=1 fallback silently mixed unweighted replicates into the
# bootstrap distribution. run_rep now returns list(row, fail, smd).
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

run_rep <- function(rep_idx) {
  set.seed(SEED_BASE + rep_idx)
  idx_b <- sample(nrow(D_IMP), replace = TRUE)
  d_b <- D_IMP[idx_b]
  wb <- tryCatch(
    weightit(PS_FORMULA, data = d_b, method = "ebal", estimand = "ATT", focal = "STATIN"),
    error = function(e) e)
  if (inherits(wb, "error")) {
    return(list(row = rep(NA_real_, 50), fail = 1L, smd = NA_real_))
  }
  wv <- as.numeric(wb$weights)
  d_b[, w := wv]
  smd_rep <- max_smd_w(d_b, wv)

  n_per <- d_b$M
  idx_long <- rep(seq_len(nrow(d_b)), n_per)
  pm <- d_b[idx_long, .(arm, w, ev_ap, ev_d)]
  pm[, month := sequence(n_per)]
  pm[, last_m := month == rep(n_per, n_per)]
  pm[, ap_m := as.integer(ev_ap == 1 & last_m)]
  pm[, death_m := as.integer(ev_d == 1 & last_m)]

  out <- c(rep(NA_real_, 48), NA_real_, NA_real_)
  if (nrow(pm) > 100 && sum(pm$ap_m) > 3 && sum(pm$death_m) > 3) {
    tryCatch({
      f_ap <- glm(ap_m ~ arm * ns(month, df = 2, Boundary.knots = c(1, 12)),
                  family = binomial(), weights = w, data = pm)
      f_d  <- glm(death_m ~ arm * ns(month, df = 2, Boundary.knots = c(1, 12)),
                  family = binomial(), weights = w, data = pm)
      cf <- list()
      for (a in c("NON_USER", "STATIN")) {
        nd <- data.frame(arm = factor(a, levels = c("NON_USER", "STATIN")), month = 1:12)
        h_ap <- predict(f_ap, nd, type = "response")
        h_d  <- predict(f_d,  nd, type = "response")
        S <- 1; cfa <- numeric(12); cfd <- numeric(12)
        for (j in 1:12) {
          cfa[j] <- (if (j > 1) cfa[j - 1] else 0) + h_ap[j] * S
          cfd[j] <- (if (j > 1) cfd[j - 1] else 0) + h_d[j] * S
          S <- S * (1 - h_ap[j] - h_d[j])
        }
        cf[[a]] <- list(ap = cfa, d = cfd)
      }
      out <- c(cf$NON_USER$ap * 1e5, cf$STATIN$ap * 1e5,
               cf$NON_USER$d * 1e5,  cf$STATIN$d * 1e5,
               cf$STATIN$ap[12] / max(cf$NON_USER$ap[12], 1e-10),
               cf$STATIN$d[12]  / max(cf$NON_USER$d[12],  1e-10))
    }, error = function(e) {})
  }
  list(row = out, fail = 0L, smd = smd_rep)   # PATCH T2c: return with fail flag + max SMD
}

col_names <- c(paste0("apNU_", 1:12), paste0("apST_", 1:12),
               paste0("dNU_", 1:12),  paste0("dST_", 1:12), "rr_ap", "rr_death")
stopifnot(length(col_names) == 50)

cl <- makeCluster(N_WORKERS)
clusterEvalQ(cl, { library(bit64); library(WeightIt); library(data.table); library(splines) })
clusterExport(cl, c("run_rep", "SEED_BASE", "col_names", "ps_vars", "max_smd_w"))

cat(sprintf("\nStart: %s | plan: grace %s x %d imp x %d boot | workers=%d | seed=%d+rep_idx\n",
            t_start, paste(GRACE_SET, collapse = "/"), N_IMP, N_BOOT, N_WORKERS, SEED_BASE))

res_all <- list()
fail_track <- integer(TOTAL)             # PATCH T2c: per-rep weightit failure flag
smd_track  <- rep(NA_real_, TOTAL)       # PATCH T2c: per-rep max |SMD| after weighting
for (GR in GRACE_SET) {
  cat(sprintf("\n========== GRACE = %d days ==========\n", GR))
  boot_res <- matrix(NA_real_, nrow = TOTAL, ncol = 50, dimnames = list(NULL, col_names))
  fail_track[] <- 0L; smd_track[] <- NA_real_    # reset per GR setting
  ess_tbl <- data.table(); info1 <- NULL
  for (i in 1:N_IMP) {
    D_IMP <- build_grace(i, GR)
    if (i == 1) {
      info1 <- list(
        n = nrow(D_IMP), n_st = sum(D_IMP$arm == "STATIN"),
        n_nu = sum(D_IMP$arm == "NON_USER"),
        ap_st = sum(D_IMP$ev_ap[D_IMP$arm == "STATIN"]),
        ap_nu = sum(D_IMP$ev_ap[D_IMP$arm == "NON_USER"]),
        d_st  = sum(D_IMP$ev_d[D_IMP$arm == "STATIN"]),
        d_nu  = sum(D_IMP$ev_d[D_IMP$arm == "NON_USER"]))
      cat(sprintf("EVENTS GR%d (date-based, imp-invariant): AP12mo ST=%d NU=%d | death12mo ST=%d NU=%d\n",
                  GR, info1$ap_st, info1$ap_nu, info1$d_st, info1$d_nu))
    } else {
      stopifnot(nrow(D_IMP) == info1$n)   # date-based build is imp-invariant
    }
    # full-data weights once per imp: ESS record (point-estimate context only)
    wf <- tryCatch(weightit(PS_FORMULA, data = D_IMP, method = "ebal",
                            estimand = "ATT", focal = "STATIN"), error = function(e) NULL)
    if (!is.null(wf)) {
      w <- as.numeric(wf$weights); st <- D_IMP$arm == "STATIN"
      ess_tbl <- rbind(ess_tbl, data.table(
        imp = i, ess_st = sum(w[st])^2 / sum(w[st]^2),
        ess_nu = sum(w[!st])^2 / sum(w[!st]^2)))
    } else cat(sprintf("WARNING: full-data weightit failed at GR=%d imp=%d\n", GR, i))
    clusterExport(cl, c("D_IMP", "PS_FORMULA"))
    reps <- (i - 1) * N_BOOT + 1:N_BOOT
    res <- parLapply(cl, reps, run_rep)
    for (k in seq_along(reps)) {
      boot_res[reps[k], ] <- res[[k]]$row            # PATCH T2c
      fail_track[reps[k]] <- res[[k]]$fail
      smd_track[reps[k]]  <- res[[k]]$smd
    }
    rr_v <- boot_res[reps, "rr_ap"]; rr_v <- rr_v[is.finite(rr_v)]
    cat(sprintf("GR%d imp%d done | %.1f min elapsed | valid so far=%d | RR_ap median so far=%.3f\n",
                GR, i, as.numeric(difftime(Sys.time(), t_start, units = "mins")),
                sum(is.finite(boot_res[, "rr_ap"])),
                if (length(rr_v)) median(rr_v) else NA))
    partial <- res_all
    partial[[as.character(GR)]] <- list(boot_res = boot_res, ess = ess_tbl, info1 = info1)
    saveRDS(partial, "grace_fresh_weightit.rds")   # checkpoint
  }
  res_all[[as.character(GR)]] <- list(boot_res = boot_res, ess = ess_tbl, info1 = info1,
                                      fail = sum(fail_track),                   # PATCH T2c
                                      smd_med = median(smd_track, na.rm = TRUE),
                                      smd_max = max(smd_track, na.rm = TRUE))
  gc()
}
stopCluster(cl)

# ---- summary ----
cat("\n=== SUMMARY ===\n")
sm <- rbindlist(lapply(GRACE_SET, function(GR) {
  r <- res_all[[as.character(GR)]]
  rr <- r$boot_res[, "rr_ap"]; rr_v <- rr[is.finite(rr) & rr > 0 & rr < 100]
  if (!length(rr_v)) rr_v <- NA_real_
  rrd <- r$boot_res[, "rr_death"]; rrd_v <- rrd[is.finite(rrd) & rrd > 0 & rrd < 100]
  if (!length(rrd_v)) rrd_v <- NA_real_
  fST <- r$boot_res[, "apST_12"]; fNU <- r$boot_res[, "apNU_12"]
  lo <- quantile(rr_v, 0.025); hi <- quantile(rr_v, 0.975)
  data.table(
    grace_days = GR,
    n_total = r$info1$n, n_statin = r$info1$n_st, n_nonuser = r$info1$n_nu,
    ap_events_12mo_statin = r$info1$ap_st, ap_events_12mo_nonuser = r$info1$ap_nu,
    death_events_12mo_statin = r$info1$d_st, death_events_12mo_nonuser = r$info1$d_nu,
    rr_ap_median = median(rr_v), rr_ap_ci_lo = lo, rr_ap_ci_hi = hi,
    rr_death_median = median(rrd_v),
    rr_death_ci_lo = quantile(rrd_v, 0.025), rr_death_ci_hi = quantile(rrd_v, 0.975),
    cif12_ap_statin_per100k = median(fST[is.finite(fST)]),
    cif12_ap_nonuser_per100k = median(fNU[is.finite(fNU)]),
    n_valid = sum(is.finite(rr) & rr > 0 & rr < 100), n_total_reps = TOTAL,
    n_weightit_fail = r$fail,                                # PATCH T2c (stored per GR)
    max_smd_median = r$smd_med,                              # PATCH T2c
    max_smd_max = r$smd_max,                                 # PATCH T2c
    ess_statin = median(r$ess$ess_st), ess_nonuser = median(r$ess$ess_nu),
    crosses_null = lo < 1 & hi > 1,
    stable = (r$info1$ap_st >= 10 & r$info1$ap_nu >= 10 &
              median(r$ess$ess_st) >= 50 & median(r$ess$ess_nu) >= 50 &
              sum(is.finite(rr) & rr > 0 & rr < 100) / TOTAL >= 0.95),
    seed_base = SEED_BASE,
    spec = "AUTHORITATIVE: full imp-frame base, per-GR landmark rebuild, 5x200, weights re-estimated per rep, AJ CIF with death as competing event")
}))
print(sm)
fwrite(sm, "grace_fresh_weightit_summary.csv")
saveRDS(res_all, "grace_fresh_weightit.rds")
cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

sink(type = "message"); sink(type = "output"); close(con)
