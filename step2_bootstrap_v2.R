# ============================================================================
# Step 2 (FIXED): Full bootstrap with fresh weightit ebal + interaction model
# RR = 12-month CIF ratio (NOT model coefficient)
# Month: 1-based (1 to 12), as in original v4 implementation
# Seed: 20260625 (added for reproducibility; this run supersedes the no-seed version)
# ============================================================================
library(mice); library(WeightIt); library(dplyr); library(data.table); library(splines)

set.seed(20260625)   # global seed; sample() in bootstrap is the only stochastic step

# sink log to file (avoid PowerShell pipe buffer blocking on long runs)
con <- file("step2_bootstrap_v2_run.log", open = "wt")
sink(con, type = "output", split = FALSE)
sink(con, type = "message")

imp <- readRDS("imp_statin_htg_v4.rds")
frozen <- readRDS("frozen_cohort_ids.rds")
frozen_ids <- frozen$person_id

ps_vars <- c("age","male","white","hispanic","smoke_ever","alcohol",
             "ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
             "n_inpatient_ed","n_outpatient","n_visits_total",
             "htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
             "dementia","as_spond","ra","sle","pulmonary","cancer","neuro",
             "aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
             "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
             "antithromb","nsaid")

n_boot <- 200
n_months <- 12
total_reps <- 5 * n_boot

cat(sprintf("Bootstrap plan: 5 imp × %d boot = %d replicates\n", n_boot, total_reps))

# 预加载5个imputation数据
imp_data_list <- vector("list", 5)
for (i in 1:5) {
  d <- as.data.table(complete(imp, i))
  d <- d[person_id %in% frozen_ids, ]
  d[, arm := factor(arm, levels = c("NON_USER", "STATIN"))]
  for (v in c(ps_vars, "ap_event", "death_event", "time_days")) {
    if (v %in% names(d) && inherits(d[[v]], "integer64")) {
      d[, (v) := as.numeric(get(v))]
    }
  }
  imp_data_list[[i]] <- d
  cat(sprintf("Imputation %d loaded: N=%d (ST=%d, NU=%d)\n",
              i, nrow(d), sum(d$arm=="STATIN"), sum(d$arm=="NON_USER")))
}

# CIF计算函数 (完全复现原始v4的cif_curves)
calc_cif <- function(f_ap, f_d, n_months = 12) {
  out <- list()
  for (a in c("NON_USER", "STATIN")) {
    nd <- data.frame(arm = factor(a, levels = c("NON_USER", "STATIN")), month = 1:n_months)
    h_ap <- predict(f_ap, nd, type = "response")
    h_d  <- predict(f_d,  nd, type = "response")
    S <- 1
    cif_ap <- numeric(n_months)
    cif_d  <- numeric(n_months)
    surv_cs <- 1
    cif_ap_cs <- numeric(n_months)
    for (j in 1:n_months) {
      cif_ap[j] <- (if (j > 1) cif_ap[j-1] else 0) + h_ap[j] * S
      cif_d[j]  <- (if (j > 1) cif_d[j-1]  else 0) + h_d[j]  * S
      S <- S * (1 - h_ap[j] - h_d[j])
      surv_cs <- surv_cs * (1 - h_ap[j])
      cif_ap_cs[j] <- 1 - surv_cs
    }
    out[[a]] <- list(ap = cif_ap, d = cif_d, ap_cs = cif_ap_cs)
  }
  out
}

# 结果矩阵
n_cols <- 2*n_months + 2*n_months + 1 + 1 + 1 + n_months + n_months
col_names <- c(
  paste0("apNU_", 1:12), paste0("apST_", 1:12),
  paste0("dNU_", 1:12),  paste0("dST_", 1:12),
  "rr_ap", "rr_apcs", "rr_death",
  paste0("mrd_", 1:12), paste0("rd_", 1:12)
)
boot_results <- matrix(NA_real_, nrow = total_reps, ncol = length(col_names))
colnames(boot_results) <- col_names

ps_formula <- reformulate(ps_vars, response = "arm")

start_time <- Sys.time()
cat(sprintf("Start time: %s\n", start_time))

rep_idx <- 0
n_wt_fail <- 0L                       # PATCH T2c: weightit failure counter (discard-and-count)
smd_wt <- rep(NA_real_, total_reps)   # PATCH T2c: max |SMD| after weighting, per rep

# PATCH T2c: per-rep balance diagnostic (universal rule: every script producing
# an effect estimate must also report the max SMD after weighting for that run)
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

for (imp_i in 1:5) {
  d_imp <- imp_data_list[[imp_i]]

  for (b in 1:n_boot) {
    rep_idx <- rep_idx + 1

    # checkpoint: save partial results every 100 reps (does not alter analysis logic)
    if (rep_idx > 1 && rep_idx %% 100 == 0) {
      saveRDS(boot_results, "bootstrap_fresh_weightit_partial.rds")
      gc()   # explicit cleanup to prevent memory accumulation
    }

    # bootstrap 样本
    idx_b <- sample(nrow(d_imp), replace = TRUE)
    d_b <- d_imp[idx_b]

    # weightit ebal (ATT, focal = STATIN)
    # PATCH 2026-09-05 (T2c): failure = DISCARD replicate + count.
    # The original code silently fell back to w=1 (unweighted), mixing
    # unweighted replicates into the bootstrap distribution; audit of the
    # authoritative run: weightit_convergence_audit_main.csv
    wb <- tryCatch(
      weightit(ps_formula, data = d_b, method = "ebal",
               estimand = "ATT", focal = "STATIN"),
      error = function(e) e)
    if (inherits(wb, "error")) {
      n_wt_fail <- n_wt_fail + 1
      next   # row stays NA: replicate discarded and counted
    }
    d_b[, w := as.numeric(wb$weights)]
    smd_wt[rep_idx] <- max_smd_w(d_b, as.numeric(wb$weights))

    # person-month long table (1-based month, 1 to 12)
    d_b[, idx := .I]
    d_b[, t := pmin(time_days, 365)]
    d_b[, ev_ap := as.integer(ap_event == 1 & time_days <= 365)]
    d_b[, ev_d  := as.integer(death_event == 1 & time_days <= 365)]
    d_b[, M := pmin(pmax(1L, ceiling(t / (365/12))), 12L)]

    # uncount (用 rep 实现)
    n_per <- d_b$M
    n_total <- sum(n_per)
    idx_long <- rep(d_b$idx, n_per)
    pm <- d_b[idx_long, .(idx, arm, w, ev_ap, ev_d)]
    pm[, month := sequence(n_per)]
    pm[, last_m := month == rep(n_per, n_per)]
    pm[, ap_m := as.integer(ev_ap == 1 & last_m)]
    pm[, death_m := as.integer(ev_d == 1 & last_m)]

    # ---- 拟合模型 WITH interaction ----
    rr_ap <- NA_real_
    rr_apcs <- NA_real_
    rr_death <- NA_real_
    mrd_vec <- numeric(12)
    ap_cif_st <- ap_cif_nu <- d_cif_st <- d_cif_nu <- numeric(12)

    if (nrow(pm) > 100 && sum(pm$ap_m) > 3 && sum(pm$death_m) > 3) {
      tryCatch({
        f_ap <- glm(ap_m ~ arm * ns(month, df = 2, Boundary.knots = c(1, 12)),
                    family = binomial(), weights = w, data = pm)
        f_d  <- glm(death_m ~ arm * ns(month, df = 2, Boundary.knots = c(1, 12)),
                    family = binomial(), weights = w, data = pm)

        cc <- calc_cif(f_ap, f_d)
        ap_cif_nu <- cc$NON_USER$ap
        ap_cif_st <- cc$STATIN$ap
        d_cif_nu  <- cc$NON_USER$d
        d_cif_st  <- cc$STATIN$d

        rr_ap    <- ap_cif_st[12] / max(ap_cif_nu[12], 1e-10)
        rr_apcs  <- cc$STATIN$ap_cs[12] / max(cc$NON_USER$ap_cs[12], 1e-10)
        rr_death <- d_cif_st[12]  / max(d_cif_nu[12],  1e-10)

        rate_st <- diff(c(0, ap_cif_st))
        rate_nu <- diff(c(0, ap_cif_nu))
        mrd_vec <- (rate_st - rate_nu) * 1e5

      }, error = function(e) {})
    }

    # 存入结果
    boot_results[rep_idx, paste0("apNU_", 1:12)] <- ap_cif_nu * 1e5
    boot_results[rep_idx, paste0("apST_", 1:12)] <- ap_cif_st * 1e5
    boot_results[rep_idx, paste0("dNU_", 1:12)]  <- d_cif_nu * 1e5
    boot_results[rep_idx, paste0("dST_", 1:12)]  <- d_cif_st * 1e5
    boot_results[rep_idx, "rr_ap"]    <- rr_ap
    boot_results[rep_idx, "rr_apcs"]  <- rr_apcs
    boot_results[rep_idx, "rr_death"] <- rr_death
    boot_results[rep_idx, paste0("mrd_", 1:12)] <- mrd_vec
    boot_results[rep_idx, paste0("rd_", 1:12)]  <- (ap_cif_st - ap_cif_nu) * 1e5

    # 进度
    if (b %% 50 == 0) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      eta <- elapsed / rep_idx * (total_reps - rep_idx)
      cat(sprintf("  Imp %d / Rep %d (%.1f%%) | elapsed=%.1f min | ETA=%.1f min | RR_ap=%.3f\n",
                  imp_i, b, 100*rep_idx/total_reps, elapsed, eta,
                  rr_ap))
    }
  }
}

end_time <- Sys.time()
cat(sprintf("\nDone! Total time: %.1f minutes\n", difftime(end_time, start_time, units = "mins")))
cat(sprintf("weightit failures (discarded): %d / %d\n", n_wt_fail, total_reps))   # PATCH T2c
cat(sprintf("max|SMD| after weighting: median=%.2e max=%.2e | any >0.01: %d\n",   # PATCH T2c
            median(smd_wt, na.rm = TRUE), max(smd_wt, na.rm = TRUE),
            sum(smd_wt > 0.01, na.rm = TRUE)))

# ---- 汇总 ----
rr_ap_valid <- boot_results[, "rr_ap"]
rr_ap_valid <- rr_ap_valid[is.finite(rr_ap_valid) & rr_ap_valid > 0 & rr_ap_valid < 100]
cat(sprintf("Valid RR_ap: %d / %d\n", length(rr_ap_valid), total_reps))
cat(sprintf("RR_ap median: %.3f\n", median(rr_ap_valid)))
cat(sprintf("RR_ap SD: %.3f\n", sd(rr_ap_valid)))
cat(sprintf("RR_ap 95%% CI: %.3f - %.3f\n",
            quantile(rr_ap_valid, 0.025), quantile(rr_ap_valid, 0.975)))

rr_d_valid <- boot_results[, "rr_death"]
rr_d_valid <- rr_d_valid[is.finite(rr_d_valid) & rr_d_valid > 0 & rr_d_valid < 100]
cat(sprintf("RR_death median: %.3f\n", median(rr_d_valid)))
cat(sprintf("RR_death 95%% CI: %.3f - %.3f\n",
            quantile(rr_d_valid, 0.025), quantile(rr_d_valid, 0.975)))

# 保存
saveRDS(boot_results, "bootstrap_fresh_weightit.rds")
cat("Saved to bootstrap_fresh_weightit.rds\n")

# close sink
sink(type = "message")
sink(type = "output")
close(con)

# remove partial checkpoint file on success
if (file.exists("bootstrap_fresh_weightit_partial.rds")) {
  file.remove("bootstrap_fresh_weightit_partial.rds")
}