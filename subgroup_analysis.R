# ============================================================================
# Subgroup analysis: Age (<65 vs >=65) and Sex (Male vs Female)
# Method: within-subgroup entropy balancing (WeightIt/ebal, ATT, focal=STATIN)
#         + weighted discrete-time survival (pooled logistic, arm × ns(month, df=2))
#         + competing-risk CIF + RR
# Bootstrap: 5 MICE × 200 = 1000 reps per subgroup
# Seed: 20260625 (same as main analysis)
# Stratifying variable removed from ps_vars within each subgroup
# ============================================================================
library(mice); library(WeightIt); library(dplyr); library(data.table); library(splines)

set.seed(20260625)   # global seed

# sink log to file
con <- file("subgroup_analysis_run.log", open = "wt")
sink(con, type = "output", split = FALSE)
sink(con, type = "message")

# ============================================================================
# 1. Load data
# ============================================================================
imp <- readRDS("imp_statin_htg_v4.rds")
frozen <- readRDS("frozen_cohort_ids.rds")
frozen_ids <- frozen$person_id

ps_vars_full <- c("age","male","white","hispanic","smoke_ever","alcohol",
                  "ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
                  "n_inpatient_ed","n_outpatient","n_visits_total",
                  "htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
                  "dementia","as_spond","ra","sle","pulmonary","cancer","neuro",
                  "aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
                  "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
                  "antithromb","nsaid")

n_boot <- 200
n_months <- 12

# Preload 5 imputations (frozen cohort only)
imp_data_list <- vector("list", 5)
for (i in 1:5) {
  d <- as.data.table(complete(imp, i))
  d <- d[person_id %in% frozen_ids, ]
  d[, arm := factor(arm, levels = c("NON_USER", "STATIN"))]
  for (v in c(ps_vars_full, "ap_event", "death_event", "time_days")) {
    if (v %in% names(d) && inherits(d[[v]], "integer64")) {
      d[, (v) := as.numeric(get(v))]
    }
  }
  # create subgroup indicators
  d[, age_grp := ifelse(age < 65, "<65", ">=65")]
  imp_data_list[[i]] <- d
  cat(sprintf("Imputation %d loaded: N=%d\n", i, nrow(d)))
}

# ============================================================================
# 2. CIF calculation function (same as main analysis)
# ============================================================================
calc_cif <- function(f_ap, f_d, n_months = 12) {
  out <- list()
  for (a in c("NON_USER", "STATIN")) {
    nd <- data.frame(arm = factor(a, levels = c("NON_USER", "STATIN")), month = 1:n_months)
    h_ap <- predict(f_ap, nd, type = "response")
    h_d  <- predict(f_d,  nd, type = "response")
    S <- 1
    cif_ap <- numeric(n_months)
    cif_d  <- numeric(n_months)
    for (j in 1:n_months) {
      cif_ap[j] <- (if (j > 1) cif_ap[j-1] else 0) + h_ap[j] * S
      cif_d[j]  <- (if (j > 1) cif_d[j-1]  else 0) + h_d[j]  * S
      S <- S * (1 - h_ap[j] - h_d[j])
    }
    out[[a]] <- list(ap = cif_ap, d = cif_d)
  }
  out
}

# ============================================================================
# 3. Bootstrap function for one subgroup
# ============================================================================
run_subgroup_bootstrap <- function(imp_data_list, subgroup_filter, ps_vars_sub,
                                    subgroup_label, n_boot = 200) {
  cat(sprintf("\n{'label': '%s', 'status': 'starting'}\n", subgroup_label))
  total_reps <- 5 * n_boot

  # check sample size and events per subgroup
  d_check <- imp_data_list[[1]][eval(subgroup_filter)]
  n_total <- nrow(d_check)
  n_st <- sum(d_check$arm == "STATIN")
  n_nu <- sum(d_check$arm == "NON_USER")
  n_ap <- sum(d_check$ap_event)
  n_death <- sum(d_check$death_event)
  epv_covar <- n_ap / length(ps_vars_sub)
  epv_model <- n_ap / 5   # arm + ns(month,df=2) + interaction ≈ 5 params
  cat(sprintf("  N=%d (ST=%d, NU=%d)  AP=%d  Death=%d  EPV(covar=%.1f, model=%.1f)\n",
              n_total, n_st, n_nu, n_ap, n_death, epv_covar, epv_model))
  if (epv_covar < 10) cat(sprintf("  WARNING: EPV(covar) < 10, CI may be wide\n"))

  col_names <- c(
    paste0("apNU_", 1:12), paste0("apST_", 1:12),
    paste0("dNU_", 1:12),  paste0("dST_", 1:12),
    "rr_ap", "rr_death",
    paste0("mrd_", 1:12)
  )
  boot_results <- matrix(NA_real_, nrow = total_reps, ncol = length(col_names))
  colnames(boot_results) <- col_names

  ps_formula <- reformulate(ps_vars_sub, response = "arm")
  start_time <- Sys.time()
  rep_idx <- 0
  fail_count <- 0

  for (imp_i in 1:5) {
    d_imp <- imp_data_list[[imp_i]][eval(subgroup_filter)]

    for (b in 1:n_boot) {
      rep_idx <- rep_idx + 1

      # checkpoint every 100 reps
      if (rep_idx > 1 && rep_idx %% 100 == 0) {
        saveRDS(boot_results,
                sprintf("subgroup_%s_partial.rds", gsub("[^a-zA-Z0-9]", "_", subgroup_label)))
        gc()
      }

      # bootstrap sample
      idx_b <- sample(nrow(d_imp), replace = TRUE)
      d_b <- d_imp[idx_b]

      # weightit ebal (ATT, focal = STATIN)
      tryCatch({
        w_b <- weightit(ps_formula, data = d_b, method = "ebal",
                         estimand = "ATT", focal = "STATIN")
        d_b[, w := as.numeric(w_b$weights)]
      }, error = function(e) {
        d_b[, w := 1.0]
      })

      # person-month long table
      d_b[, idx := .I]
      d_b[, t := pmin(time_days, 365)]
      d_b[, ev_ap := as.integer(ap_event == 1 & time_days <= 365)]
      d_b[, ev_d  := as.integer(death_event == 1 & time_days <= 365)]
      d_b[, M := pmin(pmax(1L, ceiling(t / (365/12))), 12L)]

      n_per <- d_b$M
      idx_long <- rep(d_b$idx, n_per)
      pm <- d_b[idx_long, .(idx, arm, w, ev_ap, ev_d)]
      pm[, month := sequence(n_per)]
      pm[, last_m := month == rep(n_per, n_per)]
      pm[, ap_m := as.integer(ev_ap == 1 & last_m)]
      pm[, death_m := as.integer(ev_d == 1 & last_m)]

      rr_ap <- NA_real_
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
          rr_death <- d_cif_st[12]  / max(d_cif_nu[12],  1e-10)

          rate_st <- diff(c(0, ap_cif_st))
          rate_nu <- diff(c(0, ap_cif_nu))
          mrd_vec <- (rate_st - rate_nu) * 1e5

        }, error = function(e) {})
      } else {
        fail_count <- fail_count + 1
      }

      boot_results[rep_idx, paste0("apNU_", 1:12)] <- ap_cif_nu * 1e5
      boot_results[rep_idx, paste0("apST_", 1:12)] <- ap_cif_st * 1e5
      boot_results[rep_idx, paste0("dNU_", 1:12)]  <- d_cif_nu * 1e5
      boot_results[rep_idx, paste0("dST_", 1:12)]  <- d_cif_st * 1e5
      boot_results[rep_idx, "rr_ap"]    <- rr_ap
      boot_results[rep_idx, "rr_death"] <- rr_death
      boot_results[rep_idx, paste0("mrd_", 1:12)] <- mrd_vec

      if (b %% 50 == 0) {
        elapsed <- difftime(Sys.time(), start_time, units = "mins")
        eta <- elapsed / rep_idx * (total_reps - rep_idx)
        cat(sprintf("  %s | Imp %d / Rep %d (%.0f%%) | elapsed=%.1f min | ETA=%.1f min | RR_ap=%.3f\n",
                    subgroup_label, imp_i, b, 100*rep_idx/total_reps, elapsed, eta, rr_ap))
      }
    }
  }

  end_time <- Sys.time()
  cat(sprintf("  %s DONE: %.1f minutes | fails=%d\n",
              subgroup_label, difftime(end_time, start_time, units = "mins"), fail_count))

  # summarize
  rr_ap_valid <- boot_results[, "rr_ap"]
  rr_ap_valid <- rr_ap_valid[is.finite(rr_ap_valid) & rr_ap_valid > 0 & rr_ap_valid < 100]
  rr_d_valid <- boot_results[, "rr_death"]
  rr_d_valid <- rr_d_valid[is.finite(rr_d_valid) & rr_d_valid > 0 & rr_d_valid < 100]

  cat(sprintf("  %s Summary: AP RR=%.3f (%.3f-%.3f) | Death RR=%.3f (%.3f-%.3f) | valid=%d/%d\n",
              subgroup_label,
              median(rr_ap_valid), quantile(rr_ap_valid, 0.025), quantile(rr_ap_valid, 0.975),
              median(rr_d_valid), quantile(rr_d_valid, 0.025), quantile(rr_d_valid, 0.975),
              length(rr_ap_valid), total_reps))

  boot_results
}

# ============================================================================
# 4. Run 4 subgroups
# ============================================================================
# Age <65: remove "age" from ps_vars
ps_vars_age <- setdiff(ps_vars_full, "age")
# Sex: remove "male" from ps_vars
ps_vars_sex <- setdiff(ps_vars_full, "male")

subgroups <- list(
  list(label = "Age<65",     filter = quote(age_grp == "<65"),  ps_vars = ps_vars_age),
  list(label = "Age>=65",    filter = quote(age_grp == ">=65"), ps_vars = ps_vars_age),
  list(label = "Male",       filter = quote(male == 1),         ps_vars = ps_vars_sex),
  list(label = "Female",     filter = quote(male == 0),         ps_vars = ps_vars_sex)
)

all_results <- list()
for (sg in subgroups) {
  res <- run_subgroup_bootstrap(imp_data_list, sg$filter, sg$ps_vars, sg$label, n_boot)
  all_results[[sg$label]] <- res
  saveRDS(all_results, "subgroup_all_results.rds")
  cat(sprintf("Saved partial results (completed: %s)\n", paste(names(all_results), collapse=", ")))
}

# ============================================================================
# 5. Interaction P-values (full cohort, main-analysis weights)
# ============================================================================
cat("\n=== Interaction P-values ===\n")

# Use imp[1] with main-analysis weights
d_full <- imp_data_list[[1]]
ps_formula_full <- reformulate(ps_vars_full, response = "arm")
w_full <- weightit(ps_formula_full, data = d_full, method = "ebal",
                   estimand = "ATT", focal = "STATIN")
d_full[, w := as.numeric(w_full$weights)]

# person-month long table
d_full[, idx := .I]
d_full[, t := pmin(time_days, 365)]
d_full[, ev_ap := as.integer(ap_event == 1 & time_days <= 365)]
d_full[, ev_d  := as.integer(death_event == 1 & time_days <= 365)]
d_full[, M := pmin(pmax(1L, ceiling(t / (365/12))), 12L)]

n_per <- d_full$M
idx_long <- rep(d_full$idx, n_per)
pm_full <- d_full[idx_long, .(idx, arm, w, ev_ap, ev_d, age_grp, male)]
pm_full[, month := sequence(n_per)]
pm_full[, last_m := month == rep(n_per, n_per)]
pm_full[, ap_m := as.integer(ev_ap == 1 & last_m)]
pm_full[, death_m := as.integer(ev_d == 1 & last_m)]

# Age interaction: arm × age_grp
pm_full[, age_grp := factor(age_grp, levels = c("<65", ">=65"))]
f_age_int <- glm(ap_m ~ arm * age_grp + ns(month, df=2, Boundary.knots=c(1,12)),
                 family = binomial(), weights = w, data = pm_full)
# extract armSTATIN:age_grp>=65 interaction term
age_int_term <- coef(summary(f_age_int))
age_int_row <- grep("armSTATIN:age_grp>=65", rownames(age_int_term), value = TRUE)
if (length(age_int_row) > 0) {
  z <- age_int_term[age_int_row, "z value"]
  p_age <- 2 * pnorm(-abs(z))
} else {
  p_age <- NA
}
cat(sprintf("  Age interaction P = %.4f\n", p_age))

# Sex interaction: arm × male
f_sex_int <- glm(ap_m ~ arm * factor(male) + ns(month, df=2, Boundary.knots=c(1,12)),
                 family = binomial(), weights = w, data = pm_full)
sex_int_term <- coef(summary(f_sex_int))
sex_int_row <- grep("armSTATIN:factor\\(male\\)1", rownames(sex_int_term), value = TRUE)
if (length(sex_int_row) > 0) {
  z <- sex_int_term[sex_int_row, "z value"]
  p_sex <- 2 * pnorm(-abs(z))
} else {
  p_sex <- NA
}
cat(sprintf("  Sex interaction P = %.4f\n", p_sex))

# ============================================================================
# 6. Final summary table
# ============================================================================
cat("\n=== FINAL SUBGROUP RESULTS ===\n")
cat(sprintf("%-12s | %-8s | %-25s | %-25s | %s\n",
            "Subgroup", "Outcome", "RR (95% CI)", "Events ST/NU", "EPV"))

summary_rows <- list()
for (sg_label in names(all_results)) {
  boot_mat <- all_results[[sg_label]]
  rr_ap_v <- boot_mat[, "rr_ap"]
  rr_ap_v <- rr_ap_v[is.finite(rr_ap_v) & rr_ap_v > 0 & rr_ap_v < 100]
  rr_d_v <- boot_mat[, "rr_death"]
  rr_d_v <- rr_d_v[is.finite(rr_d_v) & rr_d_v > 0 & rr_d_v < 100]

  # get event counts from imp[1]
  sg <- subgroups[[which(sapply(subgroups, function(x) x$label) == sg_label)]]
  d_sg <- imp_data_list[[1]][eval(sg$filter)]
  n_ap_st <- sum(d_sg$ap_event[d_sg$arm == "STATIN"])
  n_ap_nu <- sum(d_sg$ap_event[d_sg$arm == "NON_USER"])
  n_d_st <- sum(d_sg$death_event[d_sg$arm == "STATIN"])
  n_d_nu <- sum(d_sg$death_event[d_sg$arm == "NON_USER"])
  epv <- (n_ap_st + n_ap_nu) / length(sg$ps_vars)

  # AP
  ap_rr <- median(rr_ap_v)
  ap_lo <- quantile(rr_ap_v, 0.025)
  ap_hi <- quantile(rr_ap_v, 0.975)
  cat(sprintf("%-12s | %-8s | %.3f (%.3f-%.3f)        | %d/%d            | %.1f\n",
              sg_label, "AP", ap_rr, ap_lo, ap_hi, n_ap_st, n_ap_nu, epv))

  # Death
  d_rr <- median(rr_d_v)
  d_lo <- quantile(rr_d_v, 0.025)
  d_hi <- quantile(rr_d_v, 0.975)
  cat(sprintf("%-12s | %-8s | %.3f (%.3f-%.3f)        | %d/%d            | \n",
              sg_label, "Death", d_rr, d_lo, d_hi, n_d_st, n_d_nu))

  summary_rows[[paste0(sg_label, "_AP")]] <- data.frame(
    subgroup = sg_label, outcome = "AP",
    rr = ap_rr, ci_lo = ap_lo, ci_hi = ap_hi,
    events_statin = n_ap_st, events_nonuser = n_ap_nu,
    epv = epv, stringsAsFactors = FALSE
  )
  summary_rows[[paste0(sg_label, "_Death")]] <- data.frame(
    subgroup = sg_label, outcome = "Death",
    rr = d_rr, ci_lo = d_lo, ci_hi = d_hi,
    events_statin = n_d_st, events_nonuser = n_d_nu,
    epv = epv, stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
summary_df$interaction_p <- NA
summary_df[summary_df$outcome == "AP" & summary_df$subgroup %in% c("Age<65","Age>=65"), "interaction_p"] <- p_age
summary_df[summary_df$outcome == "AP" & summary_df$subgroup %in% c("Male","Female"), "interaction_p"] <- p_sex

write.csv(summary_df, "subgroup_results.csv", row.names = FALSE)
cat("\nSaved: subgroup_results.csv\n")
cat(sprintf("Interaction P-values: Age=%.4f, Sex=%.4f\n", p_age, p_sex))

# save all bootstrap results
saveRDS(all_results, "subgroup_all_results.rds")
cat("Saved: subgroup_all_results.rds\n")

# cleanup partial files
for (sg in subgroups) {
  f <- sprintf("subgroup_%s_partial.rds", gsub("[^a-zA-Z0-9]", "_", sg$label))
  if (file.exists(f)) file.remove(f)
}

# close sink
sink(type = "message")
sink(type = "output")
close(con)
cat("Done.\n")
