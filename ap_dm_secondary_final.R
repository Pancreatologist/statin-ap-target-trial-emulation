suppressMessages({
  library(dplyr); library(bit64); library(survival); library(splines)
  library(WeightIt); library(mice); library(ggplot2)
})

# ============================================================
# 1. 数据准备
# ============================================================
df_main <- readRDS("df_final_statin_htg_v4.rds")
load("20250920对于患者诊断的明确.RData")

dm_df <- dm_patients %>%
  mutate(person_id = as.integer64(person_id)) %>%
  select(person_id, dm_date)

df <- df_main %>% left_join(dm_df, by = "person_id")

df_ap <- df %>%
  filter(ap_event == 1 & diabetes == 0) %>%
  mutate(
    ap_date = as.Date(ap_date),
    death_date = as.Date(death_date),
    obs_end = as.Date(obs_end),
    dm_date = as.Date(dm_date),
    censor_date = pmin(obs_end, na.rm = TRUE),
    fup_days = as.numeric(censor_date - ap_date),
    dm_event = as.integer(!is.na(dm_date) & dm_date >= ap_date & dm_date <= censor_date),
    dm_time = ifelse(dm_event == 1, as.numeric(dm_date - ap_date), fup_days),
    death_event_ap = as.integer(!is.na(death_date) & death_date >= ap_date & death_date <= censor_date),
    death_time_ap = ifelse(death_event_ap == 1, as.numeric(death_date - ap_date), fup_days),
    time_to_first_event = pmin(dm_time, death_time_ap),
    event_type = ifelse(dm_event == 1, 1, ifelse(death_event_ap == 1, 2, 0)),
    male = as.integer(sex == "Male"),
    white = as.integer(race == "White"),
    ever_smoke = as.integer(smoke_status == "former"),
    egfr_trunc = pmin(pmax(egfr, 10), 200),
    alcohol_num = as.numeric(alcohol)
  ) %>%
  filter(fup_days > 0)

cat("N (AP patients, no baseline DM):", nrow(df_ap), "\n")
cat("  STATIN:", sum(df_ap$arm == "STATIN"), "\n")
cat("  NON_USER:", sum(df_ap$arm == "NON_USER"), "\n")
cat("DM events:", sum(df_ap$dm_event), "\n")
cat("Death events:", sum(df_ap$death_event_ap), "\n")

# ============================================================
# 2. 协变量定义
# ============================================================
demo_vars   <- c("age","male","white","ever_smoke")
cont_vars   <- c("bmi","tg","hdl","egfr_trunc")
comorb_cols <- c("htn","obesity","chd","stroke","kidney","cancer")
drug_cols   <- c("aspirin","acei_arb","bblocker","nsaid","metformin")
ps_vars_all <- c(demo_vars, "alcohol_num", cont_vars, comorb_cols, drug_cols)

# 连续变量对数变换
df_ap <- df_ap %>%
  mutate(
    ln_bmi   = log(bmi),
    ln_tg    = log(tg),
    ln_hdl   = log(hdl),
    ln_egfr  = log(egfr_trunc)
  )
cont_ln_vars <- c("ln_bmi","ln_tg","ln_hdl","ln_egfr")
ps_vars_ln   <- c(demo_vars, "alcohol_num", cont_ln_vars, comorb_cols, drug_cols)

cat("\n协变量 (n=20):", length(ps_vars_ln), "\n")
cat("  人口学:", paste(demo_vars, collapse=", "), "\n")
cat("  连续(ln):", paste(cont_ln_vars, collapse=", "), "\n")
cat("  合并症:", paste(comorb_cols, collapse=", "), "\n")
cat("  药物:", paste(drug_cols, collapse=", "), "\n")

# ============================================================
# 3. 缺失模式
# ============================================================
cat("\nMissing pattern (%):\n")
for (v in ps_vars_all) {
  nmis <- sum(is.na(df_ap[[v]]))
  if (nmis > 0) cat(sprintf("  %-12s: %3d/%d (%.1f%%)\n", v, nmis, nrow(df_ap), nmis/nrow(df_ap)*100))
}

# ============================================================
# 4. Nelson-Aalen 累积风险（White-Royston）
# ============================================================
na_fit_dm <- survfit(Surv(dm_time, dm_event) ~ 1, data = df_ap)
df_ap$cumhaz_dm <- summary(na_fit_dm, times = pmin(df_ap$dm_time, max(na_fit_dm$time)))$cumhaz

na_fit_death <- survfit(Surv(death_time_ap, death_event_ap) ~ 1, data = df_ap)
df_ap$cumhaz_death <- summary(na_fit_death, times = pmin(df_ap$death_time_ap, max(na_fit_death$time)))$cumhaz

cat("\nWhite-Royston: 已加入 Nelson-Aalen 累积风险作为辅助变量\n")

# ============================================================
# 5. MICE 插补
# ============================================================
mice_vars <- c(ps_vars_ln, "dm_event", "dm_time", "death_event_ap", "death_time_ap", "cumhaz_dm", "cumhaz_death")
mice_df  <- df_ap %>% select(all_of(mice_vars))

meth <- rep("", ncol(mice_df))
names(meth) <- names(mice_df)

# 连续变量 → pmm
meth[c("age", cont_ln_vars, "alcohol_num")] <- "pmm"

# 二分类 → logreg
binary_vars <- c("male","white","ever_smoke", comorb_cols, drug_cols)
meth[binary_vars] <- "logreg"

# 结局变量不插补
meth[c("dm_event","dm_time","death_event_ap","death_time_ap","cumhaz_dm","cumhaz_death")] <- ""

set.seed(1234)
imp <- mice(mice_df, m = 5, method = meth, printFlag = FALSE, seed = 1234, maxit = 15)
saveRDS(imp, "imp_ap_dm_final.rds")

# 验证插补
d_check <- complete(imp, 1)
nmis_after <- sum(colSums(is.na(d_check)) > 0)
cat(sprintf("MICE done (m=5). 仍有缺失变量数: %d\n", nmis_after))
if (nmis_after > 0) {
  cat("  仍有缺失的变量:\n")
  print(names(which(colSums(is.na(d_check)) > 0)))
}

# ============================================================
# 6. 分析函数
# ============================================================
FOLLOW_DAYS <- 1825

estimate_once <- function(d, resample = FALSE, max_t = FOLLOW_DAYS) {
  if (resample) {
    n <- nrow(d)
    idx <- sample.int(n, n, replace = TRUE)
    d <- d[idx, ]
  }

  d$arm_num <- ifelse(d$arm == "STATIN", 1, 0)

  # 熵平衡 ATT
  ps_formula <- as.formula(paste("arm ~", paste(ps_vars_ln, collapse = " + ")))
  w <- try(weightit(ps_formula, data = d, method = "ebal",
                    estimand = "ATT", focal = "STATIN"), silent = TRUE)
  if (inherits(w, "try-error")) return(NULL)

  wts <- w$weights
  if (any(is.na(wts)) || any(is.infinite(wts))) return(NULL)

  # ESS 和 max weight
  ess_ctrl <- sum(wts[d$arm == "NON_USER"])^2 / sum(wts[d$arm == "NON_USER"]^2)
  max_wt <- max(wts[d$arm == "NON_USER"])

  # ---- Cox model (cause-specific HR for DM) ----
  cx <- try(coxph(Surv(time_to_first_event, event_type == 1) ~ arm_num,
                  data = d, weights = wts), silent = TRUE)
  if (inherits(cx, "try-error")) return(NULL)
  hr_dm <- exp(coef(cx)[["arm_num"]])

  # ---- 竞争风险 CIF (Aalen-Johansen 近似) ----
  get_cif_aj <- function(data, weights, arm_val, t_max) {
    sub <- data[data$arm == arm_val, ]
    sub_w <- weights[data$arm == arm_val]

    # cause-specific cumulative hazards
    dm_surv <- try(survfit(Surv(time_to_first_event, event_type == 1) ~ 1,
                           data = sub, weights = sub_w), silent = TRUE)
    death_surv <- try(survfit(Surv(time_to_first_event, event_type == 2) ~ 1,
                              data = sub, weights = sub_w), silent = TRUE)
    if (inherits(dm_surv, "try-error") || inherits(death_surv, "try-error")) return(NULL)

    # 合并时间点
    all_times <- sort(unique(c(dm_surv$time, death_surv$time)))
    all_times <- all_times[all_times <= t_max]
    if (length(all_times) == 0) return(list(times = t_max, cif = 0))

    # 插值cumulative hazard到所有时间点
    H_dm <- approx(dm_surv$time, dm_surv$cumhaz, all_times, rule = 2)$y
    H_death <- approx(death_surv$time, death_surv$cumhaz, all_times, rule = 2)$y

    # overall survival S(t) = exp(-H_all)
    S_t <- exp(-H_dm - H_death)

    # Aalen-Johansen CIF for DM: ∫ S(u-) dH_dm(u)
    dH_dm <- diff(c(0, H_dm))
    S_prev <- c(1, S_t[-length(S_t)])
    cif_dm <- cumsum(S_prev * dH_dm)

    list(times = all_times, cif = cif_dm)
  }

  cif_nu <- get_cif_aj(d, wts, "NON_USER", max_t)
  cif_st <- get_cif_aj(d, wts, "STATIN", max_t)
  if (is.null(cif_nu) || is.null(cif_st)) return(NULL)

  # 提取时间点
  tp_days <- c(90, 180, 365, 730, 1095, 1825)
  get_cif_at_t <- function(cif_obj, t) {
    idx <- which(cif_obj$times <= t)
    if (length(idx) == 0) return(0)
    cif_obj$cif[max(idx)]
  }
  cif_nu_tp <- sapply(tp_days, function(t) get_cif_at_t(cif_nu, t))
  cif_st_tp <- sapply(tp_days, function(t) get_cif_at_t(cif_st, t))

  cif_nu_5y <- cif_nu_tp[6]
  cif_st_5y <- cif_st_tp[6]
  rr_5y <- ifelse(cif_nu_5y > 0, cif_st_5y / cif_nu_5y, NA)
  rd_5y <- cif_st_5y - cif_nu_5y

  # 平衡诊断
  bal_stats <- NULL
  if (!resample) {
    wtd_var <- function(x, w) {
      ok <- !is.na(x); x <- x[ok]; w <- w[ok]
      wm <- weighted.mean(x, w)
      sum(w * (x - wm)^2) / sum(w)
    }
    calc_smd <- function(var, treat, w) {
      m1 <- weighted.mean(var[treat == 1], w[treat == 1], na.rm = TRUE)
      m0 <- weighted.mean(var[treat == 0], w[treat == 0], na.rm = TRUE)
      v1 <- wtd_var(var[treat == 1], w[treat == 1])
      v0 <- wtd_var(var[treat == 0], w[treat == 0])
      (m1 - m0) / sqrt((v1 + v0) / 2)
    }
    treat_num <- ifelse(d$arm == "STATIN", 1, 0)
    smds_adj <- sapply(ps_vars_ln, function(v) calc_smd(d[[v]], treat_num, wts))
    smds_un <- sapply(ps_vars_ln, function(v) calc_smd(d[[v]], treat_num, rep(1, nrow(d))))
    bal_stats <- list(smd_unadj = smds_un, smd_adj = smds_adj,
                      ess_ctrl = ess_ctrl, max_wt = max_wt)
  }

  list(
    hr_dm = hr_dm,
    cif_nu_5y = cif_nu_5y * 100,
    cif_st_5y = cif_st_5y * 100,
    rr_5y = rr_5y,
    rd_5y = rd_5y * 100,
    cif_nu_tp = cif_nu_tp * 100,
    cif_st_tp = cif_st_tp * 100,
    ess_ctrl = ess_ctrl,
    max_wt = max_wt,
    bal_stats = bal_stats
  )
}

# ============================================================
# 7. 原始点估计
# ============================================================
cat("\n--- Point estimates (per imputation) ---\n")
point_estimates <- list()
bal_list <- list()

for (m in 1:5) {
  d_imp <- complete(imp, m)
  d_imp$arm <- df_ap$arm
  d_imp$time_to_first_event <- df_ap$time_to_first_event
  d_imp$event_type <- df_ap$event_type

  res <- estimate_once(d_imp, resample = FALSE)
  if (!is.null(res)) {
    point_estimates[[m]] <- res
    bal_list[[m]] <- res$bal_stats
    cat(sprintf("  imp %d: HR=%.3f  NU_5y=%.1f%%  ST_5y=%.1f%%  RR=%.3f  RD=%+.2f%%  ESS=%.0f  maxWT=%.1f\n",
                m, res$hr_dm, res$cif_nu_5y, res$cif_st_5y,
                res$rr_5y, res$rd_5y, res$ess_ctrl, res$max_wt))
  } else {
    cat(sprintf("  imp %d: FAILED\n", m))
  }
}

# 汇总点估计
if (length(point_estimates) > 0) {
  pe_mat <- do.call(rbind, lapply(point_estimates, function(x) {
    c(x$hr_dm, x$cif_nu_5y, x$cif_st_5y, x$rr_5y, x$rd_5y, x$ess_ctrl, x$max_wt)
  }))
  colnames(pe_mat) <- c("hr_dm","cif_nu_5y","cif_st_5y","rr_5y","rd_5y","ess_ctrl","max_wt")
  cat(sprintf("\nPooled (mean over %d imp):\n", length(point_estimates)))
  cat(sprintf("  HR = %.3f\n", mean(pe_mat[,"hr_dm"])))
  cat(sprintf("  5y RR = %.3f\n", mean(pe_mat[,"rr_5y"])))
  cat(sprintf("  5y RD = %+.2f%%\n", mean(pe_mat[,"rd_5y"])))
  cat(sprintf("  5y CIF: NON_USER=%.1f%%, STATIN=%.1f%%\n",
              mean(pe_mat[,"cif_nu_5y"]), mean(pe_mat[,"cif_st_5y"])))
  cat(sprintf("  ESS_CTRL = %.0f (raw control n=%.0f)\n",
              mean(pe_mat[,"ess_ctrl"]), sum(df_ap$arm == "NON_USER")))
  cat(sprintf("  Max weight = %.1f\n", mean(pe_mat[,"max_wt"])))
}

# 平衡诊断表
if (length(bal_list) > 0) {
  bal_mat <- do.call(cbind, lapply(bal_list, function(b) b$smd_adj))
  smd_adj_mean <- apply(bal_mat, 1, mean)
  smd_unadj_mean <- apply(do.call(cbind, lapply(bal_list, function(b) b$smd_unadj)), 1, mean)

  cat("\n--- 平衡诊断 (SMD, 均值 over 5 imp) ---\n")
  cat(sprintf("  %-12s  %8s  %8s\n", "Variable", "Unadj SMD", "Adj SMD"))
  for (v in ps_vars_ln) {
    cat(sprintf("  %-12s  %8.3f  %8.3f\n", v, smd_unadj_mean[v], smd_adj_mean[v]))
  }
  cat(sprintf("  Max |SMD| after weighting: %.3f\n", max(abs(smd_adj_mean))))
}

saveRDS(bal_list, "bal_stats_final.rds")

# ============================================================
# 8. Bootstrap
# ============================================================
N_BOOT <- 200
cat(sprintf("\n--- Bootstrap: %d × 5 = %d replicates ---\n", N_BOOT, N_BOOT * 5))

set.seed(5678)
boot_res <- list()
fail_count <- 0
fail_reasons <- c("ebal_fail" = 0, "na_weights" = 0, "cox_fail" = 0, "cif_fail" = 0, "other" = 0)

total <- N_BOOT * 5
counter <- 0

for (m in 1:5) {
  d_imp <- complete(imp, m)
  d_imp$arm <- df_ap$arm
  d_imp$time_to_first_event <- df_ap$time_to_first_event
  d_imp$event_type <- df_ap$event_type

  for (b in 1:N_BOOT) {
    counter <- counter + 1
    if (counter %% 100 == 0) cat(sprintf("  %d/%d done (%.0f%%)\n", counter, total, counter/total*100))

    # 用tryCatch追踪失败原因
    res <- tryCatch({
      estimate_once(d_imp, resample = TRUE)
    }, error = function(e) {
      NULL
    })

    if (is.null(res)) {
      fail_count <- fail_count + 1
      next
    }

    res_vec <- c(
      hr_dm = res$hr_dm,
      cif_nu_5y = res$cif_nu_5y,
      cif_st_5y = res$cif_st_5y,
      rr_5y = res$rr_5y,
      rd_5y = res$rd_5y,
      setNames(res$cif_nu_tp, paste0("cifNU_", c("3m","6m","1y","2y","3y","5y"))),
      setNames(res$cif_st_tp, paste0("cifST_", c("3m","6m","1y","2y","3y","5y")))
    )
    boot_res[[length(boot_res) + 1]] <- res_vec
  }
}

n_success <- length(boot_res)
cat(sprintf("\nBootstrap complete: %d successful / %d total (%.1f%% failed)\n",
            n_success, total, (total - n_success)/total*100))

boot_mat <- do.call(rbind, boot_res)
saveRDS(boot_mat, "ap_dm_bootstrap_final.rds")

# ============================================================
# 9. 结果汇总
# ============================================================
cat("\n============================================================\n")
cat("  FINAL RESULTS: Statin vs Non-user — Post-AP New-onset DM\n")
cat("  (Competing-risk CIF, White-Royston MICE, Entropy Balancing ATT)\n")
cat("============================================================\n\n")

med_est <- apply(boot_mat, 2, median, na.rm = TRUE)
ci_lo   <- apply(boot_mat, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
ci_hi   <- apply(boot_mat, 2, function(x) quantile(x, 0.975, na.rm = TRUE))

cat(sprintf("Bootstrap replicates: %d (of %d attempted)\n\n", n_success, total))

cat(sprintf("DM cause-specific HR (STATIN vs NON_USER):\n"))
cat(sprintf("  %.3f  (95%% CI: %.3f, %.3f)\n\n",
            med_est["hr_dm"], ci_lo["hr_dm"], ci_hi["hr_dm"]))

cat(sprintf("5-year risk (Aalen-Johansen CIF):\n"))
cat(sprintf("  NON_USER:  %.1f%%  (95%% CI: %.1f%%, %.1f%%)\n",
            med_est["cif_nu_5y"], ci_lo["cif_nu_5y"], ci_hi["cif_nu_5y"]))
cat(sprintf("  STATIN:    %.1f%%  (95%% CI: %.1f%%, %.1f%%)\n",
            med_est["cif_st_5y"], ci_lo["cif_st_5y"], ci_hi["cif_st_5y"]))
cat(sprintf("  RR:        %.3f  (95%% CI: %.3f, %.3f)\n",
            med_est["rr_5y"], ci_lo["rr_5y"], ci_hi["rr_5y"]))
cat(sprintf("  RD:        %+.2f%%  (95%% CI: %+.2f%%, %+.2f%%)\n\n",
            med_est["rd_5y"], ci_lo["rd_5y"], ci_hi["rd_5y"]))

cat("AJ-CIF by time point:\n")
tp_names <- c("3m","6m","1y","2y","3y","5y")
for (i in 1:6) {
  nu <- med_est[paste0("cifNU_", tp_names[i])]
  st <- med_est[paste0("cifST_", tp_names[i])]
  rr <- ifelse(nu > 0, st/nu, NA)
  cat(sprintf("  %3s: NU=%.1f%%  ST=%.1f%%  RR=%.3f\n", tp_names[i], nu, st, rr))
}

# 保存结果表
res_tab <- data.frame(
  metric = c("Cause-specific HR", "5-year CIF: NON_USER (%)", "5-year CIF: STATIN (%)",
             "5-year RR", "5-year RD (%)"),
  est_median = c(med_est["hr_dm"], med_est["cif_nu_5y"], med_est["cif_st_5y"],
                 med_est["rr_5y"], med_est["rd_5y"]),
  ci_lo = c(ci_lo["hr_dm"], ci_lo["cif_nu_5y"], ci_lo["cif_st_5y"],
            ci_lo["rr_5y"], ci_lo["rd_5y"]),
  ci_hi = c(ci_hi["hr_dm"], ci_hi["cif_nu_5y"], ci_hi["cif_st_5y"],
            ci_hi["rr_5y"], ci_hi["rd_5y"]),
  stringsAsFactors = FALSE
)
write.csv(res_tab, "ap_dm_results_final.csv", row.names = FALSE)

# ============================================================
# 10. CIF 图
# ============================================================
time_days <- c(90, 180, 365, 730, 1095, 1825)
time_months <- time_days / 30.4375

cif_df <- data.frame(
  month = rep(time_months, 2),
  arm = rep(c("NON_USER", "STATIN"), each = 6),
  est = c(med_est[paste0("cifNU_", tp_names)], med_est[paste0("cifST_", tp_names)]),
  lo = c(ci_lo[paste0("cifNU_", tp_names)], ci_lo[paste0("cifST_", tp_names)]),
  hi = c(ci_hi[paste0("cifNU_", tp_names)], ci_hi[paste0("cifST_", tp_names)])
)

p <- ggplot(cif_df, aes(month, est, color = arm, fill = arm)) +
  geom_ribbon(aes(ymin = pmax(lo, 0), ymax = pmin(hi, 50)), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#5A5A5A", "#E8A33D")) +
  scale_fill_manual(values = c("#5A5A5A", "#E8A33D")) +
  scale_x_continuous(breaks = seq(0, 60, 12), limits = c(0, 62)) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "Months after AP", y = "Cumulative incidence of diabetes (%)",
       title = "Statin vs Non-user — New-onset Diabetes after AP",
       subtitle = "ATT, Entropy Balanced, Aalen-Johansen CIF, 5×200 bootstrap",
       color = "Arm", fill = "Arm") +
  theme_bw() +
  theme(plot.title = element_text(size = 11))

ggsave("ap_dm_cif_final.png", p, width = 7, height = 5, dpi = 150)

cat("\nOutput files:\n")
cat("  ap_dm_results_final.csv   — results table\n")
cat("  ap_dm_cif_final.png       — CIF plot\n")
cat("  ap_dm_bootstrap_final.rds — bootstrap matrix\n")
cat("  imp_ap_dm_final.rds       — MICE imputations\n")
cat("  bal_stats_final.rds       — balance diagnostics\n")
cat("\nDone.\n")
