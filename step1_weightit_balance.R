# ============================================================================
# Step 1: Fresh weightit ebal across 5 MICE imputations
# 验证: 50/50 covariates SMD=0 (NU weighted mean = ST mean)
# ============================================================================
library(mice); library(WeightIt); library(dplyr); library(data.table)

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

calc_smd <- function(d, w = NULL, vars = ps_vars) {
  trt <- d$arm == "STATIN"; ctrl <- d$arm == "NON_USER"
  smd_before <- smd_after <- numeric(length(vars))
  for (i in seq_along(vars)) {
    v <- vars[i]
    x <- as.numeric(d[[v]])
    x_t <- x[trt]; x_c <- x[ctrl]
    mu_t <- mean(x_t, na.rm = T); mu_c <- mean(x_c, na.rm = T)
    var_t <- var(x_t, na.rm = T); var_c <- var(x_c, na.rm = T)
    sd_pooled <- sqrt((var_t + var_c) / 2)
    smd_before[i] <- abs(mu_t - mu_c) / sd_pooled
    if (!is.null(w)) {
      w_t <- w[trt]; w_c <- w[ctrl]
      mu_t_w <- weighted.mean(x_t, w_t, na.rm = T)
      mu_c_w <- weighted.mean(x_c, w_c, na.rm = T)
      smd_after[i] <- abs(mu_t_w - mu_c_w) / sd_pooled
    }
  }
  data.frame(variable = vars, smd_before = smd_before, smd_after = smd_after,
             stringsAsFactors = FALSE)
}

# ---- 跑5个imputation的ebal ----
cat("=== Running fresh weightit ebal on 5 imputations ===\n")
weightit_models <- vector("list", 5)
smd_list <- vector("list", 5)

for (imp_i in 1:5) {
  cat(sprintf("\n--- Imputation %d ---\n", imp_i))
  d <- as.data.table(complete(imp, imp_i))
  d <- d[person_id %in% frozen_ids, ]
  d[, arm := factor(arm, levels = c("NON_USER", "STATIN"))]
  d[, arm_n := as.numeric(arm == "STATIN")]

  # integer64 -> numeric
  for (v in c(ps_vars, "ap_event", "death_event", "time_days")) {
    if (v %in% names(d) && inherits(d[[v]], "integer64")) {
      d[, (v) := as.numeric(get(v))]
    }
  }

  cat(sprintf("N=%d, ST=%d, NU=%d\n", nrow(d), sum(d$arm=="STATIN"), sum(d$arm=="NON_USER")))

  # weightit ebal (ATT)
  fml <- as.formula(paste("arm_n ~", paste(ps_vars, collapse = " + ")))
  w_fit <- weightit(fml, data = d, method = "ebal", estimand = "ATT")
  weightit_models[[imp_i]] <- w_fit
  d[, w := as.numeric(w_fit$weights)]

  # balance check
  smd <- calc_smd(d, d$w)
  smd_list[[imp_i]] <- smd

  n_gt_01 <- sum(abs(smd$smd_after) > 0.1, na.rm = TRUE)
  max_smd <- max(abs(smd$smd_after), na.rm = TRUE)
  cat(sprintf("SMD after: max=%.4f, >0.1: %d/50\n", max_smd, n_gt_01))

  # 检查关键变量
  for (v in c("n_inpatient_ed", "n_outpatient", "n_visits_total")) {
    row <- smd[smd$variable == v, ]
    d_v <- as.numeric(d[[v]])
    mu_st_raw <- mean(d_v[d$arm=="STATIN"], na.rm=T)
    mu_nu_raw <- mean(d_v[d$arm=="NON_USER"], na.rm=T)
    mu_nu_w <- weighted.mean(d_v[d$arm=="NON_USER"], d$w[d$arm=="NON_USER"], na.rm=T)
    cat(sprintf("  %s: ST_raw=%.3f, NU_raw=%.3f, NU_weighted=%.3f, SMD_after=%.4f\n",
                v, mu_st_raw, mu_nu_raw, mu_nu_w, row$smd_after))
  }
}

# ---- 汇总5个imputation的平均SMD ----
cat("\n=== Average SMD across 5 imputations ===\n")
smd_avg <- smd_list[[1]][, c("variable", "smd_before")]
smd_after_mat <- sapply(smd_list, function(s) s$smd_after)
smd_avg$smd_after <- rowMeans(smd_after_mat)
smd_avg <- smd_avg[order(-smd_avg$smd_before), ]
smd_avg$variable <- factor(smd_avg$variable, levels = smd_avg$variable)

n_gt_01_avg <- sum(abs(smd_avg$smd_after) > 0.1, na.rm = TRUE)
max_smd_avg <- max(abs(smd_avg$smd_after), na.rm = TRUE)
cat(sprintf("Average max SMD after: %.4f\n", max_smd_avg))
cat(sprintf("Average n SMD > 0.1: %d / 50\n", n_gt_01_avg))

# ---- 保存 ----
saveRDS(list(
  weightit_models = weightit_models,
  smd_list = smd_list,
  smd_avg = smd_avg
), "weightit_ebal_fresh.rds")

cat("\nSaved to weightit_ebal_fresh.rds\n")