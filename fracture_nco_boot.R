# Fracture NCO bootstrap CI — 后台运行
library(dplyr); library(mice); library(WeightIt); library(splines)

df_master <- readRDS("df_master.rds")
df_nco    <- readRDS("df_nco.rds")
imp       <- readRDS("imp_statin_htg_v4.rds")
STEP      <- 365 / 12

df <- df_master %>% inner_join(df_nco %>% select(person_id, t0, fracture_date, mva_injury_date),
                               by = "person_id", suffix = c("", "_nco"))
df$t0_nco <- as.Date(df$t0)
df <- df %>% mutate(
  ap_date_d    = as.Date(ap_date),
  death_date_d = as.Date(death_date),
  obs_end_d    = as.Date(obs_end),
  end_date     = pmin(coalesce(ap_date_d, as.Date("2099-01-01")),
                     coalesce(death_date_d, as.Date("2099-01-01")),
                     obs_end_d, t0_nco + 365, na.rm = TRUE)
)
df$frac_event <- as.integer(!is.na(df$fracture_date) & df$fracture_date > df$t0_nco &
                           df$fracture_date <= pmin(df$end_date, df$t0_nco + 365, na.rm = TRUE))

demo_vars   <- c("age", "male", "white", "hispanic", "smoke_ever")
cont_vars   <- c("ldl", "hdl", "tc", "tg", "bmi", "calcium", "alt", "egfr",
                 "n_inpatient_ed", "n_outpatient", "n_visits_total")
comorb_cols <- c("htn", "obesity", "diabetes", "chd", "stroke", "pvd", "afib", "kidney",
                 "dementia", "as_spond", "ra", "sle", "pulmonary", "cancer", "neuro")
drug_cols   <- c("aspirin", "insulin", "metformin", "su", "tzd", "glinide", "agi", "dpp4", "glp1",
                 "sglt2", "acei_arb", "bblocker", "ccb", "diuretic", "alpha_blk", "other_htn",
                 "antithromb", "nsaid")
ps_vars_all <- c(demo_vars, "alcohol", cont_vars, comorb_cols, drug_cols)

boot_one_imp <- function(di, B = 50, seed = 42) {
  set.seed(seed); n <- nrow(di); rr <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample(n, replace = TRUE); db <- di[idx, ]
    ps_vars <- intersect(ps_vars_all, names(db))
    w <- tryCatch(weightit(reformulate(ps_vars, "arm"), data = db, method = "ebal",
                           estimand = "ATT", focal = "STATIN"), error = function(e) NULL)
    if (is.null(w) || sum(is.finite(w$weights)) == 0) { rr[b] <- NA; next }
    db$w <- w$weights; db <- db[!is.na(db$w) & is.finite(db$w), ]
    if (nrow(db) < 100) { rr[b] <- NA; next }
    db <- db %>% mutate(t = pmin(time_days, 365), ev = as.integer(frac_event == 1 & time_days <= 365),
                        M = pmin(pmax(1L, ceiling(t / STEP)), 12L)) %>% filter(M > 0)
    rid_v <- rep(db$rid, db$M); arm_v <- rep(db$arm, db$M); w_v <- rep(db$w, db$M)
    ev_v  <- rep(db$ev, db$M); m_v <- unlist(lapply(db$M, seq_len), use.names = FALSE)
    long <- data.frame(rid = rid_v, arm = arm_v, w = w_v, ev = ev_v, month = m_v)
    long <- long %>% group_by(rid) %>% mutate(last_m = month == max(month)) %>% ungroup()
    long$ev_m <- as.integer(long$ev == 1 & long$last_m)
    if (sum(long$ev_m, na.rm = TRUE) < 3) { rr[b] <- NA; next }
    f <- glm(ev_m ~ arm * ns(month, df = 2, Boundary.knots = c(1, 12)),
             family = binomial(), weights = w, data = long)
    cif <- function(a) { nd <- data.frame(arm = factor(a, levels = levels(db$arm)), month = 1:12)
      h <- predict(f, nd, type = "response"); S <- 1; ci <- 0
      for (j in 1:12) { ci <- ci + h[j] * S; S <- S * (1 - h[j]) }; ci }
    rr[b] <- cif("STATIN") / cif("NON_USER") * 1e5 / 1e5
    if (b %% 10 == 0) cat(sprintf("  Boot %d/%d\n", b, B))
  }
  rr
}

m <- 5; B <- 50; all_rr <- list()
for (i in seq_len(m)) {
  cat(sprintf("\n=== Imp %d/%d ===\n", i, m))
  di <- complete(imp, i) %>% select(person_id, arm, one_of(c(demo_vars, cont_vars, comorb_cols, drug_cols, "alcohol"))) %>%
    inner_join(df %>% select(person_id, t0_nco, end_date, frac_event, time_days), by = "person_id")
  di$arm <- factor(di$arm, levels = c("NON_USER", "STATIN")); di$rid <- seq_len(nrow(di))
  di <- di[!is.na(di$time_days) & di$time_days > 0 & !is.na(di$arm), ]
  all_rr[[i]] <- boot_one_imp(di, B = B, seed = 42 + i)
  cat(sprintf("  Imp %d: median=%.3f CI=(%.3f,%.3f)\n", i,
              median(all_rr[[i]], na.rm=T), quantile(all_rr[[i]], 0.025, na.rm=T), quantile(all_rr[[i]], 0.975, na.rm=T)))
}

all_v <- unlist(all_rr); all_v <- all_v[!is.na(all_v)]
ci_df <- data.frame(
  outcome = "骨折/跌倒 (Fracture)",
  RR_median = median(all_v),
  RR_ci_lo  = quantile(all_v, 0.025),
  RR_ci_hi  = quantile(all_v, 0.975),
  n_rep = length(all_v),
  crosses_null = quantile(all_v, 0.025) < 1 & quantile(all_v, 0.975) > 1
)
cat("\n=== FINAL ===\n"); print(ci_df)
saveRDS(ci_df, "fracture_nco_ci.rds")
write.csv(ci_df, "fracture_nco_ci.csv", row.names = FALSE)
cat("Done.\n")
