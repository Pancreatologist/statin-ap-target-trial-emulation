# ============================================================================
# S3 (fresh weightit): Fibrate vs Statin, active comparator
# focal = FIBRATE (small arm n=2386), weight to FIB
# ============================================================================
library(mice); library(WeightIt); library(dplyr); library(data.table); library(splines)

set.seed(20260625)   # added 2026-09-05 T2b: mirrors authoritative seed system (was missing when script was drafted 2026-06-28)

df_s3 <- readRDS("df_s3_statin_fibrate.rds")
imp_s3 <- readRDS("imp_s3_statin_fibrate.rds")
if (!file.exists("imp_s3_statin_fibrate.rds")) {
  cat("ERROR: imp_s3_statin_fibrate.rds not found!\n")
  quit("no")
}

ps_vars_all <- c("age","male","white","hispanic","smoke_ever","alcohol",
                  "ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
                  "n_inpatient_ed","n_outpatient","n_visits_total",
                  "htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
                  "dementia","as_spond","ra","sle","pulmonary","cancer","neuro",
                  "aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
                  "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
                  "antithromb","nsaid")

# CIF计算函数 (FIBRATE vs STATIN)
calc_cif <- function(f_ap, f_d, n_months = 12) {
  out <- list()
  for (a in c("STATIN", "FIBRATE")) {
    nd <- data.frame(arm = factor(a, levels = c("FIBRATE", "STATIN")), month = 1:n_months)
    h_ap <- predict(f_ap, nd, type = "response")
    h_d  <- predict(f_d,  nd, type = "response")
    S <- 1
    cif_ap <- numeric(n_months); cif_d <- numeric(n_months)
    for (j in 1:n_months) {
      cif_ap[j] <- (if (j > 1) cif_ap[j-1] else 0) + h_ap[j] * S
      cif_d[j]  <- (if (j > 1) cif_d[j-1]  else 0) + h_d[j]  * S
      S <- S * (1 - h_ap[j] - h_d[j])
    }
    out[[a]] <- list(ap = cif_ap, d = cif_d)
  }
  out
}

n_imp   <- imp_s3$m
n_boot  <- 200
total   <- n_imp * n_boot
cat(sprintf("S3 fresh weightit: %d imp x %d boot = %d reps\n", n_imp, n_boot, total))

col_names <- c(paste0("apFI_",1:12), paste0("apST_",1:12), "rr_ap", "rr_death")
boot_s3 <- matrix(NA_real_, nrow = total, ncol = length(col_names))
colnames(boot_s3) <- col_names

start_time <- Sys.time()
rep_idx <- 0

for (imp_i in 1:n_imp) {
  d <- as.data.table(complete(imp_s3, imp_i))
  d <- d[!is.na(time_days) & time_days > 0, ]
  d[, arm := factor(arm, levels = c("FIBRATE", "STATIN"))]
  for (v in c(ps_vars_all, "ap_event", "death_event", "time_days")) {
    if (v %in% names(d) && inherits(d[[v]], "integer64")) d[, (v) := as.numeric(get(v))]
  }

  cat(sprintf("Imp %d: N=%d (FI=%d, ST=%d)\n", imp_i, nrow(d),
              sum(d$arm=="FIBRATE"), sum(d$arm=="STATIN")))

  for (b in 1:n_boot) {
    rep_idx <- rep_idx + 1
    idx_b <- sample(nrow(d), replace = TRUE)
    d_b <- d[idx_b]

    tryCatch({
      w_b <- weightit(as.formula(paste("arm ~", paste(ps_vars_all, collapse = "+"))),
                      data = d_b, method = "ebal", estimand = "ATT", focal = "FIBRATE")
      d_b[, w := as.numeric(w_b$weights)]
    }, error = function(e) { d_b[, w := 1.0] })

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

    rr_ap <- rr_death <- NA_real_
    fi_cif <- st_cif <- numeric(12)

    if (nrow(pm) > 100 && sum(pm$ap_m) > 3) {
      tryCatch({
        f_ap <- glm(ap_m ~ arm * ns(month, df=2, Boundary.knots=c(1,12)),
                    family = binomial(), weights = w, data = pm)
        f_d  <- glm(death_m ~ arm * ns(month, df=2, Boundary.knots=c(1,12)),
                    family = binomial(), weights = w, data = pm)
        cc <- calc_cif(f_ap, f_d)
        fi_cif <- cc$FIBRATE$ap
        st_cif <- cc$STATIN$ap
        rr_ap    <- fi_cif[12] / max(st_cif[12], 1e-10)
        rr_death <- cc$FIBRATE$d[12] / max(cc$STATIN$d[12], 1e-10)
      }, error = function(e) {})
    }

    boot_s3[rep_idx, paste0("apFI_", 1:12)] <- fi_cif * 1e5
    boot_s3[rep_idx, paste0("apST_", 1:12)] <- st_cif * 1e5
    boot_s3[rep_idx, "rr_ap"]    <- rr_ap
    boot_s3[rep_idx, "rr_death"] <- rr_death

    if (b %% 50 == 0) {
      cat(sprintf("  Imp %d / Rep %d | RR_ap=%.3f\n", imp_i, b, rr_ap))
    }
  }
}

cat(sprintf("\nDone in %.1f min\n", difftime(Sys.time(), start_time, units = "mins")))

rr_v <- boot_s3[, "rr_ap"]
rr_v <- rr_v[is.finite(rr_v) & rr_v > 0 & rr_v < 100]
cat(sprintf("S3 Fibrate RR: median=%.3f, 95%%CI=(%.3f, %.3f), valid=%d/%d\n",
            median(rr_v), quantile(rr_v,0.025), quantile(rr_v,0.975), length(rr_v), total))

fi12 <- boot_s3[, "apFI_12"]; fi12 <- fi12[is.finite(fi12)]
st12 <- boot_s3[, "apST_12"]; st12 <- st12[is.finite(st12)]
cat(sprintf("12-month CIF: Fibrate=%.1f, Statin=%.1f (per 100k)\n",
            median(fi12), median(st12)))

saveRDS(boot_s3, "s3_bootstrap_fresh_weightit.rds")
cat("Saved s3_bootstrap_fresh_weightit.rds\n")