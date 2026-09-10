# ============================================================================
# 投稿版图表：大字 + 全黑 + PDF/PNG双输出
# ----------------------------------------------------------------------------
# AUTHORITATIVE sources only (2026-09-05 revision):
#   - main results:   bootstrap_fresh_weightit.rds (step2_bootstrap_v2.R)
#   - S1/S2 weights:  weightit_ebal_fresh.rds (NOT watt_* = LEGACY)
#   - grace 90/365:   grace_fresh_weightit_summary.csv / .rds (weights re-estimated
#                     per bootstrap rep; grace_weightit_results.rds = LEGACY)
#   - fracture NCO:   fracture_nco_fresh_summary.csv (frozen cohort;
#                     fracture_nco_ci.rds on the full master cohort = LEGACY)
#                     fibrate analysis CANCELLED (2026-09-06): no row, no S7.
#   - Figure S5:      cause-specific curves now EXACTLY reconstructed from stored
#                     monthly AJ CIF columns (validated vs rr_apcs, max diff 1e-13;
#                     revision_T2a_test_cs_recon.R) -- replaces ratio-scaled
#                     derived curves.
#   - Figure 3:       NO fabricated error bars; all shown CIs are real. MVA row
#                     remains NA (not estimable, 2 events, as in manuscript).
# ============================================================================
library(dplyr); library(tidyr); library(ggplot2); library(scales); library(grid)
library(mice); library(data.table); library(splines); library(stringr)
library(gridExtra)   # for tableGrob (number-at-risk)

col_statin  <- "#1565C0"
col_nonuser <- "#EF6C00"
col_exclude <- "#EF9A9A"
col_total   <- "#90CAF9"

# ---- 统一主题 ----
theme_fig <- function(base_size = 14) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text = element_text(color = "black", size = base_size * 0.85),
      axis.title = element_text(color = "black", size = base_size, face = "plain"),
      legend.text = element_text(color = "black", size = base_size * 0.85),
      legend.title = element_text(color = "black", size = base_size * 0.85, face = "plain"),
      strip.text = element_text(color = "black", size = base_size * 0.9, face = "plain"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank()
    )
}

theme_bw_fig <- function(base_size = 13) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = base_size * 0.85),
      axis.title = element_text(color = "black", size = base_size),
      legend.text = element_text(color = "black", size = base_size * 0.85),
      legend.title = element_text(color = "black", size = base_size * 0.85),
      strip.text = element_text(color = "black", size = base_size * 0.9),
      strip.background = element_rect(color = "black", linewidth = 0.5, fill = "grey95")
    )
}

# 统一保存函数
save_fig <- function(name, p, w, h) {
  ggsave(sprintf("%s.pdf", name), p, width = w, height = h)
  ggsave(sprintf("%s.png", name), p, width = w, height = h, dpi = 600)
  cat(sprintf("Saved: %s.pdf/.png\n", name))
}

# ============================================================================
# 1. 加载数据
# ============================================================================
boot_main <- readRDS("bootstrap_fresh_weightit.rds")
imp       <- readRDS("imp_statin_htg_v4.rds")
frozen     <- readRDS("frozen_cohort_ids.rds")
ebal_res  <- readRDS("weightit_ebal_fresh.rds")            # AUTHORITATIVE weights (S1/S2); watt_* = LEGACY, no longer loaded
grace_sm  <- read.csv("grace_fresh_weightit_summary.csv")   # AUTHORITATIVE grace re-run; grace_weightit_results.rds = LEGACY
grace_boot <- readRDS("grace_fresh_weightit.rds")           # monthly CIF replicates (Figure S6)
frac_fresh <- read.csv("fracture_nco_fresh_summary.csv")    # AUTHORITATIVE re-run (frozen cohort); fracture_nco_ci.rds = LEGACY
df_master  <- readRDS("df_master.rds")
# fibrate active-comparator analysis CANCELLED (2026-09-06 author ruling):
# replaced by descriptive T5 counts table; no fibrate rows/figures anywhere.

frozen_ids <- frozen$person_id
n_total_frozen <- length(frozen_ids)
n_st_frozen    <- sum(frozen$arm == "STATIN")
n_nu_frozen    <- sum(frozen$arm == "NON_USER")

# ============================================================================
# 2. Bootstrap 汇总
# ============================================================================
bs <- data.frame(
  point = apply(boot_main, 2, median, na.rm = TRUE),
  lo    = apply(boot_main, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
  hi    = apply(boot_main, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
)
rownames(bs) <- colnames(boot_main)

ap_rr  <- bs["rr_ap", "point"];    ap_lo  <- bs["rr_ap", "lo"];    ap_hi  <- bs["rr_ap", "hi"]
d_rr   <- bs["rr_death", "point"]; d_lo   <- bs["rr_death", "lo"]; d_hi   <- bs["rr_death", "hi"]

# Grace (AUTHORITATIVE re-run: weights re-estimated per bootstrap rep)
g90  <- grace_sm[grace_sm$grace_days == 90, ]
g365 <- grace_sm[grace_sm$grace_days == 365, ]
stopifnot(nrow(g90) == 1, nrow(g365) == 1)
g90_rr  <- g90$rr_ap_median;  g90_lo  <- g90$rr_ap_ci_lo;  g90_hi  <- g90$rr_ap_ci_hi
g365_rr <- g365$rr_ap_median; g365_lo <- g365$rr_ap_ci_lo; g365_hi <- g365$rr_ap_ci_hi

# Fracture NCO (AUTHORITATIVE re-run on frozen cohort)
frac_rr <- frac_fresh$rr_median; frac_lo <- frac_fresh$rr_ci_lo; frac_hi <- frac_fresh$rr_ci_hi

# ============================================================================
# Number at risk: compute from imp[1] frozen cohort (unweighted, per arm)
# M = months contributed = pmin(ceiling(pmin(time_days,365)/(365/12)), 12)
# at-risk at start of month k = sum(M >= k) per arm
# ============================================================================
d_nr <- complete(imp, 1)
d_nr <- d_nr[d_nr$person_id %in% frozen_ids, ]
d_nr$arm <- factor(d_nr$arm, levels = c("NON_USER", "STATIN"))
d_nr$t_trunc <- pmin(as.numeric(d_nr$time_days), 365)
d_nr$M <- pmin(pmax(1L, ceiling(d_nr$t_trunc / (365/12))), 12L)
nr_st <- sapply(1:12, function(k) sum(d_nr$M[d_nr$arm == "STATIN"] >= k))
nr_nu <- sapply(1:12, function(k) sum(d_nr$M[d_nr$arm == "NON_USER"] >= k))
rm(d_nr); gc()

# Helper: add number-at-risk table below a ggplot CIF figure
add_nrisk <- function(p, nr_st, nr_nu, y_scale_mult = 0.08) {
  nr_df <- data.frame(
    month = 1:12,
    statin  = nr_st,
    nonuser = nr_nu
  )
  # Build text rows for at-risk table
  # Placed below x-axis using annotation_custom with tableGrob
  tg <- gridExtra::tableGrob(
    rbind(c("Month", as.character(1:12)),
          c("Statin initiators", format(nr_st, big.mark = ",")),
          c("Non-initiators", format(nr_nu, big.mark = ","))),
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core = list(fg_params = list(cex = 0.7, col = "black")),
      colhead = list(fg_params = list(cex = 0.7, col = "black", fontface = "bold"))
    )
  )
  # Add padding
  tg$widths <- unit(rep(1/13, 13), "npc")
  p + annotation_custom(grob = tg,
                        xmin = 0.5, xmax = 12.5,
                        ymin = -Inf, ymax = 0) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 6, r = 10, b = 80, l = 10))
}

# ============================================================================
# FIGURE 1: Cohort Flow (STROBE-compliant with all exclusion steps)
# Upstream exclusion numbers (prior statin/AP/cancer) not separately tracked;
# shown as combined box per study_design_for_gemini.md L44
# ============================================================================
cat("\n=== Figure 1 ===\n")
total_aou <- 633000   # All of Us total participants (>633k, approximate)
total_htg_identified <- NA   # HTG identified before upstream exclusions (not tracked)
total_cohort <- nrow(df_master)   # master cohort, after upstream exclusions
no_follow <- nrow(df_master) - n_total_frozen   # excluded for no follow-up
final_analytical <- n_total_frozen   # final analytical cohort
final_st <- n_st_frozen
final_nu <- n_nu_frozen

# STROBE flow layout: 6 levels, exclusion boxes on right
y_aou    <- 12;   h_aou    <- 0.8   # All of Us
y_htg_id <- 10.3; h_htg_id <- 0.8   # HTG identified
y_up_ex  <- 8.4;  h_up_ex  <- 1.6   # Upstream exclusions box (right)
y_htg_ex <- 7.5;  h_htg_ex <- 0.8   # HTG after upstream exclusions
y_nofo   <- 5.6;  h_nofo   <- 1.2   # No follow-up exclusion box (right)
y_final  <- 4.5;  h_final  <- 0.9   # Final analytical cohort
y_arms   <- 2;    h_arms   <- 1.0   # Two arms

p1 <- ggplot() +
  # Level 1: All of Us
  annotate("rect", xmin = 0.20, xmax = 0.70, ymin = y_aou-h_aou/2, ymax = y_aou+h_aou/2,
           fill = col_total, color = "black", linewidth = 1.2) +
  annotate("text", x = 0.45, y = y_aou,
           label = "All of Us Research Program\n(n > 633,000)",
           size = 4.5, fontface = "bold", vjust = 0.5, color = "black") +
  # Arrow Aou -> HTG identified
  annotate("segment", x = 0.45, xend = 0.45, y = y_aou-h_aou/2, yend = y_htg_id+h_htg_id/2+0.05,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
           linewidth = 1, color = "black") +
  # Level 2: HTG identified
  annotate("rect", xmin = 0.20, xmax = 0.70, ymin = y_htg_id-h_htg_id/2, ymax = y_htg_id+h_htg_id/2,
           fill = "#BBDEFB", color = "black", linewidth = 1.2) +
  annotate("text", x = 0.45, y = y_htg_id,
           label = "HTG cohort identified\n(TG \u2265 500 mg/dL or ICD code)",
           size = 4.2, fontface = "bold", vjust = 0.5, color = "black") +
  # Arrow HTG identified -> upstream exclusions (right)
  annotate("segment", x = 0.70, xend = 0.70, y = y_htg_id-h_htg_id/2-0.02, yend = y_up_ex+h_up_ex/2,
           linewidth = 0.8, color = "black", linetype = "dashed") +
  annotate("segment", x = 0.70, xend = 0.82, y = y_up_ex+h_up_ex/2, yend = y_up_ex+h_up_ex/2,
           arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
           linewidth = 0.8, color = "black") +
  # Upstream exclusions box
  annotate("rect", xmin = 0.82, xmax = 0.99, ymin = y_up_ex-h_up_ex/2, ymax = y_up_ex+h_up_ex/2,
           fill = col_exclude, color = "black", linewidth = 1.2) +
  annotate("text", x = 0.905, y = y_up_ex+0.55,
           label = "Excluded (before index):", size = 3.5, fontface = "bold",
           color = "black", vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 0.905, y = y_up_ex+0.2,
           label = "\u2022 Prior statin use\n\u2022 Prior acute/chronic pancreatitis\n\u2022 Prior pancreatic cancer",
           size = 3.2, color = "black", vjust = 0.5, hjust = 0.5, lineheight = 0.9) +
  annotate("text", x = 0.905, y = y_up_ex-0.55,
           label = "(numbers not separately tracked)",
           size = 2.8, color = "black", vjust = 0.5, hjust = 0.5, fontface = "italic") +
  # Arrow HTG identified -> HTG after upstream exclusions
  annotate("segment", x = 0.45, xend = 0.45, y = y_htg_id-h_htg_id/2, yend = y_htg_ex+h_htg_ex/2+0.05,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
           linewidth = 1, color = "black") +
  # Level 3: HTG after upstream exclusions
  annotate("rect", xmin = 0.20, xmax = 0.70, ymin = y_htg_ex-h_htg_ex/2, ymax = y_htg_ex+h_htg_ex/2,
           fill = "#E3F2FD", color = "black", linewidth = 1.2) +
  annotate("text", x = 0.45, y = y_htg_ex,
           label = sprintf("HTG cohort after upstream exclusions\n(n = %s)", scales::comma(total_cohort)),
           size = 4.2, fontface = "bold", vjust = 0.5, color = "black") +
  # Arrow HTG after exclusions -> no follow-up exclusion (right)
  annotate("segment", x = 0.70, xend = 0.70, y = y_htg_ex-h_htg_ex/2-0.02, yend = y_nofo+h_nofo/2,
           linewidth = 0.8, color = "black", linetype = "dashed") +
  annotate("segment", x = 0.70, xend = 0.82, y = y_nofo+h_nofo/2, yend = y_nofo+h_nofo/2,
           arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
           linewidth = 0.8, color = "black") +
  # No follow-up exclusion box
  annotate("rect", xmin = 0.82, xmax = 0.99, ymin = y_nofo-h_nofo/2, ymax = y_nofo+h_nofo/2,
           fill = col_exclude, color = "black", linewidth = 1.2) +
  annotate("text", x = 0.905, y = y_nofo+0.35,
           label = "Excluded:", size = 3.5, fontface = "bold",
           color = "black", vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 0.905, y = y_nofo+0.0,
           label = "no follow-up time\nafter landmark",
           size = 3.2, color = "black", vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 0.905, y = y_nofo-0.38,
           label = sprintf("(n = %s)", scales::comma(no_follow)),
           size = 3.5, fontface = "bold", color = "black", vjust = 0.5, hjust = 0.5) +
  # Arrow HTG after exclusions -> Final
  annotate("segment", x = 0.45, xend = 0.45, y = y_htg_ex-h_htg_ex/2, yend = y_final+h_final/2+0.05,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
           linewidth = 1, color = "black") +
  # Level 4: Final analytical cohort
  annotate("rect", xmin = 0.20, xmax = 0.70, ymin = y_final-h_final/2, ymax = y_final+h_final/2,
           fill = "#E3F2FD", color = "black", linewidth = 1.8) +
  annotate("text", x = 0.45, y = y_final,
           label = sprintf("Final analytical cohort\n(n = %s)", scales::comma(final_analytical)),
           size = 4.5, fontface = "bold", vjust = 0.5, color = "black") +
  # Arrows Final -> two arms
  annotate("segment", x = 0.35, xend = 0.25, y = y_final-h_final/2, yend = y_arms+h_arms/2+0.05,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
           linewidth = 1, color = col_statin) +
  annotate("segment", x = 0.55, xend = 0.65, y = y_final-h_final/2, yend = y_arms+h_arms/2+0.05,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
           linewidth = 1, color = col_nonuser) +
  annotate("text", x = 0.25, y = (y_final-h_final/2 + y_arms+h_arms/2)/2 + 0.2,
           label = "Statin\ninitiators", size = 3.5, color = col_statin, fontface = "bold") +
  annotate("text", x = 0.65, y = (y_final-h_final/2 + y_arms+h_arms/2)/2 + 0.2,
           label = "Non-\ninitiators", size = 3.5, color = col_nonuser, fontface = "bold") +
  # Level 5: Two arms
  annotate("rect", xmin = 0.10, xmax = 0.40, ymin = y_arms-h_arms/2, ymax = y_arms+h_arms/2,
           fill = col_statin, color = "black", linewidth = 1.2) +
  annotate("text", x = 0.25, y = y_arms,
           label = sprintf("Statin initiation\n(n = %s)", scales::comma(final_st)),
           size = 4.0, fontface = "bold", color = "white", vjust = 0.5) +
  annotate("rect", xmin = 0.50, xmax = 0.80, ymin = y_arms-h_arms/2, ymax = y_arms+h_arms/2,
           fill = col_nonuser, color = "black", linewidth = 1.2) +
  annotate("text", x = 0.65, y = y_arms,
           label = sprintf("Non-initiation\n(n = %s)", scales::comma(final_nu)),
           size = 4.0, fontface = "bold", color = "white", vjust = 0.5) +
  scale_x_continuous(limits = c(0, 1.0)) +
  scale_y_continuous(limits = c(1, 13)) +
  theme_void(base_size = 14) +
  theme(plot.margin = margin(t = 12, r = 15, b = 12, l = 10))

save_fig("Figure 1", p1, 9, 12)

# ============================================================================
# FIGURE 2: AP CIF
# ============================================================================
cat("\n=== Figure 2 ===\n")
cif_df <- data.frame(
  month = 1:12,
  arm = factor(rep(c("STATIN", "NON_USER"), each = 12),
               levels = c("STATIN", "NON_USER")),
  est = c(bs[paste0("apST_", 1:12), "point"], bs[paste0("apNU_", 1:12), "point"]),
  lo  = c(bs[paste0("apST_", 1:12), "lo"],    bs[paste0("apNU_", 1:12), "lo"]),
  hi  = c(bs[paste0("apST_", 1:12), "hi"],    bs[paste0("apNU_", 1:12), "hi"])
)

p2 <- ggplot(cif_df, aes(x = month, y = est, color = arm, fill = arm)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  scale_color_manual(
    values = c("STATIN" = col_statin, "NON_USER" = col_nonuser),
    labels = c("Statin initiators", "Non-initiators")) +
  scale_fill_manual(
    values = c("STATIN" = col_statin, "NON_USER" = col_nonuser),
    labels = c("Statin initiators", "Non-initiators")) +
  scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, NA)) +
  labs(x = "Month", y = "Cumulative incidence\nper 100 000", color = "", fill = "") +
  theme_fig(14) +
  theme(legend.position = "bottom")

p2 <- add_nrisk(p2, nr_st, nr_nu)
save_fig("Figure 2", p2, 8, 7)

# ============================================================================
# FIGURE 3: Forest Plot
# ============================================================================
cat("\n=== Figure 3 ===\n")

simba_row <- data.frame(
  label = "SIMBA RCT (ITT): acute pancreatitis",
  rr = 1.07, rr_lo = 0.43, rr_hi = 2.66,
  group = "RCT", has_ci = TRUE, stringsAsFactors = FALSE
)

# Rows built conditionally; all shown CIs are REAL bootstrap CIs (no fabrication).
# MVA: not estimable (2 events) -> rr/CI stay NA, point and error bar simply not drawn.
# Fibrate row DELETED (2026-09-06 author ruling): active-comparator analysis
# cancelled, replaced by descriptive T5 counts table.
forest_rows <- list()
forest_rows[["g365"]] <- data.frame(
  label = "  Grace period: 365 days",
  rr = g365_rr, rr_lo = g365_lo, rr_hi = g365_hi,
  group = "Grace period", has_ci = TRUE, stringsAsFactors = FALSE)
forest_rows[["g90"]] <- data.frame(
  label = "  Grace period: 90 days",
  rr = g90_rr, rr_lo = g90_lo, rr_hi = g90_hi,
  group = "Grace period", has_ci = TRUE, stringsAsFactors = FALSE)
forest_rows[["mva"]] <- data.frame(
  label = "Negative control: motor vehicle accident*",
  rr = NA_real_, rr_lo = NA_real_, rr_hi = NA_real_,
  group = "Negative control", has_ci = FALSE, stringsAsFactors = FALSE)
forest_rows[["frac"]] <- data.frame(
  label = "Negative control: fracture/fall",
  rr = frac_rr, rr_lo = frac_lo, rr_hi = frac_hi,
  group = "Negative control", has_ci = TRUE, stringsAsFactors = FALSE)
forest_rows[["death"]] <- data.frame(
  label = "All-cause mortality",
  rr = d_rr, rr_lo = d_lo, rr_hi = d_hi,
  group = "Main", has_ci = TRUE, stringsAsFactors = FALSE)
forest_rows[["primary"]] <- data.frame(
  label = "Primary: Acute pancreatitis",
  rr = ap_rr, rr_lo = ap_lo, rr_hi = ap_hi,
  group = "Main", has_ci = TRUE, stringsAsFactors = FALSE)
forest_data <- do.call(rbind, forest_rows)

forest_data$label <- factor(forest_data$label, levels = rev(forest_data$label))
forest_data$group <- factor(forest_data$group,
  levels = c("Main", "Negative control", "Grace period", "RCT"))

simba_row$label <- factor(simba_row$label, levels = c(levels(forest_data$label), "SIMBA RCT (ITT): acute pancreatitis"))

all_forest <- rbind(
  forest_data[forest_data$label != "Primary: Acute pancreatitis", ],
  simba_row,
  forest_data[forest_data$label == "Primary: Acute pancreatitis", ]
)
all_forest$label <- factor(all_forest$label)
all_forest$group <- factor(all_forest$group,
  levels = c("Main", "Negative control", "Grace period", "RCT"))
grp_breaks <- c("Main", "Negative control", "Grace period", "RCT")
grp_breaks <- grp_breaks[grp_breaks %in% as.character(unique(all_forest$group))]

x_axis_min <- 0.01; x_axis_max <- 10
all_forest$rr_lo_clipped <- pmax(all_forest$rr_lo, x_axis_min)
all_forest$lo_trimmed <- all_forest$rr_lo < x_axis_min & all_forest$has_ci

p3 <- ggplot(all_forest, aes(x = rr, y = label)) +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.6) +
  geom_errorbar(aes(xmin = rr_lo_clipped, xmax = rr_hi, color = group,
                    linetype = ifelse(has_ci, "solid", "dotted")),
                width = 0.3, linewidth = 0.9) +
  geom_point(aes(color = group, shape = group), size = 3.8, na.rm = TRUE) +
  scale_x_log10(limits = c(x_axis_min, x_axis_max),
                breaks = c(0.01, 0.1, 0.5, 1, 2, 5, 10),
                labels = c("0.01", "0.1", "0.5", "1", "2", "5", "10")) +
  geom_segment(data = subset(all_forest, lo_trimmed),
               aes(x = x_axis_min, xend = x_axis_min * 1.3, y = label, yend = label, color = group),
               arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
               linewidth = 0.9, show.legend = FALSE) +
  scale_color_manual(
    values = c("Main" = col_statin, "Negative control" = "#2E7D32",
               "Grace period" = "#E65100",
               "RCT" = "black"),
    breaks = grp_breaks) +
  scale_shape_manual(
    values = c("Main" = 16, "Negative control" = 17,
               "Grace period" = 15,
               "RCT" = 0),
    breaks = grp_breaks) +
  scale_linetype_identity() +
  labs(x = "Risk ratio (log scale)", y = "", color = "", shape = "") +
  theme_fig(13) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(hjust = 0, color = "black", size = 11),
    axis.text.x = element_text(color = "black")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE))

save_fig("Figure 3", p3, 10, 7)

# ============================================================================
# FIGURE S1: Love Plot
# ============================================================================
cat("\n=== Figure S1 ===\n")

ps_vars <- c("age","male","white","hispanic","smoke_ever","alcohol",
             "ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
             "n_inpatient_ed","n_outpatient","n_visits_total",
             "htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
             "dementia","as_spond","ra","sle","pulmonary","cancer","neuro",
             "aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
             "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
             "antithromb","nsaid")

# AUTHORITATIVE fresh ebal weights (frozen cohort, imp1) -- same construction as
# generate_tables_fresh.R Table 1, so figure and table SMDs are bitwise consistent
d1 <- complete(imp, 1)
d1 <- d1[d1$person_id %in% frozen_ids, ]
d1$arm <- factor(d1$arm, levels = c("NON_USER", "STATIN"))
d1$w <- as.numeric(ebal_res$weightit_models[[1]]$weights)
stopifnot(length(d1$w) == nrow(d1))

calc_smd <- function(d, w = NULL, vars = ps_vars, arm_col = "arm") {
  trt <- d[[arm_col]] == "STATIN"; ctrl <- d[[arm_col]] == "NON_USER"
  smd_before <- rep(NA, length(vars)); smd_after <- rep(NA, length(vars))
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
  data.frame(variable = vars, smd_before = smd_before, smd_after = smd_after)
}

smd_full <- calc_smd(d1, d1$w)
# GATE: fresh ebal must be converged (residual SMD ~ 1e-6; LEGACY watt showed 0.10-0.12)
stopifnot(max(smd_full$smd_after, na.rm = TRUE) < 0.001)
smd_full <- smd_full[order(-smd_full$smd_before), ]
smd_full$variable <- factor(smd_full$variable, levels = smd_full$variable)

covar_groups <- list(
  Demographics = c("age","male","white","hispanic","smoke_ever","alcohol"),
  Lipids = c("ldl","hdl","tc","tg"),
  `Vitals/Renal` = c("bmi","calcium","alt","egfr"),
  `Healthcare utilisation` = c("n_inpatient_ed","n_outpatient","n_visits_total"),
  Comorbidities = c("htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
                     "dementia","as_spond","ra","sle","pulmonary","cancer","neuro"),
  Medications = c("aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
                   "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
                   "antithromb","nsaid")
)
smd_full$group <- "Other"
for (g in names(covar_groups)) {
  smd_full$group[smd_full$variable %in% covar_groups[[g]]] <- g
}

group_colors <- c(
  Demographics = "#E53935", Lipids = "#8E24AA", `Vitals/Renal` = "#1E88E5",
  `Healthcare utilisation` = "#FB8C00", Comorbidities = "#43A047",
  Medications = "#00ACC1", Other = "#78909C"
)

smd_long <- tidyr::pivot_longer(smd_full,
  cols = c(smd_before, smd_after),
  names_to = "time", values_to = "smd")
smd_long$time <- factor(smd_long$time,
  levels = c("smd_before", "smd_after"),
  labels = c("Unadjusted", "Adjusted"))

p4 <- ggplot(smd_full, aes(y = variable)) +
  geom_vline(xintercept = 0.1, color = "black", linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  geom_segment(aes(x = smd_before, xend = smd_after, y = variable, yend = variable),
               color = "grey60", linewidth = 0.5, alpha = 0.8) +
  geom_point(data = smd_long,
             aes(x = smd, color = group, shape = time), size = 2.5, stroke = 0.7) +
  scale_color_manual(values = group_colors, name = "Covariate group") +
  scale_shape_manual(
    values = c("Unadjusted" = 1, "Adjusted" = 16),
    name = "") +
  labs(x = "Absolute standardized mean difference", y = "") +
  theme_bw_fig(11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.y = element_text(color = "black", size = 9)
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 2),
    shape = guide_legend(order = 2, nrow = 1)
  )

save_fig("Figure S1", p4, 10, 12)

# FIGURE S7 (fibrate CIF): DELETED 2026-09-06 -- active-comparator analysis
# cancelled by author ruling; fibrate appears nowhere in figures or tables.

# ============================================================================
# FIGURE S5: Competing Risk vs Cause-Specific
# Cause-specific curves read DIRECTLY from the T2g model output
# (figs5_cs_monthly_summary.csv; revision_T2g_figs5_cs_curves.R, 2026-09-05):
# within every bootstrap replicate the cause-specific model was refitted and
# its monthly CIF stored -- no reconstruction, no ratio scaling. The T2g run
# reproduced all 75 authoritative bootstrap columns with max_abs_diff = 0
# (figs5_cs_consistency_check.csv); rr_apcs reproduced the authoritative
# bootstrap columns exactly. Resulting RR values are study results and are
# not embedded here; they are read from the summary CSV below.
# ============================================================================
cat("\n=== Figure S5 ===\n")
cs_sm <- read.csv("figs5_cs_monthly_summary.csv",
                  stringsAsFactors = FALSE)

# GATE: file Competing-risk rows must equal authoritative bootstrap medians
cr_chk <- cs_sm[cs_sm$method == "Competing-risk", ]
stopifnot(nrow(cr_chk) == 24)
gate_cr <- function(arm, cols) {
  k <- cr_chk$arm == arm
  stopifnot(max(abs(cr_chk$est[k] - bs[paste0(cols, 1:12), "point"])) < 1e-6,
            max(abs(cr_chk$lo[k]  - bs[paste0(cols, 1:12), "lo"]))    < 1e-6,
            max(abs(cr_chk$hi[k]  - bs[paste0(cols, 1:12), "hi"]))    < 1e-6)
}
gate_cr("NON_USER", "apNU_")
gate_cr("STATIN",   "apST_")

cif_compare <- cs_sm[, c("month", "arm", "method", "est", "lo", "hi")]
cif_compare$arm <- factor(cif_compare$arm, levels = c("STATIN", "NON_USER"))
# CR drawn first (bottom layer, thick solid), CS drawn second (top, thick dashed with gaps)
cif_compare$method <- factor(cif_compare$method,
  levels = c("Competing-risk", "Cause-specific"))

# Split data for explicit layer control
cr_cmp <- cif_compare[cif_compare$method == "Competing-risk", ]
cs_cmp <- cif_compare[cif_compare$method == "Cause-specific", ]

# Strategy: CR solid lw=1.6 (BOTTOM), CS dashed lw=2.8 dash="64" (6pt dash/4pt gap, TOP)
# CS dash segments stick out from CR solid sides; CS gaps reveal CR solid below
p_s2 <- ggplot() +
  # Subtle CI ribbons for CR only
  geom_ribbon(data = cr_cmp, aes(x = month, ymin = lo, ymax = hi, fill = arm),
              alpha = 0.08, color = NA, show.legend = FALSE) +
  # Layer 1: CR thick solid (lw=1.6, BOTTOM)
  geom_line(data = cr_cmp, aes(x = month, y = est, color = arm, linetype = method),
            linewidth = 1.6) +
  # Layer 2: CS thick dashed (lw=2.8, dash="64", TOP -- dashes stick out, gaps show solid)
  geom_line(data = cs_cmp, aes(x = month, y = est, color = arm, linetype = method),
            linewidth = 2.8) +
  scale_color_manual(values = c("STATIN" = col_statin, "NON_USER" = col_nonuser),
                     labels = c("Statin initiators", "Non-initiators")) +
  scale_linetype_manual(values = c("Competing-risk" = "solid", "Cause-specific" = "64"),
                        labels = c("Competing-risk", "Cause-specific")) +
  scale_fill_manual(values = c("STATIN" = col_statin, "NON_USER" = col_nonuser),
                    guide = "none") +
  scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, NA)) +
  labs(x = "Month", y = "Cumulative incidence\nper 100 000", color = "Arm", linetype = "Method") +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
    linetype = guide_legend(order = 2, override.aes = list(color = "black", linewidth = 1.2))) +
  theme_fig(14) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm")
  )

save_fig("Figure S5", p_s2, 8, 5)

# ============================================================================
# FIGURE S2: SMD Balance (adjusted only, grouped by domain)
# ============================================================================
cat("\n=== Figure S2 ===\n")
smd_adj <- smd_full[, c("variable", "smd_after", "group")]
smd_adj$group <- factor(smd_adj$group,
  levels = c("Demographics", "Lipids", "Vitals/Renal",
             "Healthcare utilisation", "Comorbidities", "Medications"))
smd_adj <- smd_adj[order(smd_adj$group, smd_adj$smd_after), ]
smd_adj$variable <- factor(smd_adj$variable, levels = smd_adj$variable)

p_s3 <- ggplot(smd_adj, aes(x = smd_after, y = variable, color = group)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0.1, color = "black", linetype = "dashed", linewidth = 0.6) +
  geom_point(size = 2.2) +
  scale_color_manual(values = group_colors, name = "Covariate domain") +
  labs(x = "Absolute standardized mean difference (adjusted)", y = "") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_bw_fig(10) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(color = "black", size = 9, face = "plain", angle = -90),
    strip.background = element_rect(color = "black", fill = "grey95", linewidth = 0.4),
    axis.text.y = element_text(color = "black", size = 8)
  )

save_fig("Figure S2", p_s3, 7, 10)

# ============================================================================
# FIGURE S3: Monthly Rate Difference
# ============================================================================
cat("\n=== Figure S3 ===\n")
rd_df <- data.frame(
  month = 1:12,
  est = bs[paste0("mrd_", 1:12), "point"],
  lo  = bs[paste0("mrd_", 1:12), "lo"],
  hi  = bs[paste0("mrd_", 1:12), "hi"]
)

p_s4 <- ggplot(rd_df, aes(x = month, y = est)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.7, stroke = 0.9,
                  color = col_statin, shape = 19, linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(x = "Month", y = "Rate difference per 100 000") +
  theme_fig(14)

save_fig("Figure S3", p_s4, 8, 5)

# ============================================================================
# FIGURE S4: Death CIF
# ============================================================================
cat("\n=== Figure S4 ===\n")
death_cif <- data.frame(
  arm = factor(rep(c("NON_USER", "STATIN"), each = 12),
               levels = c("STATIN", "NON_USER")),
  month = rep(1:12, 2),
  est = c(bs[paste0("dNU_", 1:12), "point"], bs[paste0("dST_", 1:12), "point"]),
  lo  = c(bs[paste0("dNU_", 1:12), "lo"],    bs[paste0("dST_", 1:12), "lo"]),
  hi  = c(bs[paste0("dNU_", 1:12), "hi"],    bs[paste0("dST_", 1:12), "hi"])
)

p_s5 <- ggplot(death_cif, aes(x = month, y = est, color = arm, fill = arm)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  scale_color_manual(values = c("STATIN" = col_statin, "NON_USER" = col_nonuser),
                     labels = c("Statin initiators", "Non-initiators")) +
  scale_fill_manual(values = c("STATIN" = col_statin, "NON_USER" = col_nonuser),
                    labels = c("Statin initiators", "Non-initiators")) +
  scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, NA)) +
  labs(x = "Month", y = "Cumulative incidence\nper 100 000", color = "", fill = "") +
  theme_fig(14) +
  theme(legend.position = "bottom")

# Figure S4 (all-cause mortality): numbers at risk REMOVED per author request (2026-09-09)
# Original: p_s5 <- add_nrisk(p_s5, nr_st, nr_nu) -- bottom numbers clutter the figure
save_fig("Figure S4", p_s5, 8, 5)

# ============================================================================
# FIGURE S6: Grace sensitivity (AUTHORITATIVE re-run: monthly medians + percentile
# CIs from grace_fresh_weightit.rds replicate matrices; weights re-estimated per rep)
# ============================================================================
cat("\n=== Figure S6 ===\n")
grace_panel <- function(B, lab) {   # B: replicate matrix with apST_*/apNU_* per-100k
  data.frame(
    month = 1:12,
    statin_mean = apply(B[, paste0("apST_", 1:12)], 2, median, na.rm = TRUE) / 1e5,
    statin_lci  = apply(B[, paste0("apST_", 1:12)], 2, quantile, 0.025, na.rm = TRUE) / 1e5,
    statin_uci  = apply(B[, paste0("apST_", 1:12)], 2, quantile, 0.975, na.rm = TRUE) / 1e5,
    nonuser_mean = apply(B[, paste0("apNU_", 1:12)], 2, median, na.rm = TRUE) / 1e5,
    nonuser_lci  = apply(B[, paste0("apNU_", 1:12)], 2, quantile, 0.025, na.rm = TRUE) / 1e5,
    nonuser_uci  = apply(B[, paste0("apNU_", 1:12)], 2, quantile, 0.975, na.rm = TRUE) / 1e5,
    grace = lab
  )
}
g90_monthly  <- grace_panel(grace_boot[["90"]]$boot_res, "90 days")
g180_monthly <- grace_panel(boot_main, "180 days")
g365_monthly <- grace_panel(grace_boot[["365"]]$boot_res, "365 days")
cif_grace <- rbind(g90_monthly, g180_monthly, g365_monthly)
cif_grace$grace <- factor(cif_grace$grace, levels = c("90 days", "180 days", "365 days"))

grace_colors <- c("90 days" = "#2E7D32", "180 days" = col_statin, "365 days" = "#D32F2F")

p_s6 <- ggplot(cif_grace, aes(x = month, group = grace)) +
  facet_wrap(~grace, ncol = 3) +
  geom_ribbon(aes(ymin = statin_lci, ymax = statin_uci, fill = grace), alpha = 0.15) +
  geom_ribbon(aes(ymin = nonuser_lci, ymax = nonuser_uci, fill = grace), alpha = 0.15) +
  geom_line(aes(y = statin_mean, color = grace, linetype = "Statin initiators"), linewidth = 1.1) +
  geom_line(aes(y = nonuser_mean, color = grace, linetype = "Non-initiators"), linewidth = 1.1) +
  scale_color_manual(values = grace_colors, name = "Grace period") +
  scale_fill_manual(values = grace_colors, name = "Grace period") +
  scale_linetype_manual(values = c("Statin initiators" = "solid", "Non-initiators" = "dashed"), name = "") +
  labs(x = "Months since landmark", y = "Cumulative incidence") +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  theme_bw_fig(12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(color = "black", size = 11, face = "plain")
  )

save_fig("Figure S6", p_s6, 14, 5)

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Figures saved (PDF + PNG, large font, black text) | fibrate figures: none (cancelled) ===\n")
