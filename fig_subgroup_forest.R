# ============================================================================
# Subgroup forest plot: Age (<65 vs >=65) and Sex (Male vs Female)
# Data: subgroup_results.csv (5 MICE x 200 bootstrap per subgroup, seed 20260625)
#       + overall RR from bootstrap_fresh_weightit.rds (main analysis)
# Output: Figure subgroup_forest.png/.pdf
# ============================================================================
library(ggplot2); library(data.table)

# --- read subgroup results (script-generated, no manual numbers) ---
sg <- fread("subgroup_results.csv")

# --- read overall RR from authoritative main-analysis RDS ---
boot <- readRDS("bootstrap_fresh_weightit.rds")
rr_ap <- boot[, "rr_ap"]; rr_ap <- rr_ap[is.finite(rr_ap) & rr_ap > 0 & rr_ap < 100]
rr_d  <- boot[, "rr_death"]; rr_d <- rr_d[is.finite(rr_d) & rr_d > 0 & rr_d < 100]
overall <- data.frame(
  subgroup = "Overall", outcome = c("AP", "Death"),
  rr = c(median(rr_ap), median(rr_d)),
  ci_lo = c(quantile(rr_ap, 0.025), quantile(rr_d, 0.025)),
  ci_hi = c(quantile(rr_ap, 0.975), quantile(rr_d, 0.975)),
  events_statin = NA_integer_, events_nonuser = NA_integer_,
  epv = NA_real_, interaction_p = NA_real_
)

# --- combine ---
d <- rbindlist(list(as.data.table(overall), sg))
d[, label := fcase(
  subgroup == "Overall", "Overall",
  subgroup == "Age<65",  "Age <65 yr",
  subgroup == "Age>=65", "Age \u226565 yr",
  subgroup == "Male",    "Male",
  subgroup == "Female",  "Female"
)]
d[, outcome_lab := fcase(outcome == "AP", "Acute pancreatitis",
                          outcome == "Death", "All-cause mortality")]

# row position within each panel (1=bottom Overall ... 5=top Female)
d[, row := fcase(subgroup == "Overall", 1,
                 subgroup == "Age<65",  2,
                 subgroup == "Age>=65", 3,
                 subgroup == "Male",    4,
                 subgroup == "Female",  5)]

# interaction P per panel (same value for both rows of a factor)
p_age <- unique(sg[outcome == "AP" & subgroup %in% c("Age<65","Age>=65"), interaction_p])
p_sex <- unique(sg[outcome == "AP" & subgroup %in% c("Male","Female"), interaction_p])
d[, p_int := fcase(grepl("Age", subgroup), sprintf("P(int)=%.2f", p_age),
                   grepl("Male|Female", subgroup), sprintf("P(int)=%.2f", p_sex),
                   default = NA_character_)]

# events display
d[, ev := ifelse(is.na(events_statin), "",
                 sprintf("%d/%d", events_statin, events_nonuser))]

# RR (95% CI) text; near-zero lower bound shown as <0.01
d[, rr_txt := sprintf("%.2f (%s-%.2f)", rr,
                      ifelse(ci_lo < 0.01, "<0.01", sprintf("%.2f", ci_lo)), ci_hi)]

# EPV flag for sparse subgroups
d[, sparse := !is.na(epv) & epv < 10]

# --- plot ---
# interaction P annotations: top-left corner of each panel (clear of event-count column)
grp_ann <- data.frame(
  outcome_lab = rep(c("Acute pancreatitis", "All-cause mortality"), each = 2),
  y = rep(c(6.35, 5.7), 2),
  lab = rep(c(sprintf("Age, P for interaction = %.2f", p_age),
              sprintf("Sex, P for interaction = %.2f", p_sex)), 2)
)

p <- ggplot(d, aes(x = rr, y = row)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_pointrange(aes(xmin = pmax(pmin(ci_lo, 8), 0.02), xmax = pmin(ci_hi, 8)),
                  size = 0.55, linewidth = 0.55,
                  colour = "#1565C0", shape = 16) +
  geom_text(aes(x = 0.045, label = rr_txt),
            hjust = 0, vjust = -1.1, size = 2.6, colour = "black") +
  geom_text(aes(x = 8, label = ev), hjust = 1, size = 2.6, colour = "black") +
  geom_text(data = grp_ann, aes(x = 0.022, y = y, label = lab),
            hjust = 0, size = 2.6, colour = "grey20", fontface = "italic",
            inherit.aes = FALSE) +
  facet_wrap(~ outcome_lab, ncol = 2) +
  scale_x_continuous(trans = "log10",
                     breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8),
                     labels = c("0.05", "0.1", "0.25", "0.5", "1", "2", "4", "8"),
                     limits = c(0.02, 10)) +
  scale_y_continuous(breaks = 1:5,
                     labels = c("Overall", "Age <65 yr", "Age \u226565 yr",
                                "Male", "Female"),
                     limits = c(0.4, 6.55)) +
  labs(x = "12-month cumulative-incidence RR (log scale)",
       y = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(colour = "black"),
        plot.margin = margin(5, 5, 5, 5))

ggsave("Figure_subgroup_forest.png", p,
       width = 190, height = 100, units = "mm", dpi = 600)
ggsave("Figure_subgroup_forest.pdf", p,
       width = 190, height = 100, units = "mm")

# console check values used in plot
print(d[, .(outcome_lab, label, rr, ci_lo, ci_hi, ev, p_int)])
