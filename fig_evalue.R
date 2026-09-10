# ============================================================================
# Figure: E-value plot for unmeasured confounding robustness
# Shows E-value (point) and E-value (CI bound nearest null) for AP and Death
# ============================================================================
library(ggplot2); library(scales)

# E-value data (computed with the VanderWeele-Ding method).
# Point estimates and CI-bound E-values for each outcome are study results;
# they are left as placeholders here so that this figure script does not
# embed unpublished findings. Replace with your own computed values.
evalue_df <- data.frame(
  outcome = c("Acute pancreatitis\n(primary)", "All-cause mortality\n(secondary)"),
  # Placeholders: replace with your own computed E-values (point estimate and
  # the CI bound nearest 1) so this figure script does not embed study results.
  e_point = c(NA_real_, NA_real_),
  e_ci    = c(NA_real_, NA_real_),
  stringsAsFactors = FALSE
)
evalue_df$outcome <- factor(evalue_df$outcome, levels = rev(evalue_df$outcome))

col_ap    <- "#1565C0"   # statin blue
col_death <- "#EF6C00"   # non-user orange (reused for death to distinguish)
point_colors <- c("Acute pancreatitis\n(primary)" = col_ap,
                  "All-cause mortality\n(secondary)" = col_death)

p_evalue <- ggplot(evalue_df, aes(x = e_point, y = outcome)) +
  geom_vline(xintercept = 1, color = "grey50", linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 1.5, color = "grey70", linetype = "dotted", linewidth = 0.5) +
  geom_segment(aes(x = e_ci, xend = e_point, y = outcome, yend = outcome),
               linewidth = 1.2, color = "grey40") +
  geom_point(aes(x = e_point, color = outcome), shape = 16, size = 4) +
  geom_point(aes(x = e_ci, color = outcome), shape = 1, size = 3.5, stroke = 1.2) +
  geom_text(aes(x = e_point, label = sprintf("%.2f", e_point)),
            hjust = -0.5, vjust = -0.8, size = 4, color = "black", fontface = "bold") +
  geom_text(aes(x = e_ci, label = sprintf("%.2f", e_ci)),
            hjust = -0.5, vjust = 1.5, size = 3.5, color = "grey40") +
  scale_color_manual(values = point_colors, guide = "none") +
  scale_x_continuous(limits = c(0.8, 2.5), breaks = c(1, 1.5, 2, 2.5)) +
  labs(x = "E-value", y = "") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14)
  )

ggsave("Figure E-value.pdf", p_evalue, width = 8, height = 4)
ggsave("Figure E-value.png", p_evalue, width = 8, height = 4, dpi = 600)
cat("Saved: Figure E-value.pdf/.png\n")
