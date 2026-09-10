library(dplyr)
library(mice)

df_master <- readRDS("df_master.rds")
imp <- readRDS("imp_statin_htg_v4.rds")

cat("=== 1. df_master ===\n")
cat("Total:", nrow(df_master), "\n")
cat("time_days <= 0 or NA:", sum(df_master$time_days <= 0 | is.na(df_master$time_days)), "\n")
cat("time_days > 0:", sum(df_master$time_days > 0, na.rm = TRUE), "\n")

cat("\n=== 2. imp$data ===\n")
cat("Total:", nrow(imp$data), "\n")
cat("time_days <= 0 or NA:", sum(imp$data$time_days <= 0 | is.na(imp$data$time_days)), "\n")
cat("time_days > 0:", sum(imp$data$time_days > 0, na.rm = TRUE), "\n")

cat("\n=== 3. Who differs? ===\n")
m_ids <- df_master$person_id
i_ids <- imp$data$person_id
in_m_not_i <- setdiff(m_ids, i_ids)
in_i_not_m <- setdiff(i_ids, m_ids)
cat("In df_master but not in imp:", length(in_m_not_i), "\n")
cat("In imp but not in df_master:", length(in_i_not_m), "\n")

# df_master 中被 MICE 排除的 3037 人，time_days 情况
df_excl <- df_master[df_master$person_id %in% in_m_not_i, ]
cat("\nExcluded from MICE (n=", nrow(df_excl), "):\n", sep="")
cat("  time_days <= 0 or NA:", sum(df_excl$time_days <= 0 | is.na(df_excl$time_days)), "\n")
cat("  time_days > 0:", sum(df_excl$time_days > 0, na.rm = TRUE), "\n")
cat("  arm distribution:\n")
print(table(df_excl$arm))

# imp 中 time_days <= 0 的 32 人是谁？
imp_t0 <- imp$data[imp$data$time_days <= 0 | is.na(imp$data$time_days), ]
cat("\nIn imp but time_days <= 0 (n=", nrow(imp_t0), "):\n", sep="")
cat("  arm distribution:\n")
print(table(imp_t0$arm))
cat("  person_ids in df_master?\n")
cat("  ", sum(imp_t0$person_id %in% df_master$person_id), "\n")
cat("  Their time_days in df_master:\n")
df_match <- df_master[df_master$person_id %in% imp_t0$person_id, ]
print(summary(df_match$time_days))

# 最终决策：分析集应该用谁？
# 应该用 imp$data 中有 time_days > 0 的人（有效随访人群）
# 还是 imp$data 全部（含 MICE 插补了 time_days 但 time<=0 的人）
cat("\n=== DECISION ===\n")
cat("Survival analysis uses person-month file — only rows with time_days > 0 contribute.\n")
cat("But the baseline cohort (for Table 1, N reporting) should use imp$data with valid follow-up.\n")
cat("\nFinal analytical cohort should be: imp$data where time_days > 0\n")
cat("  N =", sum(imp$data$time_days > 0, na.rm = TRUE), "\n")
cat("  statin =", sum(imp$data$arm == "STATIN" & imp$data$time_days > 0, na.rm = TRUE), "\n")
cat("  non-user =", sum(imp$data$arm == "NON_USER" & imp$data$time_days > 0, na.rm = TRUE), "\n")

# 冻结这个队列
frozen_ids <- imp$data$person_id[imp$data$time_days > 0 & !is.na(imp$data$time_days)]
cat("\nFrozen cohort: person_id exported to frozen_cohort_ids.rds\n")
cat("  N_total =", length(frozen_ids), "\n")
cat("  N_statin =", sum(imp$data$arm[imp$data$time_days > 0 & !is.na(imp$data$time_days)] == "STATIN"), "\n")
cat("  N_nonuser =", sum(imp$data$arm[imp$data$time_days > 0 & !is.na(imp$data$time_days)] == "NON_USER"), "\n")

frozen <- data.frame(
  person_id = frozen_ids,
  arm = imp$data$arm[imp$data$time_days > 0 & !is.na(imp$data$time_days)],
  time_days = imp$data$time_days[imp$data$time_days > 0 & !is.na(imp$data$time_days)]
)
saveRDS(frozen, "frozen_cohort_ids.rds")
write.csv(frozen, "frozen_cohort_ids.csv", row.names = FALSE)
