library(dplyr)
library(mice)
library(bit64)

imp <- readRDS("imp_statin_htg_v4.rds")

frozen_ids <- imp$data$person_id[imp$data$time_days > 0 & !is.na(imp$data$time_days)]
frozen <- data.frame(
  person_id = frozen_ids,
  arm = imp$data$arm[imp$data$time_days > 0 & !is.na(imp$data$time_days)],
  time_days = imp$data$time_days[imp$data$time_days > 0 & !is.na(imp$data$time_days)],
  stringsAsFactors = FALSE
)

cat("Frozen person_id class after data.frame:", class(frozen$person_id), "\n")
cat("head:", head(frozen$person_id), "\n")

# 用 tibble 或直接存 list 来保持 integer64
frozen_list <- list(
  person_id = frozen_ids,
  arm = imp$data$arm[imp$data$time_days > 0 & !is.na(imp$data$time_days)],
  time_days = imp$data$time_days[imp$data$time_days > 0 & !is.na(imp$data$time_days)]
)
class(frozen_list) <- "data.frame"
attr(frozen_list, "row.names") <- seq_along(frozen_ids)

cat("\nFinal frozen person_id class:", class(frozen_list$person_id), "\n")
cat("N total:", length(frozen_list$person_id), "\n")
cat("N statin:", sum(frozen_list$arm == "STATIN"), "\n")
cat("N non-user:", sum(frozen_list$arm == "NON_USER"), "\n")

saveRDS(frozen_list, "frozen_cohort_ids.rds")
cat("\nSaved to frozen_cohort_ids.rds (person_id preserved as integer64)\n")

# 验证：读出来能和 imp$data 匹配
frozen_read <- readRDS("frozen_cohort_ids.rds")
cat("\nVerification:\n")
cat("read back class:", class(frozen_read$person_id), "\n")
cat("intersect with imp$data:", length(intersect(frozen_read$person_id, imp$data$person_id)), "\n")
cat("imp data nrow:", nrow(imp$data), "\n")
