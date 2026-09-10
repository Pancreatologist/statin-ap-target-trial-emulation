# =====================================================================
# 第 1-3 步：MICE 插补  ->  熵平衡(ATT)  ->  平衡诊断
# 数据：df_final_statin_htg_v4.rds（更好的 ICD9 匹配）
# 估计量：ATT(他汀使用者中的效应)
# =====================================================================

# ---- 依赖 ----
need <- c("mice","MatchThem","WeightIt","cobalt","dplyr")
for (p in need) if (!requireNamespace(p, quietly=TRUE))
  install.packages(p, repos="https://cloud.r-project.org")
library(mice); library(MatchThem); library(WeightIt); library(cobalt); library(dplyr)

# ---- 读取数据 ----
df_final <- readRDS("df_final_statin_htg_v4.rds")

# ---- 排除无随访时间的人 ----
df_final <- df_final %>% filter(!is.na(time_days))

# =====================================================================
# 0) 变量准备
# =====================================================================
# 连续变量 1%-99% 缩尾（处理极端异常值）
winsorize <- function(x, probs = c(0.01, 0.99)) {
  q <- quantile(x, probs, na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

cont_to_winsorize <- c("ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
                       "n_inpatient_ed","n_outpatient","n_visits_total")
for (v in cont_to_winsorize) {
  if (v %in% names(df_final)) {
    df_final[[v]] <- winsorize(df_final[[v]])
  }
}

# 人口学变量简化
df_final <- df_final %>%
  mutate(
    hispanic   = as.integer(ethnicity == "Hispanic or Latino"),
    smoke_ever = as.integer(smoke_status %in% c("former","current")),
    male       = as.integer(sex == "Male"),
    white      = as.integer(race == "White")
  )

# 协变量清单(PS 模型用)
demo_vars  <- c("age","male","white","hispanic","smoke_ever")
cont_vars  <- c("ldl","hdl","tc","tg","bmi","calcium","alt","egfr",
                "n_inpatient_ed","n_outpatient","n_visits_total")
comorb_cols<- c("htn","obesity","diabetes","chd","stroke","pvd","afib","kidney",
                "dementia","as_spond","ra","sle","pulmonary","cancer","neuro")
drug_cols  <- c("aspirin","insulin","metformin","su","tzd","glinide","agi","dpp4","glp1",
                "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn",
                "antithromb","nsaid")
ps_vars <- c(demo_vars, "alcohol", cont_vars, comorb_cols, drug_cols)

# 建模子集(含结局做插补辅助变量)
dat <- df_final %>%
  select(person_id, arm, ap_event, death_event, time_days, all_of(ps_vars)) %>%
  mutate(
    arm        = factor(arm, levels=c("NON_USER","STATIN")),
    across(all_of(c(comorb_cols, drug_cols, demo_vars[-1])), ~as.integer(.)),
    alcohol    = factor(alcohol, ordered=TRUE)
  )

cat('建模子集维度:', nrow(dat), 'x', ncol(dat), '\n')
cat('arm分布:\n')
print(table(dat$arm))
cat('AP事件数 by arm:\n')
print(table(dat$arm, dat$ap_event))
cat('\n需插补变量缺失数:\n')
print(colSums(is.na(dat[, c("ldl","hdl","tc","tg","bmi","calcium","alt","egfr","alcohol")])))

# 检查哪些二分类变量方差为0
cat('\n=== 方差检查（二分类变量）===\n')
binary_vars <- c(comorb_cols, drug_cols, demo_vars[-1])
zero_var <- c()
for (v in binary_vars) {
  if (v %in% names(dat)) {
    n1 <- sum(dat[[v]] == 1, na.rm = TRUE)
    n0 <- sum(dat[[v]] == 0, na.rm = TRUE)
    if (n1 == 0 || n0 == 0) {
      cat(sprintf('  零方差: %s (n1=%d, n0=%d)\n', v, n1, n0))
      zero_var <- c(zero_var, v)
    }
  }
}
cat('方差检查完成\n')
if (length(zero_var) > 0) {
  cat('剔除零方差变量:', zero_var, '\n')
  comorb_cols <- setdiff(comorb_cols, zero_var)
  drug_cols <- setdiff(drug_cols, zero_var)
  demo_vars <- setdiff(demo_vars, zero_var)
  ps_vars <- setdiff(ps_vars, zero_var)
}

# =====================================================================
# 1) MICE 多重插补
# =====================================================================
imp_file <- "imp_statin_htg_v4.rds"

if (file.exists(imp_file)) {
  cat('\n>> 加载已保存的 MICE 结果...\n')
  imp <- readRDS(imp_file)
  cat('MICE 结果已加载\n')
} else {
  ini <- mice(dat, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix

  pred[, "person_id"] <- 0
  pred["person_id", ] <- 0
  meth["person_id"]   <- ""

  complete_cat <- c("arm","ap_event","death_event","time_days",
                    comorb_cols, drug_cols, demo_vars[-1])
  meth[complete_cat] <- ""

  for (v in c("ldl","hdl","tc","tg","bmi","calcium","alt","egfr")) meth[v] <- "pmm"
  meth["alcohol"] <- "polr"

  set.seed(20260625)
  cat('\n>> 运行 MICE(m=5)... 大样本可能需数分钟\n')
  imp <- mice(dat, m = 5, maxit = 5, method = meth, predictorMatrix = pred,
              printFlag = FALSE)
  cat('MICE 完成。logged events(若有):\n'); print(head(imp$loggedEvents))
  
  saveRDS(imp, imp_file)
  cat('已保存 MICE 结果到', imp_file, '\n')
}

# =====================================================================
# 2) 熵平衡(ATT)
# =====================================================================
ps_formula <- as.formula(paste("arm ~", paste(ps_vars, collapse = " + ")))

cat('\n>> weightthem: 熵平衡(ATT)...\n')
w.att <- weightthem(
  ps_formula,
  datasets  = imp,
  approach  = "within",
  method    = "ebal",
  estimand  = "ATT"
)

# =====================================================================
# 3) 平衡诊断
# =====================================================================
bt <- bal.tab(w.att, un = TRUE, thresholds = c(m = 0.1),
              stats = "mean.diffs", abs = TRUE)
cat('\n========== 平衡诊断(跨插补平均 SMD) ==========\n'); print(bt)

bal <- bt$Balance
adj_col <- grep("Mean.Diff.Adj", colnames(bal), value = TRUE)

if (length(adj_col) > 0) {
  smd_adj <- bal[[adj_col[1]]]
  over <- bal[!is.na(smd_adj) & smd_adj > 0.1, , drop = FALSE]
  cat('\n--- 加权后仍 SMD>0.1 的变量 ---\n')
  if (nrow(over) == 0) cat('全部 ≤0.1，平衡达标 ✓\n') else print(over)
}

# Love plot
lp <- love.plot(w.att, stats = "mean.diffs", abs = TRUE,
                thresholds = c(m = 0.1), drop.distance = TRUE,
                var.order = "unadjusted", line = TRUE,
                title = "Covariate Balance: Statin vs Non-user (ATT, Entropy Balancing, v4)")
print(lp)
ggplot2::ggsave("loveplot_statin_htg_v4.png", lp, width = 8, height = 12, dpi = 150)

# 权重分布
ws <- unlist(lapply(complete(w.att, "all"), function(d) d$weights))
cat('\n--- 熵平衡权重分布 ---\n'); print(summary(ws))
cat('权重 >10 的比例:', round(mean(ws > 10), 4), '\n')

saveRDS(w.att, "watt_statin_htg_v4.rds")
saveRDS(bt,    "baltab_statin_htg_v4.rds")
cat('\n已保存 watt_statin_htg_v4.rds / baltab_statin_htg_v4.rds / loveplot_statin_htg_v4.png\n')
