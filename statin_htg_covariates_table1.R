# =====================================================================
# 协变量提取定稿 + 基线 Table 1
# 研究：HTG 人群，他汀 vs non-user，结局急性胰腺炎(landmark, 源B, grace=180)
# 产出：person-level 宽表 df  +  分臂 Table 1
# 已并入三处清理：饮酒锁 concept_id=40771103 / LDL 去摩尔单位 / 吸烟双题合成
# 平台：All of Us Researcher Workbench(R)
# =====================================================================

library(bigrquery); library(glue); library(dplyr); library(tidyr)
cdr     <- Sys.getenv("WORKSPACE_CDR")
billing <- Sys.getenv("GOOGLE_CLOUD_PROJECT")
stopifnot(nchar(cdr) > 0, nchar(billing) > 0)
run_bq <- function(sql) bq_table_download(bq_project_query(billing, sql), bigint = "integer64")

GRACE_DAYS <- 180

# ---------------- 队列 CTE(内联) ----------------
COHORT_CTE <- glue("
statin_ing AS (SELECT concept_id FROM `{cdr}.concept`
  WHERE vocabulary_id='RxNorm' AND concept_class_id='Ingredient'
    AND LOWER(concept_name) IN ('atorvastatin','simvastatin','rosuvastatin','pravastatin',
                                'lovastatin','pitavastatin','fluvastatin')),
statin_drugs AS (SELECT DISTINCT ca.descendant_concept_id concept_id
  FROM `{cdr}.concept_ancestor` ca JOIN statin_ing i ON ca.ancestor_concept_id=i.concept_id),
htg_concepts AS (SELECT DISTINCT cr.concept_id_2 concept_id
  FROM `{cdr}.concept` c JOIN `{cdr}.concept_relationship` cr
    ON c.concept_id=cr.concept_id_1 AND cr.relationship_id='Maps to'
  WHERE (c.vocabulary_id='ICD10CM' AND (c.concept_code LIKE 'E78.1%' OR c.concept_code LIKE 'E78.2%'))
     OR (c.vocabulary_id='ICD9CM'  AND (c.concept_code='272.1' OR c.concept_code='272.4'))),
tg_concepts AS (SELECT concept_id FROM `{cdr}.concept`
  WHERE domain_id='Measurement' AND standard_concept='S' AND LOWER(concept_name) LIKE '%triglyceride%'),
panc_excl_concepts AS (SELECT DISTINCT cr.concept_id_2 concept_id
  FROM `{cdr}.concept` c JOIN `{cdr}.concept_relationship` cr
    ON c.concept_id=cr.concept_id_1 AND cr.relationship_id='Maps to'
  WHERE (c.vocabulary_id='ICD10CM' AND (c.concept_code LIKE 'K85%' OR c.concept_code LIKE 'K86.0%'
          OR c.concept_code LIKE 'K86.1%' OR c.concept_code LIKE 'C25%'))
     OR (c.vocabulary_id='ICD9CM' AND (c.concept_code LIKE '577%' OR c.concept_code LIKE '157%'))),
ap_concepts AS (SELECT DISTINCT cr.concept_id_2 concept_id
  FROM `{cdr}.concept` c JOIN `{cdr}.concept_relationship` cr
    ON c.concept_id=cr.concept_id_1 AND cr.relationship_id='Maps to'
  WHERE (c.vocabulary_id='ICD10CM' AND c.concept_code LIKE 'K85%')
     OR (c.vocabulary_id='ICD9CM' AND c.concept_code LIKE '577.0%')),
htg_by_code AS (SELECT co.person_id, MIN(co.condition_start_date) htg_date
  FROM `{cdr}.condition_occurrence` co JOIN htg_concepts h ON co.condition_concept_id=h.concept_id GROUP BY 1),
htg_by_lab AS (SELECT m.person_id, MIN(m.measurement_date) htg_date
  FROM `{cdr}.measurement` m JOIN tg_concepts t ON m.measurement_concept_id=t.concept_id
  WHERE m.value_as_number >= 500 GROUP BY 1),
htg_B AS (SELECT person_id, MIN(htg_date) htg_date FROM (
    SELECT * FROM htg_by_code UNION ALL SELECT * FROM htg_by_lab) GROUP BY 1),
statin_first AS (SELECT de.person_id, MIN(de.drug_exposure_start_date) first_statin
  FROM `{cdr}.drug_exposure` de JOIN statin_drugs d ON de.drug_concept_id=d.concept_id GROUP BY 1),
panc_first AS (SELECT co.person_id, MIN(co.condition_start_date) first_panc
  FROM `{cdr}.condition_occurrence` co JOIN panc_excl_concepts pe ON co.condition_concept_id=pe.concept_id GROUP BY 1),
ap_qualifying AS (SELECT co.person_id, MIN(co.condition_start_date) ap_date
  FROM `{cdr}.condition_occurrence` co JOIN ap_concepts a ON co.condition_concept_id=a.concept_id
  JOIN `{cdr}.visit_occurrence` v ON co.visit_occurrence_id=v.visit_occurrence_id
  WHERE v.visit_concept_id IN (9201,9203,262) GROUP BY 1),
eligible AS (SELECT hp.person_id, hp.htg_date index_date,
         DATE_ADD(hp.htg_date, INTERVAL {GRACE_DAYS} DAY) landmark, sf.first_statin
  FROM htg_B hp JOIN `{cdr}.observation_period` op ON hp.person_id=op.person_id
  LEFT JOIN statin_first sf ON hp.person_id=sf.person_id
  LEFT JOIN panc_first  pf ON hp.person_id=pf.person_id
  WHERE DATE_DIFF(hp.htg_date, op.observation_period_start_date, DAY) >= 365
    AND (sf.first_statin IS NULL OR sf.first_statin >= hp.htg_date)
    AND (pf.first_panc  IS NULL OR pf.first_panc  >= hp.htg_date)),
assigned AS (SELECT e.person_id, e.index_date, e.landmark AS t0,
         CASE WHEN e.first_statin IS NOT NULL AND e.first_statin < e.landmark
                   AND e.first_statin >= e.index_date THEN 'STATIN' ELSE 'NON_USER' END arm
  FROM eligible e LEFT JOIN ap_qualifying aq ON e.person_id=aq.person_id
  WHERE aq.ap_date IS NULL OR aq.ap_date >= e.landmark),
cohort AS (SELECT person_id, index_date, t0, arm FROM assigned)
")

# =====================================================================
# 1) 队列 + 人口学  -> df_base
# =====================================================================
df_base <- run_bq(glue("WITH {COHORT_CTE}
SELECT c.person_id, c.index_date, c.t0, c.arm,
  DATE_DIFF(c.t0, DATE(p.birth_datetime), DAY)/365.25 AS age,
  gc.concept_name  AS sex,
  rc.concept_name  AS race,
  ec.concept_name  AS ethnicity
FROM cohort c
JOIN `{cdr}.person` p USING(person_id)
LEFT JOIN `{cdr}.concept` gc ON p.gender_concept_id=gc.concept_id
LEFT JOIN `{cdr}.concept` rc ON p.race_concept_id=rc.concept_id
LEFT JOIN `{cdr}.concept` ec ON p.ethnicity_concept_id=ec.concept_id"))
cat('df_base:', nrow(df_base), '行\n')

# =====================================================================
# 2) 血脂/检验  -> df_lab (t0 前最近一次值；LDL 仅 mg/dL，去摩尔单位)
# =====================================================================
lab_map_sql <- "
lab_map AS (
  SELECT concept_id, 'ldl' label FROM UNNEST([3028288,3009966,3028437]) concept_id  -- 13457-7/18262-6/2089-1 (mg/dL)
  UNION ALL SELECT concept_id,'hdl' FROM UNNEST([3007070]) concept_id               -- 2085-9
  UNION ALL SELECT concept_id,'tc'  FROM UNNEST([3027114]) concept_id               -- 2093-3
  UNION ALL SELECT concept_id,'tg'  FROM `{cdr}.concept`
            WHERE domain_id='Measurement' AND standard_concept='S' AND LOWER(concept_name) LIKE '%triglyceride%'
  UNION ALL SELECT concept_id,'bmi' FROM `{cdr}.concept`
            WHERE domain_id='Measurement' AND standard_concept='S' AND LOWER(concept_name) LIKE '%body mass index%'
  UNION ALL SELECT concept_id,'calcium' FROM `{cdr}.concept`
            WHERE domain_id='Measurement' AND standard_concept='S' AND LOWER(concept_name) LIKE '%calcium%'
  UNION ALL SELECT concept_id,'alt' FROM `{cdr}.concept`
            WHERE domain_id='Measurement' AND standard_concept='S' AND LOWER(concept_name) LIKE '%alanine aminotransferase%'
  UNION ALL SELECT concept_id,'egfr' FROM `{cdr}.concept`
            WHERE domain_id='Measurement' AND standard_concept='S' AND LOWER(concept_name) LIKE '%glomerular filtration%'
)"
df_lab <- run_bq(glue("WITH {COHORT_CTE}, ", glue(lab_map_sql), ",
lab_vals AS (
  SELECT cb.person_id, lm.label, mm.value_as_number val,
         ROW_NUMBER() OVER (PARTITION BY cb.person_id, lm.label ORDER BY mm.measurement_date DESC) rn
  FROM cohort cb
  JOIN `{cdr}.measurement` mm ON cb.person_id=mm.person_id
       AND mm.measurement_date < cb.t0 AND mm.value_as_number IS NOT NULL
  JOIN lab_map lm ON mm.measurement_concept_id=lm.concept_id)
SELECT person_id,
  MAX(IF(label='ldl', val, NULL))     AS ldl,
  MAX(IF(label='hdl', val, NULL))     AS hdl,
  MAX(IF(label='tc',  val, NULL))     AS tc,
  MAX(IF(label='tg',  val, NULL))     AS tg,
  MAX(IF(label='bmi', val, NULL))     AS bmi,
  MAX(IF(label='calcium', val, NULL)) AS calcium,
  MAX(IF(label='alt', val, NULL))     AS alt,
  MAX(IF(label='egfr', val, NULL))    AS egfr
FROM lab_vals WHERE rn=1 GROUP BY person_id"))
cat('df_lab:', nrow(df_lab), '行\n')

# =====================================================================
# 3) 基线合并症 flags  -> df_comorbid
# =====================================================================
comorb_defs <- list(
  htn        = c("I10%"),
  obesity    = c("E66%"),
  diabetes   = c("E11%"),
  chd        = c("I20%","I21%","I22%","I23%","I24%","I25%"),
  stroke     = c("I60%","I61%","I62%","I63%","I64%","I69%"),
  pvd        = c("I70%","I73%"),
  afib       = c("I48%"),
  kidney     = c("N17%","N18%","N19%"),
  dementia   = c("F00%","F01%","F02%","F03%","G30%"),
  as_spond   = c("M45%"),
  ra         = c("M05%","M06%"),
  sle        = c("M32%"),
  pulmonary  = c("J40%","J41%","J42%","J43%","J44%","J45%","J47%"),
  cancer     = c("C%"),
  neuro      = c("G%")
)
make_cond_map <- function(defs){
  paste(vapply(names(defs), function(lab){
    pats <- paste0("c.concept_code LIKE '", defs[[lab]], "'", collapse=" OR ")
    glue("SELECT DISTINCT cr.concept_id_2 concept_id, '{lab}' label
          FROM `{cdr}.concept` c JOIN `{cdr}.concept_relationship` cr
          ON c.concept_id=cr.concept_id_1 AND cr.relationship_id='Maps to'
          WHERE c.vocabulary_id='ICD10CM' AND ({pats})")
  }, character(1)), collapse=" UNION ALL ")
}
comorb_cols <- paste0("MAX(IF(m.label='", names(comorb_defs), "',1,0)) AS ", names(comorb_defs), collapse=",\n  ")
df_comorbid <- run_bq(glue("WITH {COHORT_CTE},
cond_map AS ({make_cond_map(comorb_defs)})
SELECT cb.person_id,\n  {comorb_cols}
FROM cohort cb
LEFT JOIN `{cdr}.condition_occurrence` co ON cb.person_id=co.person_id AND co.condition_start_date < cb.t0
LEFT JOIN cond_map m ON co.condition_concept_id=m.concept_id
GROUP BY cb.person_id"))
cat('df_comorbid:', nrow(df_comorbid), '行\n')

# =====================================================================
# 4) 基线前药物 flags  -> df_drug
# =====================================================================
drug_defs <- list(
  aspirin   = c("aspirin"),
  metformin = c("metformin"),
  su        = c("glipizide","glimepiride","glyburide","gliclazide","chlorpropamide","tolbutamide","tolazamide"),
  tzd       = c("pioglitazone","rosiglitazone"),
  glinide   = c("repaglinide","nateglinide"),
  agi       = c("acarbose","miglitol"),
  dpp4      = c("sitagliptin","saxagliptin","linagliptin","alogliptin"),
  glp1      = c("dulaglutide","exenatide","liraglutide","lixisenatide","semaglutide","albiglutide","tirzepatide"),
  sglt2     = c("canagliflozin","dapagliflozin","empagliflozin","ertugliflozin"),
  acei_arb  = c("lisinopril","enalapril","ramipril","benazepril","captopril","quinapril","perindopril",
                "fosinopril","trandolapril","moexipril","losartan","valsartan","olmesartan","irbesartan",
                "candesartan","telmisartan","azilsartan","eprosartan"),
  bblocker  = c("metoprolol","atenolol","carvedilol","bisoprolol","propranolol","nebivolol","labetalol",
                "nadolol","sotalol","timolol","acebutolol","betaxolol"),
  ccb       = c("amlodipine","diltiazem","verapamil","nifedipine","felodipine","nicardipine","isradipine","nisoldipine"),
  diuretic  = c("hydrochlorothiazide","chlorthalidone","furosemide","spironolactone","bumetanide","torsemide",
                "metolazone","indapamide","triamterene","amiloride","eplerenone","chlorothiazide"),
  alpha_blk = c("doxazosin","prazosin","terazosin"),
  other_htn = c("hydralazine","clonidine","minoxidil","methyldopa","guanfacine"),
  antithromb= c("clopidogrel","prasugrel","ticagrelor","warfarin","apixaban","rivaroxaban","dabigatran",
                "edoxaban","dipyridamole","cilostazol","enoxaparin","heparin","fondaparinux"),
  nsaid     = c("ibuprofen","naproxen","diclofenac","celecoxib","meloxicam","indomethacin","ketorolac",
                "etodolac","nabumetone","piroxicam","sulindac","ketoprofen","oxaprozin","diflunisal")
)
make_drug_map <- function(defs){
  blocks <- vapply(names(defs), function(lab){
    nm <- paste0("'", defs[[lab]], "'", collapse=",")
    glue("SELECT DISTINCT ca.descendant_concept_id concept_id, '{lab}' label
          FROM `{cdr}.concept` ing JOIN `{cdr}.concept_ancestor` ca ON ca.ancestor_concept_id=ing.concept_id
          WHERE ing.vocabulary_id='RxNorm' AND ing.concept_class_id='Ingredient'
            AND LOWER(ing.concept_name) IN ({nm})")
  }, character(1))
  # 胰岛素特殊：成分名含 insulin
  ins <- glue("SELECT DISTINCT ca.descendant_concept_id concept_id, 'insulin' label
               FROM `{cdr}.concept` ing JOIN `{cdr}.concept_ancestor` ca ON ca.ancestor_concept_id=ing.concept_id
               WHERE ing.vocabulary_id='RxNorm' AND ing.concept_class_id='Ingredient'
                 AND LOWER(ing.concept_name) LIKE '%insulin%'")
  paste(c(blocks, ins), collapse=" UNION ALL ")
}
drug_labels <- c(names(drug_defs), "insulin")
drug_cols <- paste0("MAX(IF(dm.label='", drug_labels, "',1,0)) AS ", drug_labels, collapse=",\n  ")
df_drug <- run_bq(glue("WITH {COHORT_CTE},
drug_map AS ({make_drug_map(drug_defs)})
SELECT cb.person_id,\n  {drug_cols}
FROM cohort cb
LEFT JOIN `{cdr}.drug_exposure` de ON cb.person_id=de.person_id AND de.drug_exposure_start_date < cb.t0
LEFT JOIN drug_map dm ON de.drug_concept_id=dm.concept_id
GROUP BY cb.person_id"))
cat('df_drug:', nrow(df_drug), '行\n')

# =====================================================================
# 5) survey 饮酒/吸烟  -> df_survey
#    饮酒：锁 concept_id=40771103 + value 前缀 DrinkFrequencyPastYear_ (去 COPE)
#    吸烟：100CigsLifetime_(ever) + Smoke_(当前) 合成 never/former/current
# =====================================================================
df_survey <- run_bq(glue("WITH {COHORT_CTE},
alc AS (
  SELECT person_id, value_source_value AS alc_raw FROM `{cdr}.observation`
  WHERE observation_concept_id = 40771103
    AND value_source_value LIKE 'DrinkFrequencyPastYear_%'
  QUALIFY ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY observation_datetime DESC)=1),
ever AS (
  SELECT DISTINCT person_id, TRUE AS ever_smoke FROM `{cdr}.observation`
  WHERE value_source_value='100CigsLifetime_Yes'),
curr AS (
  SELECT DISTINCT person_id, TRUE AS curr_smoke FROM `{cdr}.observation`
  WHERE value_source_value IN ('Smoke_EveryDay','Smoke_SomeDays'))
SELECT cb.person_id, a.alc_raw, e.ever_smoke, c2.curr_smoke
FROM cohort cb
LEFT JOIN alc  a ON cb.person_id=a.person_id
LEFT JOIN ever e ON cb.person_id=e.person_id
LEFT JOIN curr c2 ON cb.person_id=c2.person_id"))
cat('df_survey:', nrow(df_survey), '行\n')

# =====================================================================
# 6) CCI (Charlson) —— 下载 t0 前 ICD10 源码，用 comorbidity 包计算
# =====================================================================
cci_codes <- run_bq(glue("WITH {COHORT_CTE}
SELECT DISTINCT cb.person_id, src.concept_code AS icd10
FROM cohort cb
JOIN `{cdr}.condition_occurrence` co ON cb.person_id=co.person_id AND co.condition_start_date < cb.t0
JOIN `{cdr}.concept` src ON co.condition_source_concept_id=src.concept_id AND src.vocabulary_id='ICD10CM'"))
cat('cci_codes:', nrow(cci_codes), '行(person×ICD10)\n')

if (!requireNamespace("comorbidity", quietly=TRUE))
  install.packages("comorbidity", repos="https://cloud.r-project.org")
library(comorbidity)
cm <- comorbidity(cci_codes, id="person_id", code="icd10",
                  map="charlson_icd10_quan", assign0=FALSE)
df_cci <- data.frame(person_id = cm$person_id,
                     cci = score(cm, weights="charlson", assign0=FALSE))

# =====================================================================
# 7) 合并所有 -> df  + 派生复合变量
# =====================================================================
df <- df_base |>
  left_join(df_lab,      by="person_id") |>
  left_join(df_comorbid, by="person_id") |>
  left_join(df_drug,     by="person_id") |>
  left_join(df_survey,   by="person_id") |>
  left_join(df_cci,      by="person_id")

df <- df |>
  mutate(
    # 吸烟三分类
    smoke_status = case_when(
      is.na(ever_smoke)            ~ "never",
      ever_smoke & !is.na(curr_smoke) ~ "current",
      ever_smoke                   ~ "former",
      TRUE                         ~ "never"),
    # 饮酒有序等级
    alcohol = factor(case_when(
      alc_raw=="DrinkFrequencyPastYear_Never"          ~ "never",
      alc_raw=="DrinkFrequencyPastYear_MonthlyOrLess"  ~ "monthly_or_less",
      alc_raw=="DrinkFrequencyPastYear_2to4PerMonth"   ~ "2to4_month",
      alc_raw=="DrinkFrequencyPastYear_2to3PerWeek"    ~ "2to3_week",
      alc_raw=="DrinkFrequencyPastYear_4orMorePerWeek" ~ "ge4_week",
      TRUE ~ NA_character_),
      levels=c("never","monthly_or_less","2to4_month","2to3_week","ge4_week"), ordered=TRUE),
    # 复合变量
    ascvd            = pmax(chd, stroke, pvd, na.rm=TRUE),
    antihtn_any      = pmax(acei_arb, bblocker, ccb, diuretic, alpha_blk, other_htn, na.rm=TRUE),
    oral_hypoglyc    = pmax(su, metformin, tzd, glinide, agi, dpp4, sglt2, na.rm=TRUE),
    arm              = factor(arm, levels=c("NON_USER","STATIN"))
  ) |>
  # 简单生理范围清洗(去离群/单位错标)
  mutate(
    tg  = ifelse(tg  > 0 & tg  < 10000, tg,  NA),
    ldl = ifelse(ldl > 0 & ldl < 1000, ldl, NA),
    bmi = ifelse(bmi >= 10 & bmi <= 80, bmi, NA),
    age = ifelse(age >= 18 & age <= 100, age, NA)
  )

cat('\n========== df 维度 ==========\n'); cat(nrow(df),'行 x',ncol(df),'列\n')
cat('两臂:\n'); print(table(df$arm))
cat('\n各连续协变量缺失率:\n')
print(round(colMeans(is.na(df[,c("age","ldl","hdl","tc","tg","bmi","calcium","alt","egfr","cci")])),3))

# =====================================================================
# 8) Table 1 (分臂)
# =====================================================================
if (!requireNamespace("tableone", quietly=TRUE))
  install.packages("tableone", repos="https://cloud.r-project.org")
library(tableone)

cont_vars <- c("age","ldl","hdl","tc","tg","bmi","calcium","alt","egfr","cci")
cat_vars  <- c("sex","race","ethnicity","smoke_status","alcohol",
               "htn","obesity","diabetes","chd","stroke","pvd","afib","kidney","dementia",
               "as_spond","ra","sle","pulmonary","cancer","neuro","ascvd",
               "aspirin","insulin","oral_hypoglyc","su","metformin","tzd","glinide","agi",
               "dpp4","glp1","sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk",
               "other_htn","antihtn_any","antithromb","nsaid")
all_vars <- c(cont_vars, cat_vars)

# 二分类 flag 转 factor 便于 Table1
flag_vars <- intersect(cat_vars, c("htn","obesity","diabetes","chd","stroke","pvd","afib",
  "kidney","dementia","as_spond","ra","sle","pulmonary","cancer","neuro","ascvd",
  "aspirin","insulin","oral_hypoglyc","su","metformin","tzd","glinide","agi","dpp4","glp1",
  "sglt2","acei_arb","bblocker","ccb","diuretic","alpha_blk","other_htn","antihtn_any","antithromb","nsaid"))
df[flag_vars] <- lapply(df[flag_vars], function(x) factor(ifelse(is.na(x),0,x), levels=c(0,1)))

tab1 <- CreateTableOne(vars=all_vars, strata="arm", data=df,
                       factorVars=cat_vars, addOverall=TRUE)
tab1_print <- print(tab1, smd=TRUE, showAllLevels=TRUE,
                    quote=FALSE, noSpaces=TRUE, printToggle=FALSE)
cat('\n========== Table 1 ==========\n'); print(tab1_print)

# 保存
write.csv(tab1_print, "table1_statin_htg.csv")
saveRDS(df, "df_statin_htg_covariates.rds")
cat('\n已保存 table1_statin_htg.csv 与 df_statin_htg_covariates.rds(本地，可 gsutil cp 到 bucket)\n')
