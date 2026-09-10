-- ============================================================================
-- supp_tg_stageA_concepts.sql
-- Stage A of the TG-change supplementary analysis: resolve the serum
-- triglyceride (TG) concept set from four fixed LOINC seeds and print a
-- review table (concept_id, concept_name, ... , per-concept record counts).
--
-- Replacement tokens (Workbench / BigQuery):
--   {cdr}  -> the AoU CDR dataset, e.g. `<PROJECT>.AllOfUs.<CDR_VERSION>`
--
-- ONLY LOINC seed CODES are hard-coded here (task-mandated):
--   2571-8, 3043-7, 12951-0, 14927-8
-- No OMOP concept_id is hard-coded.
--
-- Run this first; REVIEW the output, then copy the approved concept_ids to
-- codelists/tg_concepts_confirmed.csv before proceeding to Stage C.
-- ============================================================================

WITH
-- (1) the four LOINC seed rows themselves (their own standard-concept ids)
seed AS (
  SELECT c.concept_id, c.concept_name, c.vocabulary_id, c.concept_code,
         c.standard_concept, c.domain_id
  FROM `{cdr}.concept` c
  WHERE c.vocabulary_id = 'LOINC'
    AND c.concept_code IN ('2571-8', '3043-7', '12951-0', '14927-8')
),
-- (2) concepts mapped FROM each seed (concept_relationship)
mapped AS (
  SELECT DISTINCT
         cr.concept_id_2 AS concept_id,
         c.concept_name, c.vocabulary_id, c.concept_code,
         c.standard_concept, c.domain_id
  FROM `{cdr}.concept` seed
  JOIN `{cdr}.concept_relationship` cr
       ON cr.concept_id_1 = seed.concept_id
      AND cr.relationship_id IN ('Maps to', 'Maps to value')
  JOIN `{cdr}.concept` c ON c.concept_id = cr.concept_id_2
  WHERE seed.vocabulary_id = 'LOINC'
    AND seed.concept_code IN ('2571-8', '3043-7', '12951-0', '14927-8')
),
-- (3) full descendant subtree of each seed (concept_ancestor)
descendants AS (
  SELECT DISTINCT
         ca.descendant_concept_id AS concept_id,
         c.concept_name, c.vocabulary_id, c.concept_code,
         c.standard_concept, c.domain_id
  FROM `{cdr}.concept` seed
  JOIN `{cdr}.concept_ancestor` ca ON ca.ancestor_concept_id = seed.concept_id
  JOIN `{cdr}.concept` c ON c.concept_id = ca.descendant_concept_id
  WHERE seed.vocabulary_id = 'LOINC'
    AND seed.concept_code IN ('2571-8', '3043-7', '12951-0', '14927-8')
),
candidate AS (
  SELECT * FROM seed
  UNION DISTINCT
  SELECT * FROM mapped
  UNION DISTINCT
  SELECT * FROM descendants
)
SELECT c.concept_id, c.concept_name, c.vocabulary_id, c.concept_code,
       c.standard_concept, c.domain_id,
       COUNT(DISTINCT m.person_id) AS n_persons,
       COUNT(*)                    AS n_records
FROM candidate c
LEFT JOIN `{cdr}.measurement` m ON m.measurement_concept_id = c.concept_id
WHERE IFNULL(c.standard_concept, 'S') = 'S'
  AND IFNULL(c.domain_id, 'Measurement') = 'Measurement'
  AND c.concept_id IS NOT NULL
GROUP BY 1, 2, 3, 4, 5, 6
ORDER BY n_records DESC;

-- ============================================================================
-- Review target: a compact table (concept_id, concept_name, code, n_records).
-- Approve ONLY rows that are serum-triglyceride measurements; then save them
-- to codelists/tg_concepts_confirmed.csv and set CONFIRM_TG_CONCEPTS=TRUE in
-- supp_tg_change.R before re-running.
-- ============================================================================