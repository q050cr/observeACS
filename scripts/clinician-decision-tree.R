

# dependencies
library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(glue)
library(stringr)
library(purrr)
library(ggplot2)
library(skimr)
library(tableone)
library(Hmisc)
library(RColorBrewer)

# clinician decision tree (decisions by MM & EG)

# load data ---------------------------------------------------------------
source("./scripts/data_preprocessing_BACC.R")

dat.bacc.oz_clinician <- dat.bacc.oz %>% 
  filter(grace_score< 109)

dat.bacc.oz_clinician <- dat.bacc.oz_clinician %>% 
  mutate(clinician_decision_strict = ifelse( (grace_score < 109 | is.na(grace_score)) & 
                                               (t0_ntbnp_value < 300 | is.na(t0_ntbnp_value)) & 
                                               (t0_crp_value < 10 | is.na(t0_crp_value)) & 
                                               (COPEPTIN_COMBINED < 10 | is.na(COPEPTIN_COMBINED)) &
                                               (t0_ckdepi_value > 30 | is.na(t0_ckdepi_value)) &
                                               (t0_hb_value > 10 | is.na(t0_ckdepi_value)), "pred_surv", "pred_death")) %>% 
  mutate(punkt.bnp = ifelse(t0_ntbnp_value >= 300, 1,0 ),
         punkt.crp = ifelse(t0_crp_value >= 10, 1,0),
         punkt.copeptin = ifelse(COPEPTIN_COMBINED >= 10, 1, 0),
         punkt.gfr = ifelse(t0_ckdepi_value <= 30, 1, 0),
         punkt.hb = ifelse(t0_hb_value <= 10, 1, 0)
  ) %>% 
  rowwise() %>% 
  mutate(sum.punkte = sum(punkt.bnp,punkt.crp,punkt.copeptin, punkt.gfr, punkt.hb , na.rm = TRUE))

dat.bacc.oz_clinician <- dat.bacc.oz_clinician %>% 
  mutate(max1bm.pos = ifelse(sum.punkte<=1, "pred_surv", "pred_death"))

dat.bacc.oz_clinician %>% 
  select(grace_score, t0_ntbnp_value, t0_crp_value,  COPEPTIN_COMBINED, t0_ckdepi_value, t0_hb_value, sum.punkte, max1bm.pos ) -> testdat

#str(dat.bacc.oz_clinician$max1bm.pos)
#sum(dat.bacc.oz_clinician$grace_score > 109)

table(dat.bacc.oz_clinician$clinician_decision_strict, dat.bacc.oz_clinician$o_mortality)
table(dat.bacc.oz_clinician$max1bm.pos, dat.bacc.oz_clinician$o_mortality)
prop.table(table(dat.bacc.oz_clinician$max1bm.pos, dat.bacc.oz_clinician$o_mortality))


dat.bacc.oz_clinician %>% 
  select(rapID, max1bm.pos, t0_ntbnp_value, t0_crp_value, COPEPTIN_COMBINED, t0_ckdepi_value, t0_hb_value) -> test






# sophisticated approach did not work --------------------------------------------------
vars <- c("t0_ntbnp_value", "t0_crp_value", "COPEPTIN_COMBINED", "t0_ckdepi_value", "t0_hb_value")
thresholds <- list(thresh_ntbnp = 300, thresh_crp = 10, thresh_copeptin = 10, thresh_ckdepi = 30, thresh_hb = 10)
comparisons <- list(thresh_ntbnp = `<`, thresh_crp = `<`, thresh_copeptin = `<`, thresh_ckdepi = `>`, thresh_hb = `<`)

dat.bacc.oz_clinician$under_threshold <- apply(dat.bacc.oz_clinician[ , vars], 1, function(row) {
  sum(mapply(function(value, threshold, comparison) {
    comparison(value, threshold) | is.na(value)
  }, as.numeric(row), thresholds, comparisons)) > 1
})

comparison(2, 3)



