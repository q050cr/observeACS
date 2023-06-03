

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

pMiss <- function(x){sum(is.na(x))/length(x)*100}
base::source("./scripts/data_preprocessing1.R")

# update field from data_preprocessing1.R
save.FILE.PREPROCESS <- FALSE

# load datat
# rapid_orig <- clean_names(readxl::read_excel("../data-musti/Rapid_Rule_Out_Komplett.xlsx"))
rapid_orig1 <- clean_names(readxl::read_excel("./data/BACC/BACC_Database 1-2303-nopass.xlsx")) %>% 
  filter(stemi != 1)

# we have to calculate the OZ
rapid_orig <- rapid_orig1 %>% 
  rename(
    t0_hstnt_value = hstnt_0h,
    t1_hstnt_value = hstnt_1h,
    t2_hstnt_value = hstnt_3h) %>% 
  # the dateformat is unambiguous
  mutate(t0_hstnt_effective_date_time = excel_numeric_to_date(as.numeric(as.character(hstnt_0h_dmyhm)), 
                                                              include_time = TRUE, 
                                                              date_system = "modern"),
         t1_hstnt_effective_date_time = excel_numeric_to_date(as.numeric(as.character(hstnt_1h_dmyhm)), 
                                                              include_time = TRUE, 
                                                              date_system = "modern"),
         t2_hstnt_effective_date_time = excel_numeric_to_date(as.numeric(as.character(hstnt_3h_dmyhm)), 
                                                              include_time = TRUE, 
                                                              date_system = "modern")
  )

# clean dat finish ---------------
rapid_orig <- rapid_orig %>% 
  # recode symptombeginn to values from Heidelberg
  mutate(symptombeginn = case_when( symptom_onset_hb == "0" ~ 1,
                                    symptom_onset_hb == "1" ~ 2,
                                    symptom_onset_hb == "2" ~ 3,
                                    symptom_onset_hb == "NA" ~ 0) 
  ) %>% 
  rename(rap_id = external_id, 
         sex_f1_m0 = female,
         t0_krea_value = crea,
         grace_score = grace_db,
         h_hypertonie = hypertension,
         h_diabetes = diabetes, 
         h_cholesterin = hyperlipid, 
         aktiver_raucher = former0_current1_smoker, 
         h_familienana= fam_hx_cad,
         h_infarkt = hx_ami
  )

# subselection 
data01 <- rapid_orig %>% 
  select(rap_id,  age, sex_f1_m0, t0_krea_value, 
         # GRACE
         grace_score, 
         # CV Risk
         h_hypertonie, h_diabetes, h_cholesterin, aktiver_raucher, h_familienana,
         # History
         cad_bypass_pci, ## not distinguished in BACC (h_pci, h_cabg, h_khk)
         h_infarkt, 
         # final ACS diagnoses
         symptombeginn, # not present in BACC: cpu_admission_time, procedere_entlassung,
         # trops
         t0_hstnt_value, t0_hstnt_effective_date_time, 
         t1_hstnt_value, t1_hstnt_effective_date_time,
         t2_hstnt_value, t2_hstnt_effective_date_time
  )

## New Troponin dataframe for classification

c_trops <- paste(
  data01["t0_hstnt_value"][[1]],
  data01["t1_hstnt_value"][[1]],
  data01["t2_hstnt_value"][[1]],
  
  sep = ","
)
c_trops <- str_replace_all(c_trops, pattern=",NA", replacement = "")

time_trops <- paste(
  data01["t0_hstnt_effective_date_time"][[1]],
  data01["t1_hstnt_effective_date_time"][[1]],
  data01["t2_hstnt_effective_date_time"][[1]],
  
  sep = ","
)

# drop patients that have no information on first trop (datetime)
index <- str_starts(time_trops, pattern = "NA", negate = FALSE)
time_trops <- time_trops[!index]
c_trops <- c_trops[!index]
## missing 1st trops 
index1 <- str_starts(c_trops, pattern = "NA", negate = FALSE)
time_trops <- time_trops[!index1]
c_trops <- c_trops[!index1]

time_trops <- str_replace_all(time_trops, pattern=",NA", replacement = "")

data01 <- data01[!index, ]
data01 <- data01[!index1, ]


classification_dat <- data.frame(
  rap_id=data01$rap_id,
  c_trops=c_trops,
  time_trops=time_trops)

classification_dat <- classification_dat %>% 
  separate(c_trops, c("trop_c0", "trop_c1", "trop_c2"), sep =",") %>% 
  separate(time_trops, c("t0_hstnt", "t1_hstnt", "t2_hstnt"), sep = ",") %>% 
  mutate(across(2:4, as.numeric),
         across(5:7, as.POSIXct, tz="CET")) %>% 
  as_tibble()

## Combine with metadata

classification_dat <- classification_dat %>% 
  left_join(data01, by = c("rap_id")) %>% 
  select(-(t0_hstnt_value:t2_hstnt_effective_date_time))

# And we we need to create the necessary cols derived from those above:

classification_dat <- classification_dat %>% 
  mutate(delta_t_min=as.numeric(t1_hstnt-t0_hstnt, units="mins"),
         absolute_delta_c = abs(trop_c1 - trop_c0),
         absolute_delta_c_frac_trop0 = absolute_delta_c/ trop_c0
  ) %>% 
  relocate(symptombeginn, .before = trop_c0) %>% 
  relocate(delta_t_min:absolute_delta_c_frac_trop0, .before = age)

# Program ESC Algorithm ------------

## All Rules Applied

# Code for the different decision algorithms: 0/1h, 0/2h, 0/3h algorithm. The preferred algorithm is the 0/1h algorithm. Also important is the variable encoding for `symptombeginn`: $1 = 0-3 h$, $2 = 3-6 h$, $3 \ge 6 h$, $0 = unbekannt$. 


# SELECT ESC Algorithm and prefer 0/1h 
applied_rule <- c()
applied_rule1 <- c()

for(i in 1:nrow(classification_dat)){ 
  # seq through every observation, highest priority conditions at beginning combined with `next`statement
  
  ### ESC 0h
  if (classification_dat[[i, "symptombeginn"]] >=2 & 
      classification_dat[[i, "trop_c0"]] < 5) {
    applied_rule[i] <- "Rule-out (ESC 0h)"
    applied_rule1[i] <- "Rule-out"
    next
  }
  # if (classification_dat[[i, "symptombeginn"]] ==3 & 
  #     classification_dat[[i, "trop_c0"]] <= 14) {
  #   applied_rule[i] <- "Rule-out (ESC 0h)"
  #   applied_rule1[i] <- "Rule-out"
  #   next
  # }
  if (classification_dat[[i, "trop_c0"]] >= 52) {
    applied_rule[i] <- "Rule-in (ESC 0h)"
    applied_rule1[i] <- "Rule-in"
    next
  }
  
  # no second troponin (e.g. due to discharge against medical advice)
  if ( is.na(classification_dat[[i, "delta_t_min"]]) | 
       is.na(classification_dat[[i, "trop_c1"]])                       ### NEEDS TO BE ADDED in BACC COHORT (e.g. id BACC_1525 no second Trop but Time)
  ) {
    applied_rule[i] <- "No second Trop Measurement available"
    applied_rule1[i] <- "No second Trop Measurement available"
    next
  }
  
  # insufficient time
  if (classification_dat[[i, "delta_t_min"]] <= 30) {
    applied_rule[i] <- "Insufficient time interval between samples"
    applied_rule1[i] <- "Insufficient time interval between samples"
    next
  }
  
  ### ESC 0/1h
  if (classification_dat[[i, "delta_t_min"]] > 30 & classification_dat[[i, "delta_t_min"]] <= 90 &
      classification_dat[[i, "trop_c0"]] < 12 & classification_dat[[i, "absolute_delta_c"]] < 3) {
    applied_rule[i] <- "Rule-out (ESC 0/1h)"
    applied_rule1[i] <- "Rule-out"
    next
  }
  if (classification_dat[[i, "delta_t_min"]] > 30 & classification_dat[[i, "delta_t_min"]] <= 90 &
      classification_dat[[i, "absolute_delta_c"]] >= 5) {
    applied_rule[i] <- "Rule-in (ESC 0/1h)"
    applied_rule1[i] <- "Rule-in"
    next
  }
  # observe
  if ( (classification_dat[[i, "delta_t_min"]] > 30 & classification_dat[[i, "delta_t_min"]] <= 90) &
       ( classification_dat[[i, "trop_c0"]] >= 12 | 
         classification_dat[[i, "absolute_delta_c"]] >= 3 ) &
       classification_dat[[i, "absolute_delta_c"]] < 5 ) {
    applied_rule[i] <- "Observation zone (ESC 0/1h)"
    applied_rule1[i] <- "Observation zone"
    next
  }
  
  ### ESC 0/2h
  if (classification_dat[[i, "delta_t_min"]] > 90 & classification_dat[[i, "delta_t_min"]] <= 150 &
      classification_dat[[i, "trop_c0"]] < 14 & classification_dat[[i, "trop_c1"]] < 14 & classification_dat[[i, "absolute_delta_c"]] < 4) {
    applied_rule[i] <- "Rule-out (ESC 0/2h, Reichlin 2015)"
    applied_rule1[i] <- "Rule-out"
    next
  }
  if (classification_dat[[i, "delta_t_min"]] > 90 & classification_dat[[i, "delta_t_min"]] <= 150 &
      classification_dat[[i, "trop_c0"]] >= 52 & classification_dat[[i, "trop_c1"]] >= 52 ) {
    applied_rule[i] <- "Rule-in (ESC 0/2h, Reichlin 2015)"
    applied_rule1[i] <- "Rule-in"
    next
  } 
  if (classification_dat[[i, "delta_t_min"]] > 90 & classification_dat[[i, "delta_t_min"]] <= 150 &
      classification_dat[[i, "absolute_delta_c"]] >= 10) {
    applied_rule[i] <- "Rule-in (ESC 0/2h, Reichlin 2015)"
    applied_rule1[i] <- "Rule-in"
    next
  }   
  # observe
  if ( (classification_dat[[i, "delta_t_min"]] > 90 & classification_dat[[i, "delta_t_min"]] <= 150) &
       (classification_dat[[i, "trop_c0"]] >= 14 | classification_dat[[i, "trop_c1"]] >= 14 | classification_dat[[i, "absolute_delta_c"]] >= 4) &
       (classification_dat[[i, "trop_c0"]] < 52 | classification_dat[[i, "trop_c1"]] < 52 ) &
       classification_dat[[i, "absolute_delta_c"]] < 10 ) {
    applied_rule[i] <- "Observation zone (ESC 0/2h, Reichlin 2015)"
    applied_rule1[i] <- "Observation zone"
    next
  }
  
  ### ESC 0/3h
  if (classification_dat[[i, "delta_t_min"]] > 150 &
      classification_dat[[i, "trop_c0"]] <= 14 & 
      (classification_dat[[i, "trop_c1"]] <= 14 | 
       classification_dat[[i, "absolute_delta_c"]] <= 7) ) {
    applied_rule[i] <- "Rule-out (ESC 0/3h)"
    applied_rule1[i] <- "Rule-out"
    next
  }
  if (classification_dat[[i, "delta_t_min"]] > 150 &
      classification_dat[[i, "trop_c0"]] > 14 & 
      (classification_dat[[i, "trop_c1"]] <= 14 | 
       classification_dat[[i, "absolute_delta_c_frac_trop0"]] <= 0.2) ) {
    applied_rule[i] <- "Rule-out (ESC 0/3h)"
    applied_rule1[i] <- "Rule-out"
    next
  }
  if (classification_dat[[i, "delta_t_min"]] > 150 &
      (classification_dat[[i, "trop_c0"]] <= 14 & classification_dat[[i, "trop_c1"]] > 14) & 
      classification_dat[[i, "absolute_delta_c"]] >  7) {
    applied_rule[i] <- "Rule-in (ESC 0/3h)"
    applied_rule1[i] <- "Rule-in"
    next
  }
  if (classification_dat[[i, "delta_t_min"]] > 150 &
      (classification_dat[[i, "trop_c0"]] > 14 & classification_dat[[i, "trop_c1"]] > 14) & 
      classification_dat[[i, "absolute_delta_c_frac_trop0"]] >  0.2) {
    applied_rule[i] <- "Rule-in (ESC 0/3h)"
    applied_rule1[i] <- "Rule-in"
    next
  }
}

classified_df <- as_tibble(cbind(classification_dat, applied_rule, applied_rule1)) %>% 
  mutate_if(is.numeric, round, digits=2) 

#colnames(classified_df)
#colnames(rapid_orig)

dat.bacc.oz <- rapid_orig %>% 
  left_join(classified_df %>% 
              select(rap_id, applied_rule, applied_rule1, absolute_delta_c),
            by=c("rap_id"="rap_id")) %>% 
  filter(applied_rule1=="Observation zone")

#colnames(dat.bacc.oz)
#colnames(dat_train)


# harmonize datasets  ------------

# rename variables
dat.bacc.oz <- dat.bacc.oz %>% 
  mutate(ekg_sinus_normal = ifelse(st_depression_0h_n1==0 & ecg_schenkelblock == 0 &
                                     ecg_vensms == 0 & t_wave_inversion_0h_n1 == 0 &
                                     ecg_at == 0, 1, 0)
         ) %>% 
  rename(
    rapID = rap_id,
    o_mortality = qc_fu_death,  
    grace_score = grace_score, 
    h_diabetes = h_diabetes,
    h_hypertonie = h_hypertonie, 
    h_cholesterin = h_cholesterin, 
    aktiver_raucher = aktiver_raucher,
    h_familienana = h_familienana,  
    h_khk = cad_bypass_pci,  # assuming 'cad_bypass_pci' corresponds to 'h_khk' in df1
    h_lvdys_grad_BINARY = hf,  # here: "LV-Dysfunktion vor Aufnahme"
    vit_rr_syst = sys_bp_1,  # assuming 'sys_bp_1' corresponds to 'vit_rr_syst' in df1
    vit_herzfrequenz = heart_rate_1, 
    # 
    ekg_st_senkung = st_depression_0h_n1,  # assuming 'st_depression_0h_n1' corresponds to 'ekg_st_senkung' in df1
    ekg_schrittmacher = ecg_vensms,  # assuming 'ecg_at' corresponds to 'ekg_schrittmacher' in df1
    # labs
    t0_hstnt_value = t0_hstnt_value, 
    t0_ntbnp_value = nt_probnp_bacc,  # assuming 'nt_probnp_bacc' corresponds to 't0_ntbnp_value' in df1
    COPEPTIN_COMBINED = copeptin,
    t0_krea_value = t0_krea_value, 
    t0_ckdepi_value = gfr,
    t0_crp_value = crp, 
    t0_leuko_value = leukocytes, 
    #t0_hst_value  not included, 
    t0_na_value = natrium, 
    t0_gluc_value = glucose, 
    t0_ck_value = ck_0h, 
    t0_hb_value = hemoglobine, 
    t0_thrombo_value = thrombocytes, 
    t0_inr_value = inr
  )# %>% 
  # # missing values
  # mutate(
  #   vit_temperatur = NA,  # no apparent corresponding column in df2
  #   vit_saettigung = NA,  # no apparent corresponding column in df2
  #   vit_atemfrequenz = NA,  # no apparent corresponding column in df2
  #   delta_t_first = NA,  # no apparent corresponding column in df2
  #   t0_k_value = NA,  # no apparent corresponding column in df2
  #   t0_ldh_value = NA,  # no apparent corresponding column in df2
  #   t0_got_value = NA,  # no apparent corresponding column in df2
  #   t0_gpt_value = NA,  # no apparent corresponding column in df2
  #   t0_hkt_value = NA,  # no apparent corresponding column in df2
  #   t0_quick_value = NA,  # no apparent corresponding column in df2
  # )

# further cleaning -------------------


# convert to factor
catVars <- c(
  # demographics & history category
  "sex_f1_m0", 
  "symptombeginn","symptomtyp",
  "crowding_low0_norm1_high2","h_diabetes","h_hypertonie","h_cholesterin",
  "aktiver_raucher","h_familienana","h_infarkt","h_cabg","h_khk",
  "h_vessel_disease","h_pci","h_lvdysfunktion","h_lvdys_grad",
  
  # EKG category
  "ekg_t_negativierung","ekg_schrittmacher","ekg_atriale_tachy","ekg_block_ohne_sm",
  "ekg_st_senkung","ekg_sinusrhythmus","ekg_sinus_normal",
  
  # hospitalisation category
  "aufnahme_krankenhaus",
  # outcome category
  "mortality_30d","mortality_90d","o_mortality","o_stroke","o_cabg",
  "o_reinfarction","o_recoro","o_re_ptca","hk_within_30_days",
  # trop category 
  "t0_hstnt_categ_num",
  # scores
  "grace_low0int1high2","KHK__Killip_Class"
)

# dat.observe[catVars] <- lapply(dat.observe[catVars], factor) # not applied, otherwise problems later with replace_na()

# create some more vectors to investigate separately
demographics_cat <- c("sex_f1_m0", "symptombeginn","symptomtyp",
                      "crowding_low0_norm1_high2","h_diabetes","h_hypertonie","h_cholesterin",
                      "aktiver_raucher","h_familienana","h_infarkt","h_cabg","h_khk",
                      "h_vessel_disease","h_pci","h_lvdysfunktion","h_lvdys_grad")
ecg_cat <- c("ekg_t_negativierung","ekg_schrittmacher","ekg_atriale_tachy","ekg_block_ohne_sm",
             "ekg_st_senkung","ekg_sinusrhythmus","ekg_sinus_normal")
outcomes_cat <- c("aufnahme_krankenhaus","mortality_30d","mortality_90d","o_mortality","o_stroke","o_cabg",
                  "o_reinfarction","o_recoro","o_re_ptca", "hk_within_30_days")
scores_cat <- c("grace_low0int1high2","KHK__Killip_Class")

change2numericBACC <- c("t0_hstnt_value", "t0_ntbnp_value", "t0_tsh_value", "t0_dd_value", "delta_c_first_tnt_relativ", "o_time_re_ptca",
                        "COPEPTIN_COMBINED",
                    "o_time_recoro", "o_time_reinfarction", "o_time_cabg", "o_time_stroke", "o_time_mortality",
                    "t0_thrombo_value", "t0_hb_value", "t0_ck_value", "t0_gluc_value", "t0_ckdepi_value",
                    "t0_na_value", "t0_hst_value", "t0_leuko_value","t0_crp_value", "t0_krea_value",
                    "vit_rr_syst", "t0_inr_value", "vit_herzfrequenz")

dat.bacc.oz <- dat.bacc.oz %>% 
  mutate(across(any_of(change2numericBACC), as.numeric)) %>% 
  # first convert to numeric to prevent that "NA" gets a third factor
  mutate(across(any_of(c(catVars, outcomes_cat, demographics_cat, ecg_cat)), as.numeric)) %>% 
  mutate(across(any_of(c(catVars, outcomes_cat, demographics_cat, ecg_cat)), as.factor))

fctrs <- c("o_mortality", "sex_f1_m0", "symptombeginn", 
           "h_diabetes", "h_hypertonie", "h_cholesterin", "aktiver_raucher", 
           "h_familienana", "h_khk", "h_lvdys_grad",
           "ekg_sinus_normal", "ekg_st_senkung", "ekg_schrittmacher")

dat.bacc.oz <- dat.bacc.oz %>% 
  mutate(across(any_of(fctrs), as.numeric)) %>% 
  mutate(across(any_of(fctrs), as.factor)) %>% 
  mutate(o_mortality = factor(o_mortality, labels = c("survived", "died")))

# factors must be harmonized, different levels (1,2 vs 0,1)
#str(model_data1)
dat.bacc.oz$sex_f1_m0 <- factor(dat.bacc.oz$sex_f1_m0, labels = c(0,1))
dat.bacc.oz$symptombeginn <- factor(dat.bacc.oz$symptombeginn, labels = c(0,1,2,3))
dat.bacc.oz$h_diabetes <- factor(dat.bacc.oz$h_diabetes, labels = c(0,1))
dat.bacc.oz$h_hypertonie <- factor(dat.bacc.oz$h_hypertonie, labels = c(0,1))
dat.bacc.oz$h_cholesterin <- factor(dat.bacc.oz$h_cholesterin, labels = c(0,1))
dat.bacc.oz$aktiver_raucher <- factor(dat.bacc.oz$aktiver_raucher, labels = c(0,1))
dat.bacc.oz$h_familienana <- factor(dat.bacc.oz$h_familienana, labels = c(0,1))
dat.bacc.oz$h_khk <- factor(dat.bacc.oz$h_khk, labels = c(0,1))

dat.bacc.oz$ekg_st_senkung <- factor(dat.bacc.oz$ekg_st_senkung, labels = c(0,1))
dat.bacc.oz$ekg_schrittmacher <- factor(dat.bacc.oz$ekg_schrittmacher, labels = c(0,1))
dat.bacc.oz$ekg_sinus_normal <- factor(dat.bacc.oz$ekg_sinus_normal, labels = c(0,1))

#model_data1$symptombeginn
# NAs? -----------------------------------------------------------------
# show vals that have > 40% missing data

pMiss_bacc.oz <- apply(dat.bacc.oz,2,pMiss)

# SAVE BACC
if (save.FILE.PREPROCESS == TRUE) {
  saveRDS(object = dat.bacc.oz, file = glue("./output/Rdata/cleaned-dat/{Sys.Date()}-dat.bacc.oz-cleaned.rds"))
  #openxlsx::write.xlsx(x = dat.bacc.oz, file = glue("./output/Rdata/cleaned-dat/{Sys.Date()}-dat.bacc.oz-cleaned.xlsx"))
}



# COMBINE BACC and UKHD DATA ------------------------------------------
# Find common column names
# Add identifier column to each dataset

model_data1_ukhd <- model_data1
model_data1_ukhd$source <- "UKHD"
model_data1_ukhd$rapID <- paste0("UKHD_", model_data1_ukhd$rapID)
dat.bacc.oz$source <- "BACC"

# no lv-dysfunction in BACC, but reduced LVEF binary
dat.bacc.oz$h_lvdys_grad_BINARY <- factor(dat.bacc.oz$h_lvdys_grad_BINARY, labels=c("normal", "reduced"))
model_data1_ukhd$h_lvdys_grad_BINARY <- factor(ifelse(model_data1_ukhd$h_lvdys_grad != "good", 1, 0), labels=c("normal", "reduced"))
common_columns <- intersect(colnames(model_data1_ukhd), colnames(dat.bacc.oz))

# Perform row binding with common columns only
combined_dataset <- rbind(model_data1_ukhd[, common_columns], dat.bacc.oz[, common_columns])


if (save.FILE.PREPROCESS == TRUE) {
  saveRDS(object = combined_dataset, file = glue("./output/Rdata/cleaned-dat/{Sys.Date()}-combined_dataset_ukhd_bacc.rds"))
  #openxlsx::write.xlsx(x = combined_dataset, file = glue("./output/Rdata/cleaned-dat/{Sys.Date()}-combined_dataset_ukhd_bacc.xlsx"))
}

# sanity trop check
sanity_check <- dat.bacc.oz %>% 
  left_join(rapid_orig1 %>% select(external_id, hstnt_0h), by = c("rapID"="external_id"))

