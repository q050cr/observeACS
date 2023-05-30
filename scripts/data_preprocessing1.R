
# data description --------------------------------------------------------
'
#
# 1) dat.sav :        original data imported from SPSS data file (needed convertion from medcalc)
#
# 2) dat.observe :    original data filtered patients in "observe zone" for ACS
#
# 3) dat :            data cleaned and edited obvious NAs, dropped cols with >30% missing data  <--- USE THIS in LAST ANALYSIS for Modelling
#
# 4) dat.outcome :    "dat" split in df only containing outcome (y) variables
#
# 5) dat.independent: "dat" split in df only containing independent (X) variables
#
#
'

# std libs
library(haven)
library(dplyr)
library(lubridate)

# define fns
## fn to calculate perc missing
pMiss <- function(x){sum(is.na(x))/length(x)*100}

# import data -------------------------------------------------------------
#dat.sav <- haven::read_sav("./data/Observe_ML_Brahms.sav")

# NEW 2023-05: mastertabelle musti (different naming)
dat.sav <- readxl::read_excel(path = "./data/202305master/230308_Observe_Zone_Mastertabelle.xlsx")

# rename some cols with dplyr -----------------------------------------------
dat.sav <- dat.sav %>% 
  dplyr::rename(
    symptombeginn="symptombeginn__1___0_3h__2___3_6h__3____6h__0___unbekannt_nicht_einzugrenzen_",
    symptomtyp = "symptomtyp_1_Chest_Pain__2_Dyspnoe__3_Atyp",
    h_lvdys_grad = "h_lvdysfunktion_quant_1_leichtgradig__45_54%___2_mittelgradig__30_44%___3_hochgradig___30%_"
  ) %>% 
  # use lubridate to clean datetime
  dplyr::mutate(CPU_AdmissionDate=dmy(CPU_AdmissionDate))

# filter patients -----------------------------------------------------------
dat.observe <- dat.sav %>%
  dplyr::filter(Observe==1)

# select cols oi ------------------------------------------------------------
dat.observe %>% 
  dplyr::select(
    # demographics
    rapID, age, sex_f1_m0, 
    CPU_AdmissionDate:CPU_AdmissionTime,symptombeginn:symptomtyp,
    crowding_low0_norm1_high2,
    # history
    h_diabetes:h_lvdys_grad,
    # vital parameter
    vit_rr_diast:vit_atemfrequenz,
    # ecg
    ekg_t_negativierung:ekg_sinus_normal,
    ## ICD_Entladung ..0, AVNRT..0, Dialyse ..0, Heart_transplant..0
    # hospitalisation
    bed_time_h, aufnahme_krankenhaus,
    # outcomes
    mortality_30d:Bleed_BARC_type,
    # troponin and cinetics
    t0_hstnt_categ_num, t0_hstnt_t_diff,  # WHAT IS THIS?
    t0_hstnt_value, t1_hstnt_value, t2_hstnt_value, t3_hstnt_value,
    delta_t_first, delta_c_first_tnt_relativ, delta_t_first,
    # lab values .. c0_tnt == t0_hstnt_value - already added that
    t0_ckdepi_value:t0_inr_value, 
    COPEPTIN_COMBINED, t0_ntbnp_value, t0_dd_value, t0_tsh_value,
    # was hk performed (only for descriptive stats, excluded for modeling)
    hk_within_30_days,  # but they could have also died due to a "complication", maybe also include in a model ;)
    CPU_AdmissionDate,  # CPU_AdmissionTime only NA
    # GRACE score
    grace_score, grace_low0int1high2,  
    KHK__Killip_Class  # this Killip Class comes later in the df and has fewer NA
  ) -> dat.observe


# variable encoding -------------------------------------------------------

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



# data manipulation -------------------------------------------------------

# manipulate patient with HF=734/min
index <- which(dat.observe["vit_herzfrequenz"]>300)  # one patient with index 318 in dat.observe
dat.observe[index,"vit_herzfrequenz"] <- 134

# convert factor vars
# dat.observe[catVars] <- lapply(dat.observe[catVars], factor) # not applied, otherwise problems later with replace_na()


change2numeric <- c("t0_ntbnp_value", "t0_tsh_value", "t0_dd_value", "delta_c_first_tnt_relativ", "o_time_re_ptca",
                    "o_time_recoro", "o_time_reinfarction", "o_time_cabg", "o_time_stroke", "o_time_mortality")

dat.observe <- dat.observe %>% 
  mutate(across(all_of(change2numeric), as.numeric)) %>% 
  mutate(across(all_of(catVars), as.character)) %>% 
  mutate(across(all_of(outcomes_cat), as.factor)) %>% 
  mutate(across(all_of(demographics_cat), as.factor)) %>% 
  mutate(across(all_of(ecg_cat), as.factor)) -> test

# numeric cols
numcols <- unlist(lapply(dat.observe, is.numeric))
dat.observe_numeric <- dat.observe[ , numcols]


# NA edit -----------------------------------------------------------------

# show vals that have > 40% missing data
pMiss_observe <- apply(dat.observe,2,pMiss)
index <- pMiss_observe > 30

dat.observe %>%
  # replace_na works with "specified values" only
  tidyr::replace_na(
    list(
      # ekg_st_hebung = 0,
      # ekg_schenkelblock = 0,
      o_stroke = "0",
      o_cabg="0",
      o_reinfarction="0",
      o_recoro="0",
      o_re_ptca ="0",
      # # we also have to edit the outcome data and assume that patients where NA is observed in mortality -> no mortality
      mortality_30d = "0",
      mortality_90d="0",
      o_mortality="0",
      aufnahme_krankenhaus="0"
    )
  ) %>% 
  # to replace with values from column we use mutate and ifelse
  # here we can add time until death as those patients definitely were not
  # observed longer and would not be counted as stroke
  mutate(
    o_time_stroke = ifelse(is.na(o_time_stroke),o_time_mortality, o_time_stroke),
    o_time_cabg = ifelse(is.na(o_time_cabg),o_time_mortality, o_time_cabg),
    o_time_reinfarction = ifelse(is.na(o_time_reinfarction),o_time_mortality, o_time_reinfarction),
    o_time_recoro = ifelse(is.na(o_time_recoro),o_time_mortality, o_time_recoro),
    o_time_re_ptca = ifelse(is.na(o_time_re_ptca),o_time_mortality, o_time_re_ptca)
  ) -> dat

pMiss <- function(x){sum(is.na(x))/length(x)*100}
# show vals that have > 30% missing data
pMiss_observe <- apply(dat,2,pMiss)
index <- pMiss_observe > 30

# drop-cols ----------------------------------------------------------------
dat %>%
  mutate(h_lvdys_grad = case_when(
    h_lvdys_grad == 0 ~ "good",
    h_lvdys_grad == 0.5 ~ "slightly reduced",
    h_lvdys_grad == 1 ~ "slightly reduced",
    h_lvdys_grad == 1.5 ~ "moderately reduced",
    h_lvdys_grad == 2 ~ "moderately reduced",
    h_lvdys_grad == 2.5 ~ "severely reduced",
    h_lvdys_grad == 3 ~ "severely reduced"
  )) %>% 
  mutate(h_lvdys_grad = factor(h_lvdys_grad, levels=c("good", "slightly reduced", "moderately reduced", "severely reduced"))) %>% 
  # select cols that have not too many NAs (> 50% of full data) --> copeptin e.g. dropped
  select_if(function(col) sum(!is.na(col))/length(col)*100 > 50) -> dat


# SAVE

saveRDS(object = dat, file = paste0("./output/Rdata/cleaned-dat/", Sys.Date(), "-dat_cleaned_observe2023.rds"))


rm(list=c("dat.sav", "index", "numcols", "pMiss_observe"))


