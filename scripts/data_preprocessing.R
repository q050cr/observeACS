
# data description --------------------------------------------------------
'
#
# 1) dat.sav :        original data imported from SPSS data file (needed convertion from medcalc)
#
# 2) dat.observe :    original data filtered patients in "observe zone" for ACS
#
# 3) dat :            data cleaned and edited obvious NAs, dropped cols with >30% missing data
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
dat.sav <- haven::read_sav("./data/Observe_ML_Brahms.sav")

#2023-05: mastertabelle musti
# dat.sav <- readxl::read_excel(path = "./data/202305master/230308_Observe_Zone_Mastertabelle.xlsx")

# rename some cols with dplyr -----------------------------------------------
dat.sav <- dat.sav %>% 
  dplyr::rename(
    V_symptombeginn="V_symptombeginn__1___0_3h__2___3_6h__3____6h__0___unbekannt_nic",
    V_symptomtyp = "V_symptomtyp_1_Chest_Pain__2_Dyspnoe__3_Atyp",
    V_h_lvdys_grad = "V@AL"
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
    V_rapID, V_age, V_sex_f1_m0, 
    V_admission_month:V_admission_hour,V_symptombeginn:V_symptomtyp,
    V_crowding_low0_norm1_high2,
    # history
    V_h_diabetes:V_h_lvdysfunkt_ab_mittelgradig,
    # vital parameter
    V_vit_rr_diast:V_vit_atemfrequenz,
    # ecg
    V_ekg_t_negativierung:V_ekg_schenkelblock,
    ## ICD_Entladung ..0, AVNRT..0, Dialyse ..0, Heart_transplant..0
    # hospitalisation
    V_bed_time_h, V_aufnahme_krankenhaus,
    # outcomes
    V_mortality_30d:V_o_time_re_ptca,
    # troponin and cinetics
    V_t0_hstnt_categ_num, V_t0_hstnt_t_diff,  # WHAT IS THIS?
    V_t0_hstnt_value, V_t1_hstnt_value, V_t2_hstnt_value, V_t3_hstnt_value,
    V_delta_t_first, V_delta_c_first_tnt_relativ, V_delta_t_first,
    # lab values .. V_c0_tnt == V_t0_hstnt_value - already added that
    V_t0_ckdepi_value:V_t0_inr_value, 
    COPEPTIN_COMBINED, V_t0_ntbnp_value, V_t0_dd_value, V_t0_tsh_value,
    # was hk performed (only for descriptive stats, excluded for modeling)
    V_hk_within_30_days,  # but they could have also died due to a "complication", maybe also include in a model ;)
    CPU_AdmissionDate,  # CPU_AdmissionTime only NA
    # TIMI scores
    KHK__TIMINSTEMI, KHK__TIMISTEMI,
    # GRACE score
    V_grace_score, V_grace_low0int1high2, Killip_Class, 
    KHK__Killip_Class  # this Killip Class comes later in the df and has fewer NA
    ) -> dat.observe


# variable encoding -------------------------------------------------------

# convert to factor
catVars <- c(
  # demographics & history category
  "V_sex_f1_m0", "V_admission_month", "V_admission_week", "V_admission_weekday", 
  "V_admission_hour","V_symptombeginn","V_symptom_thoraxschmerz","V_symptomtyp",
  "V_crowding_low0_norm1_high2","V_h_diabetes","V_h_hypertonie","V_h_cholesterin",
  "V_aktiver_raucher","V_h_familienana","V_h_infarkt","V_h_cabg","V_h_khk",
  "V_h_vessel_disease","V_h_pci","V_h_lvdysfunktion","V_h_lvdys_grad",
  "V_h_lvdysfunkt_ab_mittelgradig",
  # EKG category
  "V_ekg_t_negativierung","V_ekg_schrittmacher","V_ekg_atriale_tachy","V_ekg_block_ohne_sm",
  "V_ekg_st_senkung","V_ekg_sinusrhythmus","V_ekg_sinus_normal","V_ekg_st_hebung",
  "V_ekg_schenkelblock",
  # hospitalisation category
  "V_aufnahme_krankenhaus",
  # outcome category
  "V_mortality_30d","V_mortality_90d","V_o_mortality","V_o_stroke","V_o_cabg",
  "V_o_reinfarction","V_o_recoro","V_o_re_ptca","V_hk_within_30_days",
  # trop category 
  "V_t0_hstnt_categ_num",
  # scores
  "KHK__TIMINSTEMI","KHK__TIMISTEMI","V_grace_low0int1high2","Killip_Class","KHK__Killip_Class"
  )

# dat.observe[catVars] <- lapply(dat.observe[catVars], factor) # not applied, otherwise problems later with replace_na()

# create some more vectors to investigate separately
demographics_cat <- c("V_sex_f1_m0", "V_admission_month", "V_admission_week", "V_admission_weekday", 
                      "V_admission_hour","V_symptombeginn","V_symptom_thoraxschmerz","V_symptomtyp",
                      "V_crowding_low0_norm1_high2","V_h_diabetes","V_h_hypertonie","V_h_cholesterin",
                      "V_aktiver_raucher","V_h_familienana","V_h_infarkt","V_h_cabg","V_h_khk",
                      "V_h_vessel_disease","V_h_pci","V_h_lvdysfunktion","V_h_lvdys_grad",
                      "V_h_lvdysfunkt_ab_mittelgradig")
ecg_cat <- c("V_ekg_t_negativierung","V_ekg_schrittmacher","V_ekg_atriale_tachy","V_ekg_block_ohne_sm",
               "V_ekg_st_senkung","V_ekg_sinusrhythmus","V_ekg_sinus_normal","V_ekg_st_hebung",
               "V_ekg_schenkelblock")
outcomes_cat <- c("V_aufnahme_krankenhaus","V_mortality_30d","V_mortality_90d","V_o_mortality","V_o_stroke","V_o_cabg",
                  "V_o_reinfarction","V_o_recoro","V_o_re_ptca", "V_hk_within_30_days")
scores_cat <- c("KHK__TIMINSTEMI","KHK__TIMISTEMI","V_grace_low0int1high2","Killip_Class","KHK__Killip_Class")



# data manipulation -------------------------------------------------------

# manipulate patient with HF=734/min
index <- which(dat.observe["V_vit_herzfrequenz"]>300)  # one patient with index 318 in dat.observe
dat.observe[index,"V_vit_herzfrequenz"] <- 134

# convert factor vars
# dat.observe[catVars] <- lapply(dat.observe[catVars], factor) # not applied, otherwise problems later with replace_na()

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
     # V_ekg_st_hebung = 0,
     # V_ekg_schenkelblock = 0,
      V_o_stroke = 0,
      V_o_cabg=0,
      V_o_reinfarction=0,
      V_o_recoro=0,
      V_o_re_ptca =0,
      # # we also have to edit the outcome data and assume that patients where NA is observed in mortality -> no mortality
      V_mortality_30d = 0,
      V_mortality_90d=0,
      V_o_mortality=0,
      V_aufnahme_krankenhaus=0
    )
  ) %>% 
  # to replace with values from column we use mutate and ifelse
  # here we can add time until death as those patients definitely were not
  # observed longer and would not be counted as stroke
  mutate(
    V_o_time_stroke = ifelse(is.na(V_o_time_stroke),V_o_time_mortality, V_o_time_stroke),
    V_o_time_cabg = ifelse(is.na(V_o_time_cabg),V_o_time_mortality, V_o_time_cabg),
    V_o_time_reinfarction = ifelse(is.na(V_o_time_reinfarction),V_o_time_mortality, V_o_time_reinfarction),
    V_o_time_recoro = ifelse(is.na(V_o_time_recoro),V_o_time_mortality, V_o_time_recoro),
    V_o_time_re_ptca = ifelse(is.na(V_o_time_re_ptca),V_o_time_mortality, V_o_time_re_ptca)
  ) -> dat

pMiss <- function(x){sum(is.na(x))/length(x)*100}
# show vals that have > 30% missing data
pMiss_observe <- apply(dat,2,pMiss)
index <- pMiss_observe > 30

# drop-cols ----------------------------------------------------------------
dat %>%
  # select cols that have not too many NAs (> 70% of full data) --> copeptin e.g. dropped
  select_if(function(col) sum(!is.na(col))/length(col)*100 > 70) -> dat



# split outcome data (y) and independent data (X) -------------------------

dat.outcome <- dat %>% 
  select(V_bed_time_h:V_o_time_re_ptca, V_hk_within_30_days)

dat.independent <- dat %>% 
  # remove outcome vars & some further vars not affecting a model we want to build
  select(-(V_bed_time_h:V_o_time_re_ptca) & -c(V_t0_hstnt_t_diff, V_delta_t_first, 
                                               V_hk_within_30_days,CPU_AdmissionDate, V_grace_low0int1high2,
                                               # V_ekg_schenkelblock  # only zeros
                                               )
         )

dat.independent <- dat.independent %>% 
  mutate(V_h_lvdys_grad = case_when(
    V_h_lvdys_grad == 0 ~ "good",
    V_h_lvdys_grad == 0.5 ~ "slightly reduced",
    V_h_lvdys_grad == 1 ~ "slightly reduced",
    V_h_lvdys_grad == 1.5 ~ "moderately reduced",
    V_h_lvdys_grad == 2 ~ "moderately reduced",
    V_h_lvdys_grad == 2.5 ~ "severely reduced",
    V_h_lvdys_grad == 3 ~ "severely reduced"
  )) %>% 
  mutate(V_h_lvdys_grad = factor(V_h_lvdys_grad, levels=c("good", "slightly reduced", "moderately reduced", "severely reduced")))


# remove data not needed for further analysis -----------------------------

rm(list=c("dat.sav", "index", "numcols", "pMiss_observe"))
