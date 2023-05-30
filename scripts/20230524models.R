
# used older data from EG
# new script uses data from 2023 from MY

### INFO ----------------------------------------------------------------------
# this script is sourced from `scripts/render_param_reports.R`
# scripts saves tuned models (not finalized) to: "./output/tuning_results/"

# dependencies ---------------------------------------------------------------
library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(skimr)
library(tableone)
library(pROC)
library(tidymodels)
options(tidymodels.dark = TRUE)
library(discrim)
library(finetune)
library(vetiver)
library(workflowsets)
library(baguette)
library(rules)
tidymodels_prefer()
conflicted::conflict_prefer("expand", "tidyr")

cores <- parallel::detectCores()
if (!grepl("mingw32", R.Version()$platform)) {
  library(doMC)
  registerDoMC(cores = cores/2)
} else {
  library(doParallel)
  cl <- makePSOCKcluster(cores/2)
  registerDoParallel(cl)
}


tune_grid_eval <- FALSE  # otherwise only racing methods will be used for hyperpar tuning

###
# load data ---------------------------------------------------------------

# MIRNA DAT
model_data1 <- readRDS(file = './output/Rdata/cleaned-dat/datobserve.rds') %>% 
  select(
    # ID
    V_rapID,
    # outcome
    V_o_mortality, 
    # predictors
    c("V_age", "V_sex_f1_m0", "V_symptombeginn", "V_crowding_low0_norm1_high2", "V_grace_score",
      "V_h_diabetes", "V_h_hypertonie", "V_h_cholesterin", "V_aktiver_raucher", 
      "V_h_familienana", "V_h_khk", "V_h_lvdys_grad", "V_vit_rr_syst", "V_vit_herzfrequenz", "V_vit_temperatur", 
      "V_vit_saettigung", "V_vit_atemfrequenz", "V_ekg_sinus_normal", "V_ekg_st_senkung", "V_ekg_schrittmacher", "V_ekg_schenkelblock",
      # labs
      "V_t0_hstnt_value", "V_delta_t_first", "V_t0_krea_value", "V_t0_crp_value", "V_t0_leuko_value", "V_t0_hst_value", 
      "V_t0_na_value", "V_t0_k_value", "V_t0_gluc_value", "V_t0_ck_value", 
      "V_t0_ldh_value", "V_t0_got_value", "V_t0_gpt_value", "V_t0_hb_value", 
      "V_t0_hkt_value", "V_t0_thrombo_value", "V_t0_quick_value", "V_t0_inr_value"
      )
  )

fctrs <- c("V_o_mortality", "V_sex_f1_m0", "V_symptombeginn", "V_crowding_low0_norm1_high2",
           "V_h_diabetes", "V_h_hypertonie", "V_h_cholesterin", "V_aktiver_raucher", 
           "V_h_familienana", "V_h_khk", "V_h_lvdys_grad",
           "V_ekg_sinus_normal", "V_ekg_st_senkung", "V_ekg_schrittmacher", "V_ekg_schenkelblock")
           
model_data1 <- model_data1 %>% 
  mutate(across(all_of(fctrs), as.factor)) %>% 
  mutate(V_o_mortality = factor(V_o_mortality, labels = c("survived", "died")))

pMiss <- function(x){sum(is.na(x))/length(x)*100}
pMiss_observe <- apply(model_data1,2,pMiss)
index <- pMiss_observe > 40
# drop cols with more than 40% missings
model_data1 <- model_data1[ , !index]

###
# MODEL -------------------------------------------------------------------
###
  
### Initial Split
set.seed(123)

# make control first factor for all analyses ;)
modeldat <- model_data1
  
dat_split <- rsample::initial_split(modeldat, strata = V_o_mortality)
dat_train <- training(dat_split)
dat_test <- testing(dat_split)
  
folds <- 
  vfold_cv(dat_train, strata = V_o_mortality, v = 10)
  
###
# recipe -------------------------------------------------
# A
normalized_rec <- 
  recipe(V_o_mortality ~ ., data = dat_train) %>%
  ### https://recipes.tidymodels.org/articles/Ordering.html
  #  To make sure we don’t get any unexpected results, it’s best to use 
  #  the following ordering of high-level transformations:
  #     Skewness Transformations - step_YeoJohnson()
  #     Centering, Scaling, or Normalization on Numeric Predictors
  #     Dummy Variables for Categorical Data
  update_role(V_rapID, new_role="ID") %>% 
  step_zv(all_predictors()) %>%  
  step_impute_mean(all_numeric_predictors()) %>% 
  step_impute_mode(all_nominal_predictors(), -V_o_mortality) %>% 
  step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
  step_YeoJohnson() %>% 
  step_normalize(all_numeric_predictors()) %>%  
  step_dummy(all_nominal_predictors(),-V_o_mortality)
# B
poly_rec <- 
  recipe(V_o_mortality ~ ., data = dat_train) %>%
  update_role(V_rapID, new_role="ID") %>% 
  step_zv(all_predictors()) %>%  
  step_impute_mean(all_numeric_predictors()) %>% 
  step_impute_mode(all_nominal_predictors(), -V_o_mortality) %>% 
  step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
  step_normalize(all_numeric_predictors()) %>%  
  step_poly(all_numeric_predictors()) %>% 
  step_dummy(all_nominal_predictors(),-V_o_mortality) #%>% 
#step_interact( ~all_predictors():all_predictors())
  
# C
simple_rec <- 
  recipe(V_o_mortality ~ ., data = dat_train) %>%
  update_role(V_rapID, new_role="ID") %>% 
  # ZERO VARIANCE
  step_zv(all_predictors()) %>%  
  # IMPUTE
  step_impute_mean(all_numeric_predictors()) %>% 
  step_impute_mode(all_nominal_predictors(), -V_o_mortality) %>% 
  # DECORRELATE
  step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
  step_dummy(all_nominal_predictors(),-V_o_mortality)
  
###
# specs-parsnip ----------------------------------------------------------
#parsnip::set_dependency("kknn", "glmnet", "ranger", "naivebayes") #, "kernlab", "xgboost", "nnet")
  
nearest_neighbor_kknn_spec <-
  nearest_neighbor(neighbors = tune(), weight_func = tune(), # tune() is here a placeholder!
                   dist_power = tune()) %>%
  set_engine('kknn') %>%
  set_mode('classification')
  
logistic_reg_glmnet_spec <-
  logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_mode("classification") %>% 
  set_engine('glmnet')
  
rand_forest_ranger_spec <-
  rand_forest(mtry = tune(), min_n = tune()) %>%
  set_engine('ranger') %>%
  set_mode('classification')
  
naive_Bayes_naivebayes_spec <-
  naive_Bayes(smoothness = tune(), Laplace = tune()) %>%
  set_engine('naivebayes') %>% 
  set_mode('classification')
  
svm_linear_kernlab_spec <-
  svm_linear(cost = tune(), margin = tune()) %>%
  set_engine('kernlab') %>%
  set_mode('classification')
  
svm_poly_kernlab_spec <-
  svm_poly(cost = tune(), degree = tune(), 
           scale_factor = tune(), margin = tune()) %>%
  set_engine('kernlab') %>%
  set_mode('classification')
  
svm_rbf_kernlab_spec <-
  svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>%
  set_engine('kernlab') %>%
  set_mode('classification')
  
boost_tree_xgboost_spec <-  # https://parsnip.tidymodels.org/reference/details_boost_tree_xgboost.html
  boost_tree(trees = tune(), tree_depth = tune(), 
             min_n = tune(),                              ## first three: model complexity
             loss_reduction = tune(), 
             stop_iter = tune(), 
             sample_size = tune(), mtry = tune(),         ## randomness
             learn_rate = tune(),                         ## step size
  ) %>%
  set_engine('xgboost') %>%
  set_mode('classification')
  
mlp_nnet_spec <-
  mlp(hidden_units = tune(), penalty = tune(), 
      epochs = tune()) %>%
  set_engine('nnet') %>%  #`nnet`, `keras` does not easily run on M1 Mac
  set_mode('classification')
  
##
# workflow set --------------------------------------------------------------
## A
normalized <- 
  workflow_set(   # IDs are generated {name of preproc} + {name of model}
    preproc = list(normalized = normalized_rec),
    models = list(
      KNN = nearest_neighbor_kknn_spec, 
      SVM_radial = svm_rbf_kernlab_spec, 
      SVM_poly = svm_poly_kernlab_spec, 
      SVM_linear = svm_linear_kernlab_spec,
      neural_network = mlp_nnet_spec
    )
  )
  
## B
simple <- 
  workflow_set(
    preproc = list(simple = simple_rec),
    models = list(
      #naive_bayes = naive_Bayes_naivebayes_spec,    ## ERROR when tuning "Error in pkgs$pkg[[1]] : subscript out of bounds"
      RF = rand_forest_ranger_spec,
      XGB = boost_tree_xgboost_spec
    )
  )
  
## C
with_features <- 
  workflow_set(
    preproc = list(full_quad = poly_rec),
    models = list(
      logistic_reg = logistic_reg_glmnet_spec,
      KNN = nearest_neighbor_kknn_spec
    )
  )
  
# normalized %>% extract_workflow(id = "normalized_KNN")
  
## these objects are tibbles -> row binding does not affect the state of the sets
### result is itself a workflow set
all_workflows <- 
  bind_rows(simple, normalized, with_features) %>% 
  # make workflow IDs a little more simple:
  mutate(wflow_id = gsub("(simple_)|(normalized_)", "", wflow_id))
  
###
# define-grids --------------------------------------------------------------
###
# need to finalize mtry - data dependent (unknown values - how many predictors?)
## data dependent: mtry(), sample_size(), num_terms(), num_comp()

grid_RF <- rand_forest_ranger_spec %>%   # 2 hyperparams
  extract_parameter_set_dials() %>% 
  # data dependent
  update(mtry = mtry(range = c(1, ncol(dat_train))) ) %>% 
  grid_latin_hypercube(size=100) 
  
grid_XGB <- boost_tree_xgboost_spec %>%   # 2 hyperparams
  extract_parameter_set_dials() %>% 
  # data dependent
  update(mtry = mtry(range = c(1, ncol(dat_train))) ) %>% 
  grid_latin_hypercube(size=150) 
  
grid_KNN <- nearest_neighbor_kknn_spec %>%  # 3 hyperparams
  extract_parameter_set_dials() %>%
  grid_latin_hypercube(size=100)
  
grid_SVM_radial <- svm_rbf_kernlab_spec %>%   # 3 hyperparams
  extract_parameter_set_dials() %>%
  grid_latin_hypercube(size=100)
  
grid_SVM_poly <- svm_poly_kernlab_spec %>%   # 4 hyperparams
  extract_parameter_set_dials() %>%
  grid_latin_hypercube(size=100)
  
grid_SVM_linear <- svm_linear_kernlab_spec %>%   # 2 hyperparams
  extract_parameter_set_dials() %>%
  grid_latin_hypercube(size=100)
  
grid_neural_network <- mlp_nnet_spec %>%   # 3 hyperparams
  extract_parameter_set_dials() %>% 
  update(epochs = epochs() %>% range_set(c(10, 100))) %>%   # epochs()  Range: [10, 1000] (default)
  grid_latin_hypercube(size=100)
  
grid_full_quad_logistic_reg <- logistic_reg_glmnet_spec %>%  # 2 hyperparams
  extract_parameter_set_dials() %>% 
  grid_latin_hypercube(size=100)
  
## supply grid to workflow options
# https://github.com/tidymodels/workflowsets/issues/37
all_workflows <- all_workflows %>% 
  option_add(grid = grid_RF, id = "RF") %>% 
  option_add(grid = grid_XGB, id = "XGB") %>% 
  option_add(grid = grid_KNN, id = "KNN") %>% 
  option_add(grid = grid_SVM_radial, id = "SVM_radial") %>% 
  option_add(grid = grid_SVM_poly, id = "SVM_poly") %>% 
  option_add(grid = grid_SVM_linear, id = "SVM_linear") %>% 
  option_add(grid = grid_neural_network, id = "neural_network") %>% 
  option_add(grid = grid_full_quad_logistic_reg, id = "full_quad_logistic_reg") %>% 
  option_add(grid = grid_KNN, id = "full_quad_KNN")   # same hyperparams
  
###
# tune-grid-workflows -----------------------------------------
###
  
if (tune_grid_eval == TRUE) {   # otherwise only racing methods will be used for hyperpar tuning
  set.seed(123)
  # the workflow_map() function will apply the same function to all of the workflows in the set
  ## the default is fn="tune_grid"
  grid_ctrl <- 
    control_grid(
      save_pred = TRUE,
      pkgs = NULL,
      event_level = "first",  # default
      parallel_over = "everything",
      save_workflow = TRUE
    )
  
  time1 <- Sys.time()
  grid_results <- 
    all_workflows %>% 
    workflow_map(
      # options to `tune_grid()`
      resamples = folds,
      #grid = 25,  # grid specifed above
      metrics = metric_set(roc_auc),
      control = grid_ctrl,
      # options to `workflow_map()`
      seed = 20221111,
      verbose = TRUE
    )
  time2 <- Sys.time()
  time.diff.grid <- time2-time1
  
  ## SAVE
  filename_tune_grid_results <- paste0("./output/tuning_results/", Sys.Date(), "_tune_grid_results_.rds")
  saveRDS(object = grid_results, file = filename_tune_grid_results)
  
  num_grid_models <- nrow(collect_metrics(grid_results, summarize = FALSE))
  
  
  # grid_results
  grid_results %>% 
    rank_results() %>% 
    filter(.metric == "roc_auc") %>% 
    select(model, .config, roc_auc=mean, rank) -> rankings
  
  top2models_grid <- rankings$model[1:2]
  for (j in 1:length(all_workflows$info)) {
    if (all_workflows$info[[j]]$model == top2models_grid[1]) {
      topmodel1_grid <- all_workflows$wflow_id[i]
    }
    
    if (all_workflows$info[[j]]$model == top2models_grid[2]) {
      topmodel2_grid <- all_workflows$wflow_id[j]
    }
  }
  
  # workflow-sets-plot-rank}
  autoplot(
    grid_results,
    rank_metric = "roc_auc",  # <- how to order models
    metric = "roc_auc",       # <- which metric to visualize
    select_best = TRUE     # <- one point per workflow
  ) +
    geom_text(aes(y = mean - 0.05, label = wflow_id), angle = 45, hjust = 1, size =3) +
    lims(y = c(0.6, 1.0)) +
    labs(title="Title Placeholder",
         subtitle="Grid approach")+
    ylab("ROC-AUC")+
    ggthemes::theme_few()+
    theme(legend.position = "none") -> plot_tune_grid_ranking
  
  filename_plot_tune_grid_ranking <- paste0("output/plots/tune_grid/", Sys.Date(), "_tune_grid_ranking.png")
  ggsave(filename = filename_plot_tune_grid_ranking, plot = plot_tune_grid_ranking, 
         width = 14, height = 10, 
         units = "in"  # default
  )
  
  # inspect hyperparameter results for specific model
  autoplot(grid_results, id = topmodel1_grid, metric = "roc_auc")+
    ylab("ROC-AUC")+
    labs(title= "Placeholder",
         subtitle = stringr::str_to_upper(topmodel1_grid))+
    ggthemes::theme_few() -> plot_tune_grid_hyperpars_topmodel1
  
  filename_plot_tune_grid_hyperpars_topmodel1 <- paste0("output/plots/tune_grid/", Sys.Date(), "_tune_grid_hyperpars_topmodel1.png")
  ggsave(filename = filename_plot_tune_grid_hyperpars_topmodel1, plot = plot_tune_grid_hyperpars_topmodel1, 
         width = 14, height = 10, 
         units = "in"  # default
  )
  
  autoplot(grid_results, id = topmodel2_grid, metric = "roc_auc")+
    ylab("ROC-AUC")+
    labs(title="Hyperparameter Performance - Grid Approach",
         subtitle = stringr::str_to_upper(topmodel2_grid))+
    ggthemes::theme_few() -> plot_tune_grid_hyperpars_topmodel2
  
  filename_plot_tune_grid_hyperpars_topmodel2 <- paste0("output/plots/tune_grid/", Sys.Date(), "_tune_grid_hyperpars_topmodel2.png")
  ggsave(filename = filename_plot_tune_grid_hyperpars_topmodel2, plot = plot_tune_grid_hyperpars_topmodel2, 
         width = 14, height = 10, 
         units = "in"  # default
  )
}
  
###
# tune-race-anova -----------------------------------------------------------
###
  
set.seed(123)
race_ctrl <-
  control_race(
    verbose = TRUE,
    verbose_elim = FALSE,
    allow_par = TRUE,
    save_pred = TRUE,
    burn_in = 3,
    num_ties = 10,
    alpha = 0.05,
    randomize = TRUE,
    pkgs = NULL,
    event_level = "first",
    parallel_over = "everything",
    save_workflow = TRUE
  )
  
time1 <- Sys.time()
race_results <-  all_workflows %>%
  workflow_map(
    "tune_race_anova",
    # options to `tune_race_anova()`
    resamples = folds,
    # grid = 25,  # grid specified to workflow_map above
    metrics = metric_set(roc_auc),
    control = race_ctrl,
    # options to `workflow_map()`
    seed = 20221111,
    verbose = TRUE
  )
time2 <- Sys.time()
time.diff.race <- time2-time1
  
## SAVE RACE RESULTS -----------------------------------------
filename_tune_race_results <- paste0("./output/tuning-results/", Sys.Date(), "_tune_race_results.rds")
saveRDS(object = race_results, file = filename_tune_race_results)

num_race_models <- sum(collect_metrics(race_results)$n)
  
# PLOT results ----------------------------------------------
race_results %>% 
  rank_results() %>% 
  filter(.metric == "roc_auc") %>% 
  select(model, .config, roc_auc=mean, rank) -> rankings_race
  
top2models_race <- rankings_race$model[1:2]
for (j in 1:length(all_workflows$info)) {
  if (all_workflows$info[[j]]$model == top2models_race[1]) {
    topmodel1_race <- all_workflows$wflow_id[j]
  }
  
  if (all_workflows$info[[j]]$model == top2models_race[2]) {
    topmodel2_race <- all_workflows$wflow_id[j]
  }
}
  
autoplot(
  race_results,
  rank_metric = "roc_auc",  # <- how to order models
  metric = "roc_auc",       # <- which metric to visualize
  select_best = TRUE     # <- one point per workflow
) +
  geom_text(aes(y = mean - 0.05, label = wflow_id), angle = 45, hjust = 1, size =3) +
  lims(y = c(0.6, 1.0)) +
  ylab("ROC-AUC")+
  labs(title="TITLE Placeholder",
       subtitle= paste0("Racing approach | n_grids evaluated: ", num_race_models))+
  ggthemes::theme_few()+
  theme(legend.position = "none") -> plot_tune_race_ranking
  
filename_plot_tune_race_ranking <- paste0("output/plots/tune_grid/", Sys.Date(), "_tune_race_ranking.png")
ggsave(filename = filename_plot_tune_race_ranking, plot = plot_tune_race_ranking, 
       width = 14, height = 10, 
       units = "in"  # default
)
  
if( topmodel1_race != "KNN" &  topmodel1_race != "full_quad_KNN" ) {
  #       Error in `geom_point()`:
  #         ! Problem while computing aesthetics.
  #       ℹ Error occurred in the 1st layer.
  #       Caused by error in `FUN()`:
  #         ! Objekt 'resamples' nicht gefunden
  # inspect hyperparameter results for specific model
  autoplot(race_results, id = topmodel1_race, metric = "roc_auc")+
    ylab("ROC-AUC")+
    labs(title=paste0("Hyperparameter Performance | Racing Approach:"),
         subtitle = paste0("MODEL: ", stringr::str_to_upper(topmodel1_race))) +
    ggthemes::theme_few() -> plot_tune_race_hyperpars_topmodel1
  
  filename_plot_tune_race_hyperpars_topmodel1 <- paste0("output/plots/tune_grid/", Sys.Date(), "_tune_race_hyperpars_topmodel1.png")
  ggsave(filename = filename_plot_tune_race_hyperpars_topmodel1, plot = plot_tune_race_hyperpars_topmodel1, 
         width = 14, height = 10, 
         units = "in"  # default
  )
}
  
if( topmodel2_race != "KNN" &  topmodel2_race != "full_quad_KNN" ) {
  autoplot(race_results, id = topmodel2_race, metric = "roc_auc")+
    ylab("ROC-AUC")+
    labs(title=paste0("Hyperparameter Performance | Racing Approach:"),
         subtitle = paste0("MODEL: ", stringr::str_to_upper(topmodel2_race))) +
    ggthemes::theme_few() -> plot_tune_race_hyperpars_topmodel2
  
  filename_plot_tune_race_hyperpars_topmodel2 <- paste0("output/plots/tune_grid/", Sys.Date(), "_tune_race_hyperpars_topmodel2.png")
  ggsave(filename = filename_plot_tune_race_hyperpars_topmodel2, plot = plot_tune_race_hyperpars_topmodel2, 
         width = 14, height = 10, 
         units = "in"  # default
  )
}
  
if (tune_grid_eval == TRUE) {
  matched_results <- 
    rank_results(race_results, select_best = TRUE) %>% 
    select(wflow_id, .metric, race = mean, config_race = .config) %>% 
    inner_join(
      rank_results(grid_results, select_best = TRUE) %>% 
        select(wflow_id, .metric, complete = mean, 
               config_complete = .config, model),
      by = c("wflow_id", ".metric"),
    ) %>%  
    filter(.metric == "roc_auc")
  
  matched_results %>% 
    ggplot(aes(x = complete, y = race)) + 
    geom_abline(lty = 3) + 
    geom_point() + 
    geom_text_repel(aes(label = model)) +
    coord_obs_pred() + 
    labs(x = "Complete Grid ROC-AUC", y = "Racing ROC-AUC") +
    ggthemes::theme_few() -> plot_matched_grid_race
  
  filename_plot_matched_grid_race <- paste0("output/plots/tune_grid/", Sys.Date(), "_matched_grid_race.png")
  ggsave(filename = filename_plot_matched_grid_race, plot = plot_matched_grid_race, 
         width = 14, height = 10, 
         units = "in"  # default
  )
}

## --> finalized workflows in `main.Rmd`
# last fit saved to: "./output/test_set_results/"

