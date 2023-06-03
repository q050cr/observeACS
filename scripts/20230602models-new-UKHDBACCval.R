
# this script uses the BACC cohort as the test cohort for the build algorithms

# INFO ----------------------------------------------------------------------
# this script is sourced from `scripts/render_param_reports.R`
# scripts saves tuned models (not finalized) to: "./output/tuning_results_UKHDBACC/"

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

run_model <- FALSE
tune_grid_eval <- FALSE  # otherwise only racing methods will be used for hyperpar tuning
imputation_knn <- TRUE
logged_numeric <- FALSE  # we did YeoJohnson step ( power transformation method used to stabilize variance and make data more closely follow a normal distribution)

# cross-validation
v <- 10
# mixed cohort (both UKHD & BACC in train/test)
mixed_cohort <- TRUE
if(mixed_cohort==TRUE){
  mixed_cohort_naming <- "traintestmixed"
}else{
  mixed_cohort_naming <- "trainUKHDtestBACC"
}

# load data ---------------------------------------------------------------

combined_dataset <- readRDS(file = "./output/Rdata/cleaned-dat/2023-06-02-combined_dataset_ukhd_bacc.rds")  # is "combined_dataset" from data_preprocessing_BACC.R

# BACC cohort has 62% missings in ntproBNP - drop this for modelling! 

###
# MODEL -------------------------------------------------------------------
###

## Initial Split ----
set.seed(123)
# make control first factor for all analyses ;)
modeldat <- combined_dataset

if(mixed_cohort == TRUE) {
  dat_split <- rsample::initial_split(modeldat, strata = o_mortality)
  dat_train <- training(dat_split)
  dat_test <- testing(dat_split)
} else{
  dat_split <- rsample::group_initial_split(modeldat, group = source)
  dat_train <- training(dat_split) 
  dat_test <- testing(dat_split)
}

folds <- 
  vfold_cv(dat_train, strata = o_mortality, v = v)

###
## Recipes -------------------------------------------------
# A
normalized_rec <- 
  recipe(o_mortality ~ ., data = dat_train) %>%
  ### https://recipes.tidymodels.org/articles/Ordering.html
  #  To make sure we don’t get any unexpected results, it’s best to use 
  #  the following ordering of high-level transformations:
  #     Skewness Transformations - step_YeoJohnson()
  #     Centering, Scaling, or Normalization on Numeric Predictors
  #     Dummy Variables for Categorical Data
  update_role(rapID, new_role="ID") %>% 
  update_role(source, new_role="ID") %>% 
  step_rm(t0_ntbnp_value, aktiver_raucher) %>%  ## too many missings in BACC, MNAR!!
  step_zv(all_predictors()) %>%  
  step_impute_mean(all_numeric_predictors()) %>% 
  step_impute_mode(all_nominal_predictors(), -o_mortality) %>% 
  step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
  step_YeoJohnson(all_numeric_predictors()) %>% 
  step_normalize(all_numeric_predictors()) %>%  
  step_dummy(all_nominal_predictors(),-o_mortality)

if (imputation_knn == TRUE) {
  normalized_rec <- 
    recipe(o_mortality ~ ., data = dat_train) %>%
    ### https://recipes.tidymodels.org/articles/Ordering.html
    #  To make sure we don’t get any unexpected results, it’s best to use 
    #  the following ordering of high-level transformations:
    #     Skewness Transformations - step_YeoJohnson()
    #     Centering, Scaling, or Normalization on Numeric Predictors
    #     Dummy Variables for Categorical Data
    update_role(rapID, new_role="ID") %>% 
    update_role(source, new_role="ID") %>% 
    step_rm(t0_ntbnp_value, aktiver_raucher) %>%  ## too many missings in BACC, MNAR!!
    step_zv(all_predictors()) %>%  
    step_impute_mean(t0_krea_value, t0_crp_value) %>% # only one missing krea value, 4 missing crps
    step_impute_knn(all_predictors(), neighbors = 5, 
                    impute_with = imp_vars(age, sex_f1_m0, t0_hstnt_value, ekg_sinus_normal, t0_krea_value, t0_crp_value)) %>% 
    
    step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
    step_YeoJohnson(all_numeric_predictors()) %>% 
    step_normalize(all_numeric_predictors()) %>%  
    step_dummy(all_nominal_predictors(),-o_mortality)
}

# B
poly_rec <- 
  recipe(o_mortality ~ ., data = dat_train) %>%
  update_role(rapID, new_role="ID") %>% 
  update_role(source, new_role="ID") %>% 
  step_rm(t0_ntbnp_value, aktiver_raucher) %>%  ## too many missings in BACC, MNAR!!
  step_zv(all_predictors()) %>%  
  step_impute_mean(all_numeric_predictors()) %>% 
  step_impute_mode(all_nominal_predictors(), -o_mortality) %>% 
  step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
  step_normalize(all_numeric_predictors()) %>%  
  step_poly(all_numeric_predictors()) %>% 
  step_dummy(all_nominal_predictors(),-o_mortality) #%>% 
#step_interact( ~all_predictors():all_predictors())

if (imputation_knn == TRUE) {
  poly_rec <- 
    recipe(o_mortality ~ ., data = dat_train) %>%
    update_role(rapID, new_role="ID") %>% 
    update_role(source, new_role="ID") %>% 
    step_rm(t0_ntbnp_value, aktiver_raucher) %>%  ## too many missings in BACC, MNAR!!
    step_zv(all_predictors()) %>%  
    step_impute_mean(t0_krea_value, t0_crp_value) %>% # only one missing krea value, 4 missing crps
    step_impute_knn(all_predictors(), neighbors = 5, 
                    impute_with = imp_vars(age, sex_f1_m0, t0_hstnt_value, ekg_sinus_normal, t0_krea_value, t0_crp_value)) %>% 
    
    step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
    step_normalize(all_numeric_predictors()) %>%  
    step_poly(all_numeric_predictors()) %>% 
    step_dummy(all_nominal_predictors(),-o_mortality) #%>% 
  #step_interact( ~all_predictors():all_predictors())
}

# C
simple_rec <- 
  recipe(o_mortality ~ ., data = dat_train) %>%
  update_role(rapID, new_role="ID") %>% 
  update_role(source, new_role="ID") %>% 
  step_rm(t0_ntbnp_value, aktiver_raucher) %>%  ## too many missings in BACC, MNAR!!
  # ZERO VARIANCE
  step_zv(all_predictors()) %>%  
  # IMPUTE
  step_impute_mean(all_numeric_predictors()) %>% 
  step_impute_mode(all_nominal_predictors(), -o_mortality) %>% 
  # DECORRELATE
  step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
  step_dummy(all_nominal_predictors(),-o_mortality)

if (imputation_knn == TRUE) {
  simple_rec <- 
    recipe(o_mortality ~ ., data = dat_train) %>%
    update_role(rapID, new_role="ID") %>% 
    update_role(source, new_role="ID") %>% 
    step_rm(t0_ntbnp_value, aktiver_raucher) %>%  ## too many missings in BACC, MNAR!!
    # ZERO VARIANCE
    step_zv(all_predictors()) %>%  
    # IMPUTE
    step_impute_mean(t0_krea_value, t0_crp_value) %>% # only one missing krea value, 4 missing crps
    step_impute_knn(all_predictors(), neighbors = 5, 
                    impute_with = imp_vars(age, sex_f1_m0, t0_hstnt_value, ekg_sinus_normal, t0_krea_value, t0_crp_value)) %>% 
    
    # DECORRELATE
    step_corr(all_numeric_predictors(), threshold = 0.8) %>% 
    step_dummy(all_nominal_predictors(),-o_mortality)
}

prepped <- poly_rec %>% prep(strings_as_factors = FALSE)  # otherwise ID columns --> NA
test <- bake(prepped, new_data = dat_test)

###
## SPECS-parsnip ----------------------------------------------------------
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
## Workflow Set --------------------------------------------------------------
## A
normalized <- 
  workflow_set(   # IDs are generated {name of preproc} + {name of model}
    preproc = list(normalized = normalized_rec),
    models = list(
      KNN = nearest_neighbor_kknn_spec, 
      SVM_radial = svm_rbf_kernlab_spec, 
      SVM_poly = svm_poly_kernlab_spec, 
      SVM_linear = svm_linear_kernlab_spec,
      neural_network = mlp_nnet_spec,
      logistic_reg_norm = logistic_reg_glmnet_spec
    )
  )

## B
simple <- 
  workflow_set(
    preproc = list(simple = simple_rec),
    models = list(
      #naive_bayes = naive_Bayes_naivebayes_spec,    ## ERROR when tuning "Error in pkgs$pkg[[1]] : subscript out of bounds"
      RF = rand_forest_ranger_spec,
      XGB = boost_tree_xgboost_spec,
      logistic_reg_simple = logistic_reg_glmnet_spec
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
## Grids --------------------------------------------------------------
###
# need to finalize mtry - data dependent (unknown values - how many predictors?)
## data dependent: mtry(), sample_size(), num_terms(), num_comp()

grid_RF <- rand_forest_ranger_spec %>%   # 2 hyperparams
  extract_parameter_set_dials() %>% 
  # data dependent
  update(mtry = mtry(range = c(1, ncol(dat_train)-2)) ) %>%   # ids removed
  grid_latin_hypercube(size=200) 

grid_XGB <- boost_tree_xgboost_spec %>%   # 2 hyperparams
  extract_parameter_set_dials() %>% 
  # data dependent
  update(mtry = mtry(range = c(1, ncol(dat_train)-2)) ) %>% 
  grid_latin_hypercube(size=250) 

grid_KNN <- nearest_neighbor_kknn_spec %>%  # 3 hyperparams
  extract_parameter_set_dials() %>%
  grid_latin_hypercube(size=200)

grid_SVM_radial <- svm_rbf_kernlab_spec %>%   # 3 hyperparams
  extract_parameter_set_dials() %>%
  grid_latin_hypercube(size=200)

grid_SVM_poly <- svm_poly_kernlab_spec %>%   # 4 hyperparams
  extract_parameter_set_dials() %>%
  grid_latin_hypercube(size=250)

grid_SVM_linear <- svm_linear_kernlab_spec %>%   # 2 hyperparams
  extract_parameter_set_dials() %>%
  grid_latin_hypercube(size=200)

grid_neural_network <- mlp_nnet_spec %>%   # 3 hyperparams
  extract_parameter_set_dials() %>% 
  update(epochs = epochs() %>% range_set(c(10, 100))) %>%   # epochs()  Range: [10, 1000] (default)
  grid_latin_hypercube(size=200)

grid_full_quad_logistic_reg <- logistic_reg_glmnet_spec %>%  # 2 hyperparams
  extract_parameter_set_dials() %>% 
  grid_latin_hypercube(size=200)

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
  option_add(grid = grid_full_quad_logistic_reg, id = "logistic_reg_simple") %>%   # same hyperparams
  option_add(grid = grid_full_quad_logistic_reg, id = "logistic_reg_norm")

###
## Tune-grid-workflows -----------------------------------------
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
  filename_tune_grid_results <- paste0("./output/tuning_results_UKHDBACC/", Sys.Date(), "_tune_grid_results_v_", v, ".rds")
  saveRDS(object = grid_results, file = filename_tune_grid_results)
  
  num_grid_models <- nrow(collect_metrics(grid_results, summarize = FALSE))
  
  
  # grid_results
  grid_results %>% 
    rank_results() %>% 
    filter(.metric == "roc_auc") %>% 
    select(model, .config, roc_auc=mean, rank) -> rankings
}

###
## tune-race-anova -----------------------------------------------------------
###
if(run_model == TRUE) {
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
} else{
  if(mixed_cohort==TRUE){
    ## LOAD RACE? ----
    race_results <- readRDS(file = "./output/tuning-results-UKHDBACC/cluster/2023-06-02_tune_race_results_knnImput_v_10_traintestmixed.rds")
  }else{
    race_results <- readRDS(file = "./output/tuning-results-UKHDBACC/cluster/2023-06-02_tune_race_results_knnImput_v_10_trainUKHDtestBACC.rds")
  }
  
}


## SAVE RACE RESULTS -----------------------------------------

if (run_model == TRUE & imputation_knn == TRUE & logged_numeric == TRUE) {
  filename_tune_race_results <- paste0("./output/tuning-results-UKHDBACC/cluster/", Sys.Date(), "_tune_race_results_knnImput_logged_num_v_", v, "_", mixed_cohort_naming, ".rds")
  saveRDS(object = race_results, file = filename_tune_race_results)
} else if (run_model == TRUE & imputation_knn == TRUE & logged_numeric == FALSE) {
  filename_tune_race_results <- paste0("./output/tuning-results-UKHDBACC/cluster/", Sys.Date(), "_tune_race_results_knnImput_v_", v,"_", mixed_cohort_naming, ".rds")
  saveRDS(object = race_results, file = filename_tune_race_results)
} else if(run_model == TRUE) {
  filename_tune_race_results <- paste0("./output/tuning-results-UKHDBACC/cluster/", Sys.Date(), "_tune_race_results_v_", v,"_", mixed_cohort_naming, ".rds")
  saveRDS(object = race_results, file = filename_tune_race_results)
}

num_race_models <- sum(collect_metrics(race_results)$n)

# PLOT results in .Rmd file
race_results %>% 
  rank_results() %>% 
  filter(.metric == "roc_auc") %>% 
  select(model, .config, roc_auc=mean, rank) -> rankings_race

## --> finalized workflows in `main.Rmd`
# last fit saved to: "./output/test_set_results/"

