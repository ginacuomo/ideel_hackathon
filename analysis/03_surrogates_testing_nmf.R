# -------------------------------------------------------- #
# 0. Pull data and summary plots for model stochasticity/bias etc ----
# -------------------------------------------------------- #

# Get our data that we are training our models for
res_list <- readRDS(here::here("analysis/data-derived/modelnmf_s.rds"))
for(i in seq_along(res_list)) {
  res_list[[i]]$parset <- i
}
final <- do.call(rbind, res_list)
X <- readRDS(here::here("analysis/data-derived/sims_nmf/X.rds"))

# load caret and set up our splits
library(caret)
library(tidyverse)

# Set up training and testing data
df_med <- final %>%
  group_by(pick(al:parset)) %>%
  summarise(across(s_a_1:micro210, median)) %>%
  select(-parset)

set.seed(123)
test <- df_med %>% na.omit() %>%
  ungroup %>%
  select(-linked, -mu, -fitness, -pcr, -eir)
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# -------------------------------------------------------- #
# 0. Create known monotonic relationships between risk factors ----
# -------------------------------------------------------- #

mc_list <- list(
  "s_a_1" = c(-1, 1, 0, 1, 0, 0, -1, 1, -1),
  "s_a_2" = c(-1, 1, 0, 1, 0, 0, -1, 1, -1),
  "s_a_3" = c(1, -1, 0, 1, 0, -1, 0, 1, -1),
  "s_a_4" = c(1,  0, 0, 1, 0, -1, 0, 1, -1),
  "s_a_5" = c(0, 0, 0, 0, 1, 1, 1, 1, -1),
  "s_a_6" = c(0, 0, 1, 1, 0, 0, 0, 1, -1)
)
mc_list <- lapply(mc_list, setNames, names(train %>% select(-starts_with("s"))))

mutate_tdf <- function(tdf, s_name) {

  for(i in seq_along(mc_list[[s_name]])) {

    if(mc_list[[s_name]][i] == -1) {
      tdf[, i]  <-  1 - tdf[, i]
    }

  }

  tdf

}

# -------------------------------------------------------- #
# 1. Train model to predict hyperparameters for  xgboost ----
# -------------------------------------------------------- #

# xgboost model
library(xgboost)
library(tidyr)
library(tidyverse)
library(mlbench)
library(gbm)
library(randomForest)
library(caretEnsemble)
library(elasticnet)

# Load data
# Set up training and testing data

train_xgboost <- function(s = "s_a_5", eta = c(0.05), subsample = 0.66, colsample_bytree = c(0.66),
                          nrounds = c(200), max_depth = c(6), nfold = 5, train){

  train_control <- caret::trainControl(method="repeatedcv", number=nfold, repeats = 10)
  xgb_model <- caret::train(
    x = train %>% select(-starts_with("s")),
    y = train[[s]],
    method="xgbTree",
    trControl = train_control,
    tuneGrid = expand.grid("eta" = eta,
                           "max_depth" = max_depth,
                           "min_child_weight" = c(1),
                           "subsample" = subsample,
                           "colsample_bytree" = colsample_bytree,
                           "nrounds" = nrounds,
                           "gamma" = c(0)),
    monotone_constraints = mc_list[[s]])

  xgb_params <- list(objective = "reg:squarederror",
                     eta = xgb_model$results$eta[1],
                     max_depth = xgb_model$results$max_depth[1],
                     min_child_weight = 1,
                     subsample = xgb_model$results$subsample[1],
                     colsample_bytree = xgb_model$results$colsample_bytree[1])

  # Train xgboost model with cross-validation
  xgb_cv <- xgb.cv(params = xgb_params, data = as.matrix(train %>% select(-starts_with("s"))),
                   label = train$s_a_5, nfold = nfold, verbose = 0, na.rm = TRUE,
                   monotone_constraints = mc_list[[s]],
                   nrounds = xgb_model$results$nrounds[1])

  # Extract best iteration
  best_iter <- which.min(xgb_cv$evaluation_log$test_rmse_mean)

  # stop at best round
  xgb_model <- caret::train(
    x = train %>% select(-starts_with("s")),
    y = train[[s]],
    method="xgbTree",
    trControl = train_control,
    tuneGrid = expand.grid("eta" = xgb_model$results$eta[1],
                           "max_depth" = xgb_model$results$max_depth[1],
                           "min_child_weight" = c(1),
                           "subsample" = xgb_model$results$subsample[1],
                           "colsample_bytree" = xgb_model$results$colsample_bytree[1],
                           "nrounds" = best_iter,
                           "gamma" = c(0)),
    monotone_constraints = mc_list[[s]])

  return(xgb_model)
}

a1_xgb <- train_xgboost("s_a_1", train = train)
a2_xgb <- train_xgboost("s_a_2", train = train)
a3_xgb <- train_xgboost("s_a_3", train = train)
a4_xgb <- train_xgboost("s_a_4", train = train)
a5_xgb <- train_xgboost("s_a_5", train = train)
a6_xgb <- train_xgboost("s_a_6", train = train)
saveRDS(list(a1_xgb,a2_xgb,a3_xgb,a4_xgb,a5_xgb,a6_xgb), "analysis/data-derived/xgb_nmf_list.rds")

# -------------------------------------------------------- #
# 2. Train model to predict hyperparameters for  brnn ----
# -------------------------------------------------------- #

train_brnn <- function(s = "s_a_5",neurons = 2, train = train){
  train_control <- caret::trainControl(method="repeatedcv", number=5, repeats = 10)
  caret::train(
    x = train %>% select(-starts_with("s")),
    y = train[[s]],
    method="brnn",
    trControl = train_control,
    tuneGrid = expand.grid("neurons" = neurons))
}

a1_brnn <- train_brnn("s_a_1", train = train)
a2_brnn <- train_brnn("s_a_2", train = train)
a3_brnn <- train_brnn("s_a_3", train = train)
a4_brnn <- train_brnn("s_a_4", train = train)
a5_brnn <- train_brnn("s_a_5", train = train)
a6_brnn <- train_brnn("s_a_6", train = train)
saveRDS(list(a1_brnn,a2_brnn,a3_brnn,a4_brnn,a5_brnn,a6_brnn), "analysis/data-derived/brnn_nmf_list.rds")


# -------------------------------------------------------- #
# 3. Train model to predict hyperparameters for MONMLP ----
# -------------------------------------------------------- #

train_monmlp <- function(s = "s_a_5", h2 = 3, h1 = 3, train = train){

  train_control <- caret::trainControl(method="repeatedcv", number=5, repeats = 10)

  hout <- lapply(h2, function(h2){
    caret::train(
      x = train %>% select(-starts_with("s")) %>%
        mutate_tdf(s),
      y = train[[s]],
      method="monmlp",
      trControl = train_control,
      tuneGrid = expand.grid("hidden1" = h1, "n.ensemble" = 1),
      monotone = which(mc_list[[s]] != 0),
      hidden2 = h2)
  })


  return(hout)
}


a1_monmlp <- train_monmlp("s_a_1", train = train)
a1_monmlp_best <- a1_monmlp[[which.min(vapply(a1_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a2_monmlp <- train_monmlp("s_a_2", train = train)
a2_monmlp_best <- a2_monmlp[[which.min(vapply(a2_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a3_monmlp <- train_monmlp("s_a_3", train = train)
a3_monmlp_best <- a3_monmlp[[which.min(vapply(a3_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a4_monmlp <- train_monmlp("s_a_4", train = train)
a4_monmlp_best <- a4_monmlp[[which.min(vapply(a4_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a5_monmlp <- train_monmlp("s_a_5", train = train)
a5_monmlp_best <- a5_monmlp[[which.min(vapply(a5_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a6_monmlp <- train_monmlp("s_a_6", train = train)
a6_monmlp_best <- a6_monmlp[[which.min(vapply(a6_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]


saveRDS(list(a1_monmlp,a2_monmlp,a3_monmlp,a4_monmlp,a5_monmlp,a6_monmlp), "analysis/data-derived/monmlp_nmf_list.rds")

# -------------------------------------------------------- #
# 4. Create selection model object ----
# -------------------------------------------------------- #

# Create the res selection predicition model
res_mod <- R6_res_mod$new(data = train)

# create our model lists
xgb_mod_list <- readRDS("analysis/data-derived/xgb_nmf_list.rds")
brnn_mod_list <- readRDS("analysis/data-derived/brnn_nmf_list.rds")
monmlp_mod_list <- readRDS("analysis/data-derived/monmlp_nmf_list.rds")
# pdp comparisons show this is the best smoooth relationship and no overfitting
monmlp_mod_list <- lapply(monmlp_mod_list, function(x) {
  x[[1]]
})

names(xgb_mod_list) <- paste0("s_a_",1:6)
names(monmlp_mod_list) <- paste0("s_a_",1:6)
names(brnn_mod_list) <- paste0("s_a_",1:6)

train_brnn_pred_err <- function(test_in, err){

  train_control <- caret::trainControl(method="repeatedcv", number=5, repeats = 10)
  test_in$err <- err

  model <- caret::train(
    x = test_in %>% select(-starts_with("s")) %>% select(-err),
    y = test_in$err,
    method="brnn",
    trControl = train_control,
    tuneGrid = expand.grid("neurons" = 2))

  return(model)
}


# add all our xgb models
for(n in names(xgb_mod_list)) {
  res_mod$add_model(xgb_mod_list[[n]], "xgb", n)
  err <- (train[[n]] - predict(xgb_mod_list[[n]], train %>% select(-starts_with("s")) %>% as.matrix))
  err_mod <- train_brnn_pred_err(train, err)
  res_mod$add_model_weight(err_mod, "xgb", n)
  predict_xgb <- function(object, data, s) {
    predict(object, data)
  }
  res_mod$add_model_predict_f(predict_xgb, "xgb")

}

# add all our brnn models
for(n in names(brnn_mod_list)) {
  res_mod$add_model(brnn_mod_list[[n]], "brnn", n)
  err <- (train[[n]] - predict(brnn_mod_list[[n]], train %>% select(-starts_with("s")) %>% as.matrix))
  err_mod <- train_brnn_pred_err(train, err)
  res_mod$add_model_weight(err_mod, "brnn", n)
  predict_brnn <- function(object, data, s) {
    predict(object, data)
  }
  res_mod$add_model_predict_f(predict_brnn, "brnn")

}

# add all our monmlp models
for(n in names(monmlp_mod_list)) {
  res_mod$add_model(monmlp_mod_list[[n]], "monmlp", n)
  pred <- predict(monmlp_mod_list[[n]], train %>% select(-starts_with("s")) %>% mutate_tdf(n) %>% as.matrix)
  err <- (train[[n]] - pred)
  err_mod <- train_brnn_pred_err(train, err)
  res_mod$add_model_weight(err_mod, "monmlp", n)
  predict_monlp <- function(object, data, s) {
    mc_list <- list(
      "s_a_1" = c(-1, 1, 0, 1, 0, 0, -1, 1, -1),
      "s_a_2" = c(-1, 1, 0, 1, 0, 0, -1, 1, -1),
      "s_a_3" = c(1, -1, 0, 1, 0, -1, 0, 1, -1),
      "s_a_4" = c(1,  0, 0, 1, 0, -1, 0, 1, -1),
      "s_a_5" = c(0, 0, 0, 0, 1, 1, 1, 1, -1),
      "s_a_6" = c(0, 0, 1, 1, 0, 0, 0, 1, -1)
    )

    mutate_tdf <- function(tdf, s_name) {
      for(i in seq_along(mc_list[[s_name]])) {
        if(mc_list[[s_name]][i] == -1) {
          tdf[, i]  <-  1 - tdf[, i]
        }
      }
      tdf
    }
    predict(object, mutate_tdf(data, s))
  }
  res_mod$add_model_predict_f(predict_monlp, "monmlp")

}

saveRDS(res_mod, "analysis/data-derived/res_nmf_mod.rds")


# -------------------------------------------------------- #
# 5. Create selection error models ----
# -------------------------------------------------------- #

res_mod <- readRDS(here::here("analysis/data-derived/res_nmf_mod.rds"))

# create our data set
ps_through <- final %>%
  group_by(pick(al:parset)) %>%
  summarise(across(s_a_1:micro210, median)) %>% na.omit() %>%
  ungroup %>%
  select(-linked, -mu, -fitness, -pcr, -eir) %>%
  pull(parset)

train_reps <- do.call(rbind,res_list[ps_through][train_indices]) %>%
  group_by(parset) %>%
  mutate(micro210 = median(micro210)) %>%
  ungroup %>%
  group_by(pick(al:parset, micro210)) %>%
  summarise(across(s_a_1:s_a_6, function(x){log(sd(abs(x - median(x))))})) %>%
  ungroup %>%
  select(-linked, -mu, -fitness, -parset, -eir)

# train xgb_err_models
# simpler framework here as these don't overfit to the same degree
train_xgboost_err <- function(s = "s_a_5", eta = c(0.05,0.1,0.2), subsample = 1, colsample_bytree = c(1),
                              nrounds = c(200), max_depth = c(20), train){

  # only micro210 helps here
  mc <- c("al"=0,"asaq"=0,"dhapp"=0,"art_res"=0,"ppq_res"=0,"aq_res"=0,"lu_res"=0,"ft"=0,"micro210"=-1)

  train_control <- caret::trainControl(method="repeatedcv", number=5, repeats = 10)
  xgb_model <- caret::train(
    x = train %>% select(-starts_with("s")),
    y = train[[s]],
    method="xgbTree",
    trControl = train_control,
    tuneGrid = expand.grid("eta" = eta,
                           "max_depth" = max_depth,
                           "min_child_weight" = c(1),
                           "subsample" = subsample,
                           "colsample_bytree" = colsample_bytree,
                           "nrounds" = nrounds,
                           "gamma" = c(0.2)),
    monotone_constraints = mc)

  return(xgb_model)
}

a1_xgb_err <- train_xgboost_err("s_a_1", train = train_reps, eta = 0.05)
a2_xgb_err <- train_xgboost_err("s_a_2", train = train_reps, eta = 0.05)
a3_xgb_err <- train_xgboost_err("s_a_3", train = train_reps, eta = 0.05)
a4_xgb_err <- train_xgboost_err("s_a_4", train = train_reps, eta = 0.05)
a5_xgb_err <- train_xgboost_err("s_a_5", train = train_reps, eta = 0.05)
a6_xgb_err <- train_xgboost_err("s_a_6", train = train_reps, eta = 0.05)
saveRDS(list(a1_xgb_err, a2_xgb_err, a3_xgb_err, a4_xgb_err, a5_xgb_err, a6_xgb_err),
        "analysis/data-derived/xgb_nmf_err_list.rds")

# train brnn err models
a1_brnn_err <- train_brnn("s_a_1", train = train_reps)
a2_brnn_err <- train_brnn("s_a_2", train = train_reps)
a3_brnn_err <- train_brnn("s_a_3", train = train_reps)
a4_brnn_err <- train_brnn("s_a_4", train = train_reps)
a5_brnn_err <- train_brnn("s_a_5", train = train_reps)
a6_brnn_err <- train_brnn("s_a_6", train = train_reps)
saveRDS(list(a1_brnn_err, a2_brnn_err, a3_brnn_err, a4_brnn_err, a5_brnn_err, a6_brnn_err),
        "analysis/data-derived/brnn_nmf_err_list.rds")

# train monmlp err models
train_monmlp_err <- function(s = "s_a_5", h2 = 4, h1 = 4, train = train){

  train_control <- caret::trainControl(method="repeatedcv", number=5, repeats = 10)
  mc <- c("al"=0,"asaq"=0,"dhapp"=0,"art_res"=0,"ppq_res"=0,"aq_res"=0,"lu_res"=0,"ft"=0,"micro210"=-1)

  hout <- lapply(h2, function(h2){
    caret::train(
      # flip the micro210 due to monotne constraints
      x = train %>% mutate(micro210 = 1-micro210) %>%
        select(-starts_with("s")),
      y = train[[s]],
      method="monmlp",
      trControl = train_control,
      tuneGrid = expand.grid("hidden1" = h1, "n.ensemble" = 1),
      monotone = which(mc != 0),
      hidden2 = h2)
  })

  return(hout)
}

a1_monmlp_err <- train_monmlp_err("s_a_1", train = train_reps, h2 = 5, h1 = 4)
a2_monmlp_err <- train_monmlp_err("s_a_2", train = train_reps, h2 = 5, h1 = 4)
a3_monmlp_err <- train_monmlp_err("s_a_3", train = train_reps, h2 = 5, h1 = 4)
a4_monmlp_err <- train_monmlp_err("s_a_4", train = train_reps, h2 = 5, h1 = 4)
a5_monmlp_err <- train_monmlp_err("s_a_5", train = train_reps, h2 = 5, h1 = 4)
a6_monmlp_err <- train_monmlp_err("s_a_6", train = train_reps, h2 = 5, h1 = 4)
saveRDS(list(a1_monmlp_err[[1]], a2_monmlp_err[[1]], a3_monmlp_err[[1]],
             a4_monmlp_err[[1]], a5_monmlp_err[[1]], a6_monmlp_err[[1]]),
        "analysis/data-derived/monmlp_nmf_err_list.rds")

# create our model lists
xgb_err_mod_list <- readRDS("analysis/data-derived/xgb_nmf_err_list.rds")
brnn_err_mod_list <- readRDS("analysis/data-derived/brnn_nmf_err_list.rds")
monmlp_err_mod_list <- readRDS("analysis/data-derived/monmlp_nmf_err_list.rds")

names(xgb_err_mod_list) <- paste0("s_a_",1:6)
names(monmlp_err_mod_list) <- paste0("s_a_",1:6)
names(brnn_err_mod_list) <- paste0("s_a_",1:6)


# add all our xgb models
for(n in names(xgb_err_mod_list)) {
  res_mod$add_err_model(xgb_err_mod_list[[n]], "xgb", n)
  err <- (train_reps[[n]] - predict(xgb_err_mod_list[[n]], train_reps %>% select(-starts_with("s")) %>% as.matrix))
  weight <- mean(abs(err))
  res_mod$add_err_model_weight(weight, "xgb", n)
  predict_xgb <- function(object, data, s) {
    (predict(object, data))
  }
  res_mod$add_err_model_predict_f(predict_xgb, "xgb")

}

# add all our brnn models
for(n in names(brnn_err_mod_list)) {
  res_mod$add_err_model(brnn_err_mod_list[[n]], "brnn", n)
  err <- (train_reps[[n]] - predict(brnn_err_mod_list[[n]], train_reps %>% select(-starts_with("s")) %>% as.matrix))
  weight <- mean(abs(err))
  res_mod$add_err_model_weight(weight, "brnn", n)
  predict_brnn <- function(object, data, s) {
    (predict(object, data))
  }
  res_mod$add_err_model_predict_f(predict_brnn, "brnn")
}

# add all our monmlp models
for(n in names(monmlp_err_mod_list)) {
  res_mod$add_err_model(monmlp_err_mod_list[[n]], "monmlp", n)
  pred <- predict(
    monmlp_err_mod_list[[n]],
    train_reps %>% mutate(micro210 = 1 - micro210) %>% select(-starts_with("s")) %>% as.matrix
  )
  err <- (train_reps[[n]] - pred)
  weight <- mean(abs(err))
  res_mod$add_err_model_weight(weight, "monmlp", n)
  predict_err_monlp <- function(object, data, s) {
    data[,"micro210"] = 1 - data[,"micro210"]
    (predict(object, data))
  }
  res_mod$add_err_model_predict_f(predict_err_monlp, "monmlp")
}

saveRDS(res_mod, "analysis/data-derived/res_nmf_mod.rds")


# create our test data set
test_reps <- do.call(rbind,res_list[ps_through][-train_indices]) %>%
  group_by(parset) %>%
  mutate(micro210 = median(micro210)) %>%
  ungroup %>%
  group_by(pick(al:parset, micro210)) %>%
  summarise(across(s_a_1:s_a_6, function(x){log(sd(abs(x - median(x))))})) %>%
  ungroup %>%
  select(-linked, -mu, -fitness, -parset)

# -------------------------------------------------------- #
# 6. Rebuild if underlying R6 adapts ----
# -------------------------------------------------------- #
res_mod <- readRDS(here::here("analysis/data-derived/res_nmf_mod.rds"))
new_res <- create_res_mod_from_res_mod(res_mod)
saveRDS(new_res, here::here("analysis/data-derived/res_nmf_mod.rds"))
