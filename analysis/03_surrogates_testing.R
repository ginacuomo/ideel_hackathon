# -------------------------------------------------------- #
# 0. Pull data and summary plots for model stochasticity/bias etc ----
# -------------------------------------------------------- #

# Get our data that we are training our models for
res_list <- readRDS(here::here("analysis/data-derived/model_s.rds"))
for(i in seq_along(res_list)) {
  res_list[[i]]$parset <- i
}
final <- do.call(rbind, res_list)
X <- readRDS(here::here("analysis/data-derived/sims/X.rds"))

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

train_xgboost <- function(s = "s_a_5", eta = c(0.05,0.1,0.2), subsample = 1, colsample_bytree = c(1),
                          nrounds = c(200), max_depth = c(50)){

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
                           "gamma" = c(0)),
    monotone_constraints = mc_list[[s]])

  # Evaluate performance on test set
  # pred <- predict(xgb_model, test %>% select(-starts_with("s")) %>% as.matrix())
  # rmse <- sqrt(mean((test[[s]] - pred) ^ 2))
  # plot(pred, test[[s]], xlab = paste0("XGBOOST Model Predictions of Selection Coefficient", s),
  #      ylab = "Observed Transmission Model Selection Coefficient")
  # abline(0, 1, col = "red")

  return(xgb_model)
}

a1_xgb <- train_xgboost("s_a_1")
a2_xgb <- train_xgboost("s_a_2")
a3_xgb <- train_xgboost("s_a_3")
a4_xgb <- train_xgboost("s_a_4")
a5_xgb <- train_xgboost("s_a_5")
a6_xgb <- train_xgboost("s_a_6")
saveRDS(list(a1_xgb,a2_xgb,a3_xgb,a4_xgb,a5_xgb,a6_xgb), "analysis/data-derived/xgb_list.rds")

# -------------------------------------------------------- #
# 2. Train model to predict hyperparameters for  brnn ----
# -------------------------------------------------------- #

train_control <- caret::trainControl(method="repeatedcv", number=5, repeats = 10)
train_brnn <- function(s = "s_a_5",neurons = 2){
caret::train(
  x = train %>% select(-starts_with("s")),
  y = train[[s]],
  method="brnn",
  trControl = train_control,
  tuneGrid = expand.grid("neurons" = neurons))
}

a1_brnn <- train_brnn("s_a_1")
a2_brnn <- train_brnn("s_a_2")
a3_brnn <- train_brnn("s_a_3")
a4_brnn <- train_brnn("s_a_4")
a5_brnn <- train_brnn("s_a_5")
a6_brnn <- train_brnn("s_a_6")
saveRDS(list(a1_brnn,a2_brnn,a3_brnn,a4_brnn,a5_brnn,a6_brnn), "analysis/data-derived/brnn_list.rds")


# -------------------------------------------------------- #
# 3. Train model to predict hyperparameters for MONMLP ----
# -------------------------------------------------------- #

train_monmlp <- function(s = "s_a_5", h2 = 3:6, h1 = 4:10){

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


a1_monmlp <- train_monmlp("s_a_1")
a1_monmlp_best <- a1_monmlp[[which.min(vapply(a1_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a2_monmlp <- train_monmlp("s_a_2")
a2_monmlp_best <- a2_monmlp[[which.min(vapply(a2_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a3_monmlp <- train_monmlp("s_a_3")
a3_monmlp_best <- a3_monmlp[[which.min(vapply(a3_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a4_monmlp <- train_monmlp("s_a_4")
a4_monmlp_best <- a4_monmlp[[which.min(vapply(a4_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a5_monmlp <- train_monmlp("s_a_5")
a5_monmlp_best <- a5_monmlp[[which.min(vapply(a5_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]

a6_monmlp <- train_monmlp("s_a_6")
a6_monmlp_best <- a6_monmlp[[which.min(vapply(a6_monmlp, function(x){min(x$results$RMSE)}, numeric(1)))]]


saveRDS(list(a1_monmlp,a2_monmlp,a3_monmlp,a4_monmlp,a5_monmlp,a6_monmlp), "analysis/data-derived/monmlp_list.rds")

# -------------------------------------------------------- #
# 4. Create selection model object ----
# -------------------------------------------------------- #

# Create the res selection predicition model
res_mod <- R6_res_mod$new(data = test)

# create our model lists
xgb_mod_list <- readRDS("analysis/data-derived/xgb_list.rds")
brnn_mod_list <- readRDS("analysis/data-derived/brnn_list.rds")
monmlp_mod_list <- readRDS("analysis/data-derived/monmlp_list.rds")
# pdp comparisons show this is the best smoooth relationship and no overfitting
monmlp_mod_list <- lapply(monmlp_mod_list, function(x) {
  x[[3]]
})

names(xgb_mod_list) <- paste0("s_a_",1:6)
names(monmlp_mod_list) <- paste0("s_a_",1:6)
names(brnn_mod_list) <- paste0("s_a_",1:6)

train_xgboost_err <- function(test_in, err){

  train_control <- caret::trainControl(method="repeatedcv", number=5, repeats = 10)
  test_in$err <- err
  xgb_model <- caret::train(
    x = test_in %>% select(-starts_with("s")),
    y = test_in$err,
    method="xgbTree",
    trControl = train_control,
    tuneGrid = expand.grid("eta" = 0.03,
                           "max_depth" = 50,
                           "min_child_weight" = c(1),
                           "subsample" = 1,
                           "colsample_bytree" = 1,
                           "nrounds" = 200,
                           "gamma" = c(0)))

  # Evaluate performance on test set
  # pred <- predict(xgb_model, test %>% select(-starts_with("s")) %>% as.matrix())
  # rmse <- sqrt(mean((test[[s]] - pred) ^ 2))
  # plot(pred, test[[s]], xlab = paste0("XGBOOST Model Predictions of Selection Coefficient", s),
  #      ylab = "Observed Transmission Model Selection Coefficient")
  # abline(0, 1, col = "red")

  return(xgb_model)
}


# add all our xgb models
for(n in names(xgb_mod_list)) {
  res_mod$add_model(xgb_mod_list[[n]], "xgb", n)
  rmse <- sqrt(mean((test[[n]] - predict(xgb_mod_list[[n]], test %>% select(-starts_with("s")) %>% as.matrix)) ^ 2))

  res_mod$add_model_weight(rmse, "xgb", n)
  predict_xgb <- function(object, data, s) {
    predict(object, data)
  }
  res_mod$add_model_predict_f(predict_xgb, "xgb")
}

# # add all our brnn models
for(n in names(brnn_mod_list)) {
  res_mod$add_model(brnn_mod_list[[n]], "brnn", n)
  rmse <- sqrt(mean((test[[n]] - predict(brnn_mod_list[[n]], test %>% select(-starts_with("s")) %>% as.matrix)) ^ 2))
  res_mod$add_model_weight(rmse, "brnn", n)
  predict_brnn <- function(object, data, s) {
    predict(object, data)
  }
  res_mod$add_model_predict_f(predict_brnn, "brnn")
}

# add all our monmlp models
for(n in names(monmlp_mod_list)) {
  res_mod$add_model(monmlp_mod_list[[n]], "monmlp", n)
  pred <- predict(monmlp_mod_list[[n]], test %>% select(-starts_with("s")) %>% mutate_tdf("s_a_5") %>% as.matrix)
  rmse <- sqrt(mean((test[[n]] - pred) ^ 2))
  res_mod$add_model_weight(rmse, "monmlp", n)
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

saveRDS(res_mod, "analysis/data-derived/res_mod.rds")


# -------------------------------------------------------- #
# 5. PDP Testing ----
# -------------------------------------------------------- #

# Let's produce these for each locus

pdp_plot <- function(res_mod, s, grid.resolution = 10) {

mods <- res_mod$get_models()[[s]]

# Create our partial dependence
vars <- names(train %>% select(-starts_with("s_")))
pdp_dat <- lapply(seq_along(vars), function(var){
  dat <- do.call(rbind, lapply(seq_along(mods), function(i){
    if(names(mods)[i] == "monmlp"){
      out <- pdp::partial(mods[[i]], pred.var = as.character(vars[var]),
                   train = train %>% select(-starts_with("s")) %>% mutate_tdf(s),
                   grid.resolution = grid.resolution, type = "regression") %>%
        rename(s = yhat) %>%
        mutate(model = names(mods)[i]) %>%
        setNames(c("x", "s", "model"))
      if(mc_list[[s]][var] == -1) {
        out$x <- 1-out$x
      }
      return(out)
    } else {
    pdp::partial(mods[[i]], pred.var = as.character(vars[var]),
                 train = train %>% select(-starts_with("s")), grid.resolution = 20, type = "regression") %>%
      rename(s = yhat) %>%
      mutate(model = names(mods)[i]) %>%
      setNames(c("x", "s", "model"))
    }
  }))
  dat %>% ggplot(aes(x,s,color = model)) +
    geom_line() +
    xlab(vars[var]) +
    ylab(paste0("Selection Coeffiecient:" , s)) +
    ggpubr::theme_pubclean(base_size = 14) +
    theme(axis.line = element_line()) +
    scale_color_brewer(type = "qual") +
    theme(legend.key = element_rect(fill = "white"),
          text = element_text(family = "Helvetica"))

})

ggpdp <- cowplot::plot_grid(plotlist = pdp_dat)
return(ggpdp)

}

pdp_plot(res_mod, "s_a_5", grid.resolution = 20)


# -------------------------------------------------------- #
# 6. Create selection error models ----
# -------------------------------------------------------- #

res_mod <- readRDS(here::here("analysis/data-derived/res_mod.rds"))
test_reps <- do.call(rbind,res_list[train_indices]) %>%
  select(s_a_1:micro210, al:ft, parset) %>%
  select(-eir)

# get out models
mods <- res_mod$get_models()

# train error models
err_models <- map(seq_along(mods), function(j) {

  mod_j <- mods[[j]]

  map(seq_along(mod_j), function(i) {

    s_nm <- names(mods)[j]
    mod <- mod_j[[i]]
    mod_nm <- names(mod_j)[i]
    test_df <- test_reps %>% select(names(mod$ptype), parset, s_nm)

  # make error prediction
  predf <- res_mod$get_model_predict_f()[[mod_nm]]
  pred <- predf(mod,
                test_df %>% select(names(mod$ptype)) %>% as.matrix,
                s = s_nm)

  # work out the per simulation
  test_df$error <- pred-test_df[[s_nm]]
  test_df$error <- abs(pred-test_df[[s_nm]])
  test_df_sum <- test_df %>%
    group_by(pick(al:ft,parset)) %>%
    summarise(micro210 = mean(micro210),
              !!sym(s_nm) := median(!!sym(s_nm)),
              error = sd(error, na.rm = TRUE)) %>%
    na.omit() %>%
    ungroup()

  train_indices <- sample(nrow(test_df_sum), nrow(test_df_sum) * 0.75)
  train_df <- test_df_sum[train_indices, ]
  test_df <- test_df_sum[-train_indices, ]

  # now do standard train for a brnn 2 neuron model
  train_control <- caret::trainControl(method="cv", number=20)
  tune_grid <- expand.grid("neurons" = 2)
  err_model <- caret::train(
    x = train_df %>% select(al:ft, micro210),
    y = train_df$error, method="brnn",
    trControl = train_control, tuneGrid = tune_grid)

  # Evaluate performance on test set
  pred <- predict(err_model, test_df %>% select(al:ft, micro210))
  rmse_error <- sqrt(mean((test_df$error - pred) ^ 2))

  return(list("err_model" = err_model, "rmse" = rmse_error, "mod" = mod_nm, "s" = s_nm))

})
})

# add to our model
for(i in seq_along(err_models)) {
  for(j in seq_along(err_models[[i]])){
  res_mod$add_err_model(err_models[[i]][[j]]$err_model,
                        model_name = err_models[[i]][[j]]$mod,
                        s_name = err_models[[i]][[j]]$s)
    res_mod$add_err_model_weight(weight = err_models[[i]][[j]]$rmse,
                                 model_name = err_models[[i]][[j]]$mod,
                                 s_name = err_models[[i]][[j]]$s)
    predict_brnn <- function(object, data, s) {
      predict(object, data)
    }
    res_mod$add_err_model_predict_f(predict_brnn, err_models[[i]][[j]]$mod)
  }
}
saveRDS(res_mod, here::here("analysis/data-derived/res_mod.rds"))
