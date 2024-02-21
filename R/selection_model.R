#' @title R6 Class for antimalarial resistance model
#'
#' @description Antimalarial resistance model emulator
#'
#' @importFrom R6 R6Class
#'
R6_res_mod <- R6::R6Class(
  classname = "res_mod",
  cloneable = FALSE,

  # PUBLIC METHODS
  public = list(

    # INITIALISATION
    #' @description
    #' Create a new resistance selection model using emulators.
    #' @param models Model for predicting selection coefficient
    #' @param err_models Model for predicting error of selection predictions
    #' @param model_weights Model for predicting selection coefficient weights
    #' @param err_model_weights Model for predicting selection error weights
    #' @param model_predict_f Functions to make prediction for specific model
    #' @param err_model_predict_f Functions to make prediction for specific error model
    #' @param data Data used originally for model training
    #' @param test_data Data used originally for model training
    #' @return A new `res_mod` object.
    initialize = function(models = list(),
                          err_models  = list(),
                          model_weights = list(),
                          err_model_weights = list(),
                          model_predict_f = list(),
                          err_model_predict_f = list(),
                          data = data.frame(),
                          test_data = data.frame()) {

      private$models <- models
      private$err_models <- err_models
      private$model_weights <- model_weights
      private$err_model_weights <- err_model_weights
      private$model_predict_f <- model_predict_f
      private$err_model_predict_f <- err_model_predict_f
      private$data <- data
      private$test_data <- test_data

    },

    #' @description
    #' Predict median and 95% selection coefficients and threshold times
    #' @param al AL drug use
    #' @param asaq ASAQ drug use
    #' @param dhappq DHAPPQ drug use
    #' @param art_res Art-R prevalence
    #' @param ppq_res PPQ-R prevalence
    #' @param aq_res AQ-R prevalence
    #' @param lu_res LU-R prevalence
    #' @param ft Treatment Coverage
    #' @param micro210 Microscopy 2-10 prevalence
    #' @param s_name Name of the selection coefficient to be predicted
    #' @param f1 Starting frequency of mutation
    #' @param f2 Ending frequency of mutation
    #' @param sd_scale Scaling factor for sd based estimation of CI.
    #' @return Data.frame of selection coefficients and times from 1% to 5%
    predict = function(al, asaq, dhappq, art_res, ppq_res, aq_res, lu_res, ft, micro210, s_name, f1 = 0.01, f2 = 0.05, sd_scale = 3.3) {

      # get args and turn into matrix
      dat <- cbind(al, asaq, dhappq, art_res, ppq_res, aq_res, lu_res, ft, micro210)
      ret <- as.data.frame(dat)

      # predict s
      ret[[s_name]] <- self$predict_s(dat, s_name)
      err <- self$predict_err(dat, s_name)

      # We estimated error using only 10 stochastic realisations. If we had used
      # 100 realisations, bsaed on a simple model of normally distributed variance
      # it would be 3.3 times smaller. Consequently, divide by 3.3.
      ret[[paste0(s_name, "_min")]] <- ret[[s_name]] - 1.96*(err/sd_scale)
      ret[[paste0(s_name, "_max")]] <- ret[[s_name]] + 1.96*(err/sd_scale)

      # use to calculate times
      ret[[paste0("t", "_", s_name)]] <- (log(f2 / (1 - f2)) - log(f1 / (1 - f1))) / ret[[s_name]]
      ret[[paste0("t", "_", s_name, "_max")]] <- (log(f2 / (1 - f2)) - log(f1 / (1 - f1))) / ret[[paste0(s_name, "_min")]]
      ret[[paste0("t", "_", s_name, "_min")]] <- (log(f2 / (1 - f2)) - log(f1 / (1 - f1))) / ret[[paste0(s_name, "_max")]]

      return(ret)

    },

    # Predict Ensemble s
    #' Predict selection coefficient for data frame of covariates
    #' @param dat Data frame of covariates
    #' @param s_name Name of the selection coefficient to be predicted
    predict_s = function(dat, s_name) {
      private$predict_internal(dat,
                               private$models,
                               private$model_weights,
                               private$model_predict_f,
                               s_name)
    },

    # Predict Ensemble Error
    #' Predict error in selection coefficient estimates for covariate data frame
    #' @param dat Data frame of covariates
    #' @param s_name Name of the selection coefficient to be predicted
    predict_err = function(dat, s_name) {
      exp(private$predict_error_internal(dat,
                               private$err_models,
                               private$err_model_weights,
                               private$err_model_predict_f,
                               s_name))
    },

    # Predict t
    #' Predict t between f1 and f2 and covariate data frame
    #' @param dat Data frame of covariates
    #' @param f1 Frequency at time point 1
    #' @param f2 Frequency at time point 2
    #' @param s_name Name of the relevant selection coefficient
    predict_t = function(dat, f1, f2, s_name) {

      # catch for directionality
      if(f1 >= f2) {
        stop("f1 must be less than f2")
      }

      # create results name in data frame
      # name <- paste0("t_", f1, "_", f2)

      # create s if not available
      if(!(s_name %in% names(dat))) {
        # predict s
        dat[[s_name]] <- self$predict_s(dat, s_name)
        # err <- self$predict_err(dat, s_name)
        # dat[[paste0(s_name, "_min")]] <- dat[[s_name]] - 1.96*(err/3.3)
        # dat[[paste0(s_name, "_max")]] <- dat[[s_name]] + 1.96*(err/3.3)
      }

      # create new time
      ret_t <- (log(f2 / (1 - f2)) - log(f1 / (1 - f1))) / dat[[s_name]]
      return(ret_t)
    },

    # Predict f2
    #' Predict f2 given f1 and t and covariate data frame
    #' @param dat Data frame of covariates
    #' @param f1 Frequency at time point 1
    #' @param t Duration of selection
    #' @param s_name Name of the relevant selection coefficient
    predict_f2 = function(dat, f1, t, s_name) {

      # create s if not available
      if(!(s_name %in% names(dat))) {
        # predict s
        dat[[s_name]] <- self$predict_s(dat, s_name)
        # err <- self$predict_err(dat, s_name)
        # dat[[paste0(s_name, "_min")]] <- dat[[s_name]] - 1.96*(err/3.3)
        # dat[[paste0(s_name, "_max")]] <- dat[[s_name]] + 1.96*(err/3.3)
      }

      # create new freq
      f2 = (exp(t * dat[[s_name]]) * (f1 / (1 - f1))) / (1 + exp(t * dat[[s_name]]) * (f1 / (1 - f1)))
      return(f2)
    },


    #' Add models
    #' @param model Selection prediction model
    #' @param model_name Name of model
    #' @param s_name Name of the selection coefficient the model relates to
    add_model = function(model, model_name, s_name) {
      if(is.null(private$models[[s_name]])) {
        private$models[[s_name]] <- list()
      }
      private$models[[s_name]][[model_name]] <- model
    },

    #' Add model prediction function
    #' @param f Function to be used for predictions for a specific model
    #' @param model_name Name of model
    add_model_predict_f = function(f, model_name) {
      private$model_predict_f[[model_name]] <- f
    },

    #' Add model weights
    #' @param weight Model for predicting weights (error) for each model
    #' @param model_name Name of model
    #' @param s_name Name of the selection coefficient the weights relates to
    add_model_weight = function(weight, model_name, s_name) {
      if(is.null(private$model_weights[[s_name]])) {
        private$model_weights[[s_name]] <- list()
      }
      private$model_weights[[s_name]][[model_name]] <- weight
    },

    #' Add error model
    #' @param err_model Selection error prediction model
    #' @param model_name Name of model
    #' @param s_name Name of the selection coefficient the model relates to
    add_err_model = function(err_model, model_name, s_name) {
      if(is.null(private$err_models[[s_name]])) {
        private$err_models[[s_name]] <- list()
      }
      private$err_models[[s_name]][[model_name]] <- err_model
    },

    #' Add model error prediction function
    #' @param f Function to be used for error predictions for a specific model
    #' @param model_name Name of model
    add_err_model_predict_f = function(f, model_name) {
      private$err_model_predict_f[[model_name]] <- f
    },

    #' @param weight Model for predicting weights (error) for each err model
    #' @param model_name Name of model
    #' @param s_name Name of the selection coefficient the weights relates to
    add_err_model_weight = function(weight, model_name, s_name) {
      if(is.null(private$err_model_weights[[s_name]])) {
        private$err_model_weights[[s_name]] <- list()
      }
      private$err_model_weights[[s_name]][[model_name]] <- weight
    },

    # GETTERS

    #' Get all selection prediction models
    #' @return List of models
    get_models = function() private$models,

    #' Get all selection error prediction models
    #' @return List of error models
    get_err_models = function() private$err_models,

    #' Get all model prediction functions
    #' @return List of model prediction functions
    get_model_predict_f = function() private$model_predict_f,

    #' Get all error model prediction functions
    #' @return List of error model prediction functions
    get_err_model_predict_f = function() private$err_model_predict_f,

    #' Get all selection prediction model weights
    #' @return List of mode weights
    get_model_weights = function() private$model_weights,

    #' Get all selection error prediction model weights
    #' @return List of error model weights
    get_err_model_weights = function() private$err_model_weights,

    #' Get data the models were trained on
    #' @return Training data
    get_data = function() private$data,

    #' Get the hold out test data from model training
    #' @return Training data
    get_test_data = function() private$test_data,


    # SETTERS

    #' Set all selection prediction models
    #' @param models selection model list
    set_models = function(models) { private$models <- models },

    #' Set all selection error prediction models
    #' @param err_models selection error model list
    set_err_models = function(err_models) { private$err_models <- err_models },

    #' Set all model prediction functions
    #' @param model_predict_f Model prediction function list
    set_model_predict_f = function(model_predict_f) { private$model_predict_f <- model_predict_f },

    #' Set all model error prediction functions
    #' @param err_model_predict_f Model error prediction function list
    set_err_model_predict_f = function(err_model_predict_f) { private$err_model_predict_f <- err_model_predict_f },

    #' Set all selection prediction model weights
    #' @param model_weights selection model weights (RMSE) list
    set_model_weights = function(model_weights) { private$model_weights <- model_weights },

    #' Set all selection prediction model weights
    #' @param err_model_weights selection error model weights (RMSE) list
    set_err_model_weights = function(err_model_weights) { private$err_model_weights <- err_model_weights },

    #' Set data the models were trained on
    #' @param data Training data for models
    set_data = function(data) { private$data <- data },

    #' Set hold out test data from model training
    #' @param data Training data for models
    set_test_data = function(test_data) { private$test_data <- test_data },

    #' Set s adjustment
    #' @param s_adj Adjustment for s (NMF fix)
    set_s_adj = function(s_adj) { private$s_adj <- s_adj }

  ),

  private = list(
    models = NULL,
    err_models = NULL,
    model_predict_f = NULL,
    err_model_predict_f = NULL,
    model_weights = NULL,
    err_model_weights = NULL,
    data = NULL,
    test_data = NULL,
    s_adj = 1,

    # Predict Generic
    predict_internal = function(dat, models, weights, predict_f, s_name) {

      # get our models
      model_names <- names(models[[s_name]])

      # and our err models
      model_errs <- weights[[s_name]][model_names]

      # set up our data removing NA rows
      dat <- as.data.frame(
        dat[, c("al", "asaq", "dhappq", "art_res", "ppq_res",
                "aq_res", "lu_res", "ft", "micro210")]
      )
      dat_na <- na.omit(dat)
      dat_namat <- as.matrix(dat_na)

      # make prediction for each model
      predictions <- vapply(
        lapply(model_names, function(x) {
          predict_f[[x]](models[[s_name]][[x]], dat_namat, s_name)
        }),
        as.numeric,
        FUN.VALUE = numeric(nrow(dat_namat))
      )

      # make prediction of error for each model
      errs <- vapply(
        lapply(model_names, function(x) {
          predict(model_errs[[x]], dat_namat)
        }),
        as.numeric,
        FUN.VALUE = numeric(nrow(dat_namat))
      )

      # Catch for when you are only requesting on one row of data
      if(!is.matrix(predictions)) {
        ensemb <- mean(predictions + errs)
      } else {
        ensemb <- apply(predictions + errs, MARGIN = 1, mean)
      }

      # return values with NAs in for missing data
      ret <- rep(NA, nrow(dat))
      ret[as.integer(which(apply(dat, 1, function(x){all(!is.na(x))})))] <- ensemb

      # adjust s (temp fix for NMF)
      ret <- ret * private$s_adj

      return(ret)
    },

    # Predict Generic
    predict_error_internal = function(dat, models, weights, predict_f, s_name) {

      # get our models
      model_names <- names(models[[s_name]])
      model_weights <- 1 / as.numeric(weights[[s_name]][model_names])
      normalised_weights <- model_weights / sum(model_weights)

      # set up our data removing NA rows
      dat <- as.data.frame(
        dat[, c("al", "asaq", "dhappq", "art_res", "ppq_res",
                "aq_res", "lu_res", "ft", "micro210")]
      )
      dat_na <- na.omit(dat)
      dat_namat <- as.matrix(dat_na)

      # make prediction for each model
      predictions <- vapply(
        lapply(model_names, function(x) {
          predict_f[[x]](models[[s_name]][[x]], dat_namat, s_name)
        }),
        as.numeric,
        FUN.VALUE = numeric(nrow(dat_namat))
      )

      # Catch for when you are only requesting on one row of data
      if(!is.matrix(predictions)) {
        ensemb <- weighted.mean(predictions, normalised_weights)
      } else {
        ensemb <- apply(predictions, MARGIN = 1, weighted.mean, normalised_weights)
      }

      # return values with NAs in for missing data
      ret <- rep(NA, nrow(dat))
      ret[as.integer(which(apply(dat, 1, function(x){all(!is.na(x))})))] <- ensemb

      return(ret)
    }
  )
)

#' @noRd
create_res_mod_from_res_mod <- function(res_mod){

  R6_res_mod$new(
    data = res_mod$get_data(),
    test_data = res_mod$get_test_data(),
    models = res_mod$get_models(),
    err_models = res_mod$get_err_models(),
    model_weights = res_mod$get_model_weights(),
    err_model_weights = res_mod$get_err_model_weights(),
    model_predict_f = res_mod$get_model_predict_f(),
    err_model_predict_f = res_mod$get_err_model_predict_f()
  )

}

