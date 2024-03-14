
# -------------------------------------------------------- #
# 1. Surrogating Testing ----
# -------------------------------------------------------- #

library(tidyverse)
res_mod <- readRDS("analysis/data-derived/res_nmf_mod.rds")

# Create a plot of model predictions plot for SI
s <- "s_a_5"
dat <- res_mod$get_test_data() %>% select(-starts_with("s")) %>% as.matrix
pred_brnn <- res_mod$get_model_predict_f()$brnn(res_mod$get_models()[[s]]$brnn, dat, s)
pred_monmlp <- res_mod$get_model_predict_f()$monmlp(res_mod$get_models()[[s]]$monmlp, dat, s)
pred_xgb <- res_mod$get_model_predict_f()$xgb(res_mod$get_models()[[s]]$xgb, dat, s)
pred_ensemble <- res_mod$predict_s(dat, s)
pred_df <- rbind(
  data.frame("Prediction" = pred_brnn, "Model" = "BRNN", "Observed" = test[[s]]),
  data.frame("Prediction" = pred_monmlp, "Model" = "MONMLP", "Observed" = test[[s]]),
  data.frame("Prediction" = pred_xgb, "Model" = "XGB", "Observed" = test[[s]]),
  data.frame("Prediction" = pred_ensemble, "Model" = "Weighted Ensemble", "Observed" = test[[s]])
)
pred_df$Model <- factor(pred_df$Model, levels = c("XGB", "BRNN", "MONMLP","Weighted Ensemble"))

# Create comparative observation plots of OvsE
pred_gg <- pred_df %>% ggplot(aes(y = Observed, x = Prediction, color = Model)) +
  geom_point(alpha = 0.5) +
  # geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  facet_wrap(~Model, ncol = 4) +
  coord_equal() +
  MetBrewer::scale_color_met_d("Egypt") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
pred_gg
save_figs("model_prediction_performance", pred_gg, 14, 4, font_family = "Helvetica")

# Create model performance in terms of classical model performance statistics
mod_perf <- pred_df %>% group_by(Model) %>%
  summarise(RMSE = sqrt(mean((Observed - Prediction)^2)),
            MAE = mean(abs(Observed-Prediction)),
            '1-R2' = 1-summary(lm(Observed~Prediction))$r.squared) %>%
  mutate(Model =
           factor(Model, c("XGB", "BRNN", "MONMLP", "Weighted Ensemble"))
  )

# Create the plot for this
mod_perf_gg <- mod_perf %>% pivot_longer(RMSE:`1-R2`) %>%
  ggplot(aes(name, value, fill = Model, group = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  MetBrewer::scale_fill_met_d("Egypt") +
  xlab("Model Statistic") +
  ylab("Value")
save_figs("model_prediction_summary_performance", mod_perf_gg, 6, 4, font_family = "Helvetica")


# -------------------------------------------------------- #
# 2. Relationship Demonstration ----
# -------------------------------------------------------- #

# 1. Show plots of main ft drug relationships
al <- c(0.6,0.8,1)
dhappq <- 0
art_res <- 0.01
ppq_res <- 0.001
aq_res <- 0.001
lu_res <- 0.001
ft <- c(0.2,0.4, 0.8)
micro210 <- seq(0.01, 0.8, 0.01)
s_name <- "s_a_5"
parm_grid <- expand.grid(al, dhappq, art_res, ppq_res, aq_res, lu_res, ft, micro210)
colnames(parm_grid) <- c("al", "dhappq", "art_res", "ppq_res", "aq_res", "lu_res", "ft", "micro210")
parm_grid <- parm_grid %>% mutate(asaq = 1-al, .after = 1)
parm_grid$s_a_5 <- res_mod$get_model_predict_f()$monmlp(res_mod$get_models()$s_a_5$monmlp, parm_grid, s_name)

relation_gg <- parm_grid %>% ggplot(
  aes(micro210, !!sym(s_name), color = as.factor(ft), linetype = as.factor(al), group = interaction(ft, al))
) + geom_line() +
  theme_bw() +
  coord_equal(ratio = 2) +
  MetBrewer::scale_color_met_d("Egypt", direction = -1, name = "Effective \nTreatment \nCoverage") +
  scale_linetype(name = "Proportion \nof AL used") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  ylab("Selection Coefficient") +
  xlab("Slide Prevalence (2-10)")


save_figs("emulator_relation", relation_gg, 6, 4, font_family = "Helvetica")

# -------------------------------------------------------- #
# 3. PDP Testing ----
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
