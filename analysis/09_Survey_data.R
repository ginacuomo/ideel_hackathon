# ---------------------------------------------------- o
# 1. Turn resistance data into spatial points by admins ----
# ---------------------------------------------------- o

## packages
library(tidyverse)

# read our resistance data
resdfs <- readRDS(here::here("analysis/data-derived/resistance_dfs.rds"))
allres <- do.call(rbind, resdfs)

# read in our scenario map
scenario_maps <- readRDS(here::here("analysis/data-derived/scenario_maps_full.rds"))
map_1 <- sf::st_make_valid(scenario_maps$map)

# get our admin1 map from earlier
res_coord <- sf::st_as_sf(allres %>% dplyr::filter(continent == "Africa") %>% dplyr::select(lat, long),
                          coords = c("long", "lat"), crs = sf::st_crs(scenario_maps$map))

# identify where the coordinates fall
pts <- sf::st_within(res_coord, map_1)

# some points fall just outside border joins
missingpts <- which(is.na(as.integer(pts)))
pts[missingpts] <- sf::st_nearest_feature(res_coord[missingpts, ], map_1)

# assign the allocated admin region to our res database
afrres <- allres %>% dplyr::filter(continent == "Africa") %>%
  mutate(id_1 = scenario_maps$map$id_1[as.integer(pts)])

# now we need to group by admin region, i.e. id_1 and by year break
afr_res_splits <- afrres %>%
  mutate(yrbin = cut(syear, c(1975,2004,2010,2016,2021), right = FALSE, include.lowest = TRUE)) %>%
  group_by(yrbin, gene, mut, id_1) %>%
  summarise(x = sum(x), n = sum(n)) %>%
  mutate(prev = x/n) %>%
  split(.$mut) %>%
  lapply(function(x){split(x, x$yrbin)}) %>%
  unlist(recursive = FALSE)

# ---------------------------------------------------- #
# 2. Get covariate dataset from MAP rasters ----
# ---------------------------------------------------- #

# --------------------------------o
## 2.1 First get raster and average to admin level ----
# --------------------------------o

# rasters
rasters <- c("Walking-only travel time to healthcare map without access to motorized transport",
             "Global friction surface enumerating land-based travel speed with access to motorized transport for a nominal year 2019",
             "Global friction surface enumerating land-based travel walking-only speed without access to motorized transport for a nominal year 2019",
             "Global travel time to healthcare map with access to motorized transport")

travel_ft_map <- function(iso3c,
                          raster = "Global friction surface enumerating land-based travel speed with access to motorized transport for a nominal year 2019",
                          year = 2020) {

  # Get Admin 1 polygons from geoboundaries API
  sf_poly <- map %>% filter(iso == iso3c)

  # get countries data
  data <- cart::pull_cart(iso3c = iso3c, year = year,
                          vector_rasters = list(),
                          prevalence_rasters = list(),
                          spatial_limits_rasters = list(raster)
  )

  # extract
  extracted_data <- cart::unpack_cart(sf_poly, data)

  # population summarise by iso3c
  summary_data <- extracted_data %>%
    dplyr::select(id_1, pop, matches("layer")) %>%
    tidyr::unnest(cols = c(pop, matches("layer"))) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(id_1) %>%
    dplyr::summarise(across(starts_with("layer"), function(x){weighted.mean(x, pop, na.rm = TRUE)}),
                     pop_total = sum(pop))


  # joing back together
  fc <- dplyr::left_join(sf_poly, summary_data, by = "id_1")

  return(fc)

}

# get our isos
iso3cs <- unique(map_1$iso)
travel_list <- vector("list", length(iso3cs))
for(i in seq_along(travel_list)) {
  message(i)
  travel_list[[i]] <- travel_ft_map(iso3cs[i], rasters)
}
# not on Github due to size constraints as well
saveRDS(travel_list, here::here("analysis/data-derived/travel_list.rds"))
travel_list <- readRDS(here::here("analysis/data-derived/travel_list.rds"))
all_travel <- do.call(rbind, travel_list)

# --------------------------------o
## 2.2 Combine with our map ----
# --------------------------------o

new_df <- map_1 %>% sf::st_drop_geometry() %>% select(id_1) %>%
  left_join(all_travel %>% sf::st_drop_geometry() %>% select(id_1, pop_total, matches("layer")))

# Grab our covariate parameter ranges
covars <-  readRDS(here::here("analysis/data-derived/global_covariate_ranges.rds"))
new_df <- left_join(new_df, covars)

# bring in a covariate dataset at the iso3c level from Economist for imputation help
preds <- readRDS(here::here("analysis/data-raw/economist_covariates.RDS"))
preds <- preds %>% filter(date == 18264)
pred_df <- preds %>% select(iso3c, hospital_beds_per_thousand, life_expectancy, vdem_freedom_of_expression_score:polity_democracy_score, wdi_obs_lag:wdi_urban_pop_1m_cities_pct, percent_land_area_in_tropics:total_deaths_latest_per_100k,gdpppc_ppp_imf, wdi_gini_index)

# create out complete covars for the CAR model
full_df <- new_df %>% left_join(pred_df) %>% ungroup()
map_for_CAR <- left_join(map_1, full_df)

# add in additional covars
climate <- raster::getData('worldclim', var='bio', res=2.5)
seasonalprecip <- climate$bio15
seasonalprecip_inter <- raster::extract(seasonalprecip, map_for_CAR)
map_for_CAR$seasonality <- unlist(lapply(seasonalprecip_inter, mean, na.rm = TRUE))
map_for_CAR$seasonality[is.na(map_for_CAR$seasonality)] <- 59 # average for Sao Tome nearby

saveRDS(map_for_CAR, here::here("analysis/data-derived/map_for_CAR.rds"))


# ---------------------------------------------------- #
# 3. Build CAR models ----
# ---------------------------------------------------- #

## Prep data for CAR analysis by converting to df
prep_data_for_CAR <- function(res_df, map_for_CAR){

# create the map for a specific mutation/time dataset
res_sp <- left_join(map_for_CAR %>% dplyr::filter(!is.na(Micro.2.10_mean)), res_df)

# Nearest neighbor definition
spatial_neighbors <- spdep::poly2nb(res_sp, queen=TRUE, snap=sqrt(0))
spatial_neighbors_mat <- spdep::nb2mat(spatial_neighbors, zero.policy=TRUE, style="B")

# remove islands
res_sp <- res_sp[-which(rowSums(spatial_neighbors_mat)==0), ]
spatial_neighbors <- spdep::poly2nb(res_sp, queen=TRUE, snap=sqrt(0))
spatial_neighbors_mat <- spdep::nb2mat(spatial_neighbors, zero.policy=TRUE, style="B")

# create our df to fit model to
prev_df <- res_sp %>% sf::st_drop_geometry()
prev_df$n_test <- prev_df$n
prev_df$x_test <- prev_df$x
prev_df$n_test[is.na(prev_df$n_test)] <- 1
trials <- prev_df$n_test

# Now select our data that we want to model and scale it
df_mod <- prev_df %>%
  select(layer.1, Micro.2.10_mean, ft_mean,
         AL, gdpppc_ppp_imf,
         seasonality) %>%
  scale

return(list("mod_df" = df_mod, "prev_df" = prev_df, "W" = spatial_neighbors_mat))

}

# Create objects needed for each model fit
crt_76T_prep <- prep_data_for_CAR(afr_res_splits$`crt_76T.[2016,2021]`, map_for_CAR)
mdr_86Y_prep <- prep_data_for_CAR(afr_res_splits$`mdr1_86Y.[2016,2021]`, map_for_CAR)
mdr_184F_prep <- prep_data_for_CAR(afr_res_splits$`mdr1_184F.[2016,2021]`, map_for_CAR)
mdr_CNV_prep <- prep_data_for_CAR(afr_res_splits$`mdr1_CNV.[2016,2021]`, map_for_CAR)
k13_valid_prep <- prep_data_for_CAR(afr_res_splits$`k13_valid.[2016,2021]`, map_for_CAR)
pfpm23_CNV_prep <- prep_data_for_CAR(afr_res_splits$`pfpm23_CNV.[2016,2021]`, map_for_CAR)

# function to fit CAR models
fit_CAR_model <- function(prep, burnin = 10000, n.sample = 110000,
                          n.chains = 4, n.cores = 4, thin = 10,
                          prior.tau2 = c(100, 0.01), prior.var.beta = 100){

  # build CAR model
  model <- CARBayes::S.CARleroux(
    formula = x_test ~ .,
    family = "binomial",
    trials = prep$prev_df$n_test,
    W = prep$W,
    burnin = burnin,
    n.sample = n.sample, n.chains = n.chains, n.cores = n.cores,
    thin = thin,
    prior.tau2 = prior.tau2,
    prior.var.beta = rep(prior.var.beta^2, times = ncol(prep$mod_df)+1), ## times= no. of covariates + 1
    data = prep$mod_df %>% as.data.frame() %>% mutate(x_test = round(prep$prev_df$x)), )

  prep$model <- model
  return(prep)

}


# function to check the chain convergence
print_mcmc <- function(model, thin = 10){

  if(!is.list(model$samples$beta)){
    model$samples$beta <- list(model$samples$beta)
  }
  lapply(seq_along(model$samples$beta), function(x){
    model$samples$beta[[x]] %>%
      as.data.frame() %>%
      set_names(colnames(model$X)) %>%
      mutate(chain = x) %>%
      mutate(n = seq_len(n()))
  }) %>% do.call(rbind, .) %>%
    pivot_longer(-c(n, chain)) %>%
    filter(n %in% as.integer(seq(1, max(.data$n),length.out = (max(.$n))/thin))) %>%
    ggplot(aes(n,value,color = as.factor(chain))) +
    geom_line() +
    facet_wrap(~name, scales = "free_y")

}

# fit and add models to prep and check the mcmc chains for convergence
crt_76T_prep <- fit_CAR_model(crt_76T_prep)
print_mcmc(crt_76T_prep$model)

mdr_86Y_prep <- fit_CAR_model(mdr_86Y_prep)
print_mcmc(mdr_86Y_prep$model)

# this needs a longer burnin but does converge eventually
mdr_184F_prep <- fit_CAR_model(mdr_184F_prep)
print_mcmc(mdr_184F_prep$model)

mdr_CNV_prep <- fit_CAR_model(mdr_CNV_prep)
print_mcmc(mdr_CNV_prep$model)

k13_valid_prep <- fit_CAR_model(k13_valid_prep)
print_mcmc(k13_valid_prep$model)

pfpm23_CNV_prep <- fit_CAR_model(pfpm23_CNV_prep)
print_mcmc(pfpm23_CNV_prep$model)

saveRDS(crt_76T_prep, "analysis/data-derived/car_models/crt_76T.rds")
saveRDS(mdr_86Y_prep, "analysis/data-derived/car_models/mdr1_86Y.rds")
saveRDS(mdr_184F_prep, "analysis/data-derived/car_models/mdr1_184F.rds")
saveRDS(mdr_CNV_prep, "analysis/data-derived/car_models/mdr1_CNV.rds")
saveRDS(k13_valid_prep, "analysis/data-derived/car_models/k13_valid.rds")
saveRDS(pfpm23_CNV_prep, "analysis/data-derived/car_models/pfpm23_CNV.rds")

# ---------------------------------------------------- #
# 4. Calculating resistance frequency ----
# ---------------------------------------------------- #

# create summary samples or use the whole chain and summarise that
sample_res <- function(prep, n = 1000, sample = FALSE){

  if (sample) {

  ## Combine our different chains
  mod_beta <- do.call(rbind, prep$model$samples$beta)
  mod_phi <- do.call(rbind, prep$model$samples$phi)

  # Get our covariates for each admin region
  # covars_time <- cbind(rep(1, nrow(prep$prev_df)), prep$mod_df)

  # Samples from posterior
  p_k_1 <- matrix(nrow=n, ncol=nrow(covars_time))
  samps <- sample(nrow(mod_beta), n, replace = FALSE)

  # Calculate model prediction
  for (k in seq_len(nrow(covars_time))) {
    for (iter in seq_len(n)) {
      x <- covars_time[k,] %*% mod_beta[samps[iter],] + mod_phi[samps[iter], k]
      p_k_1[iter,k] <- exp(x)/(1+exp(x))
    }
  }

  } else {

    p_k_1 <- do.call(rbind, prep$model$samples$fitted)

  }

  # admin summaries
  adm_summaries <-
    apply(p_k_1, 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %*%
    diag(1/prep$prev_df$n_test) %>%
    t() %>%
    as.data.frame()

  names(adm_summaries) <- paste0("q_",c("025","25","50","75","975"))
  return(adm_summaries)

}

# create summaries for each marker
crt_76T_est <- sample_res(crt_76T_prep) %>% mutate(id_1 = crt_76T_prep$prev_df$id_1, mut = "crt_76T")
mdr1_86Y_est <- sample_res(mdr_86Y_prep) %>% mutate(id_1 = crt_76T_prep$prev_df$id_1, mut = "mdr1_86Y")
mdr1_184F_est <- sample_res(mdr_184F_prep) %>% mutate(id_1 = crt_76T_prep$prev_df$id_1, mut = "mdr1_184F")
mdr1_CNV_est <- sample_res(mdr_CNV_prep) %>% mutate(id_1 = crt_76T_prep$prev_df$id_1, mut = "mdr1_CNV")
k13_valid_est <- sample_res(k13_valid_prep) %>% mutate(id_1 = crt_76T_prep$prev_df$id_1, mut = "k13_valid")
pfpm23_CNV_est <- sample_res(pfpm23_CNV_prep) %>% mutate(id_1 = crt_76T_prep$prev_df$id_1, mut = "pfpm23_CNV")

# bind together and save
model_est_res <- rbind(crt_76T_est, mdr1_86Y_est, mdr1_184F_est, mdr1_CNV_est, k13_valid_est, pfpm23_CNV_est)
saveRDS(model_est_res, here::here("analysis/data-derived/CAR_modelled_res.rds"))
