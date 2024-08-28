library(tidyverse)
library(brnn)
library(caret)
library(earth)
library(scam)
library(furrr)
library(future)
sf::sf_use_s2(FALSE)

# --------------------------------------------------------------------------#
# 1. Get our needed objects to start simulating African spread of hrp2 -----
# --------------------------------------------------------------------------#

# Get the map_0 object
isos <- na.omit(unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"]))
# MAP world map
map_0 <- readRDS("analysis/data-derived/admin0_sf.rds") %>%
  filter(iso %in% isos) %>% sf::st_as_sf()


# Get the selection model object
res_mod <- readRDS(here::here("analysis/data-derived/res_nmf_mod.rds"))

# Get the scenario maps
scenario_maps <- readRDS(here::here("analysis/data-derived/scenario_maps_full.rds"))

# Get our seed data
# read in our res from the CAR modelling
res_car <- readRDS(here::here("analysis/data-derived/CAR_modelled_res.rds"))
seeds <- pivot_wider(res_car, names_from = mut, values_from = starts_with("q"),
                     names_glue = "{mut}_{.value}")

# Get out admin 1 map with covars
afr_map <- readRDS(here::here("analysis/data-derived/map_for_CAR.rds"))
afr_map <- afr_map %>% filter(!is.na(Micro.2.10_mean))

# And remove islands

# Nearest neighbor definition
spatial_neighbors <- spdep::poly2nb(afr_map, queen=TRUE, snap=sqrt(0))
spatial_neighbors_mat <- spdep::nb2mat(spatial_neighbors, zero.policy=TRUE, style="B")

# remove islands
afr_map <- afr_map[-which(rowSums(spatial_neighbors_mat)==0), ]
seeds_adm1 <- afr_map %>% sf::st_drop_geometry() %>% left_join(seeds)

# Lastly expand our scenarios and modelling to include current resistance (low, mid, high)
scenarios <- rbind(
  scenario_maps$scenarios %>% mutate(res = "low"),
  scenario_maps$scenarios %>% mutate(res = "med"),
  scenario_maps$scenarios %>% mutate(res = "high")
)
map_data <- vector("list", nrow(scenarios))

# let's use the samples to get our AQ and LU res quantiles
# we need these so we can take samples from both posteriors
# multiple to estimate multigenotype frequencies under assumtpion
# of independence before calculating quantiles
crt_76T_prep <- readRDS("analysis/data-derived/car_models/crt_76T.rds")
mdr_86Y_prep <- readRDS("analysis/data-derived/car_models/mdr1_86Y.rds")
mdr_184F_prep <- readRDS("analysis/data-derived/car_models/mdr1_184F.rds")
mdr_CNV_prep <- readRDS("analysis/data-derived/car_models/mdr1_CNV.rds")

# frequency draws
crt76T_q <-
  do.call(rbind,crt_76T_prep$model$samples$fitted) %*%
  diag(1/crt_76T_prep$prev_df$n_test)

mdr86Y_q <-
  do.call(rbind,mdr_86Y_prep$model$samples$fitted) %*%
  diag(1/mdr_86Y_prep$prev_df$n_test)

mdr184F_q <-
  do.call(rbind,mdr_184F_prep$model$samples$fitted) %*%
  diag(1/mdr_184F_prep$prev_df$n_test)

mdrCNV_q <-
  do.call(rbind,mdr_CNV_prep$model$samples$fitted) %*%
  diag(1/mdr_CNV_prep$prev_df$n_test)

# convert into resistance frequency
aq_res_q <- (crt76T_q*mdr86Y_q) %>%
  apply(2, quantile, c(0.025,0.5,0.975))

# convert into resistance frequency
lu_res_q <- (mdr184F_q*mdrCNV_q) %>%
  apply(2, quantile, c(0.025,0.5,0.975))


# loop over and generate selection coefficients for each mutation
for(i in seq_along(map_data)) {

  micro210 <- case_when(scenarios$micro210[i] == "med" ~ seeds_adm1$Micro.2.10_mean,
                        scenarios$micro210[i] == "low" ~ seeds_adm1$Micro.2.10_low,
                        scenarios$micro210[i] == "high" ~ seeds_adm1$Micro.2.10_high)

  ft <- case_when(scenarios$ft[i] == "med" ~ seeds_adm1$ft_mean,
                  scenarios$ft[i] == "low" ~ seeds_adm1$ft_low,
                  scenarios$ft[i] == "high" ~ seeds_adm1$ft_high)

  al <- case_when(scenarios$al[i] == "med" ~ seeds_adm1$AL,
                  scenarios$al[i] == "low" ~ seeds_adm1$AL_min,
                  scenarios$al[i] == "high" ~ seeds_adm1$AL_max)

  asaq <- case_when(scenarios$asaq[i] == "med" ~ seeds_adm1$ASAQ,
                    scenarios$asaq[i] == "low" ~ seeds_adm1$ASAQ_min,
                    scenarios$asaq[i] == "high" ~ seeds_adm1$ASAQ_max)

  dhappq <- case_when(scenarios$dhappq[i] == "med" ~ seeds_adm1$DP,
                      scenarios$dhappq[i] == "low" ~ seeds_adm1$DP_min,
                      scenarios$dhappq[i] == "high" ~ seeds_adm1$DP_max)

  # normalise to sum to 1
  drugmat <- cbind(al, asaq, dhappq)
  drugmat <- apply(drugmat, 1, function(x){x/sum(x)}) %>% t()
  al <- drugmat[,1]
  asaq <- drugmat[,2]
  dhappq <- drugmat[,3]

  # get resistance for each scenario

  # K13 valid
  art_res <- case_when(scenarios$res[i] == "med" ~ seeds_adm1$k13_valid_q_50,
                       scenarios$res[i] == "low" ~ seeds_adm1$k13_valid_q_025,
                       scenarios$res[i] == "high" ~ seeds_adm1$k13_valid_q_975)

  # PPQ CNV
  ppq_res <- case_when(scenarios$res[i] == "med" ~ seeds_adm1$pfpm23_CNV_q_50,
                       scenarios$res[i] == "low" ~ seeds_adm1$pfpm23_CNV_q_025,
                       scenarios$res[i] == "high" ~ seeds_adm1$pfpm23_CNV_q_975)

  # AQ res
  aq_res <- case_when(scenarios$res[i] == "med" ~ aq_res_q[2,],
                      scenarios$res[i] == "low" ~ aq_res_q[1,],
                      scenarios$res[i] == "high" ~ aq_res_q[3,])


  # LU res
  lu_res <- case_when(scenarios$res[i] == "med" ~ lu_res_q[2,],
                      scenarios$res[i] == "low" ~ lu_res_q[1,],
                      scenarios$res[i] == "high" ~ lu_res_q[3,])

  # Loop over each allele and generate the selection coefficient
  map_data[[i]] <- map(
    paste0("s_a_",1:6),
    function(x){
      res_mod$predict(al, asaq, dhappq, art_res, ppq_res, aq_res, lu_res, ft,
                      micro210, x, f1 = 0.02, f2 = 0.10) # NB the f is not used and just needs to be provided
    }) %>%
    purrr::reduce(left_join, by = c("al","asaq","dhappq","art_res","ppq_res","aq_res","lu_res","ft","micro210"))

  map_data[[i]] <- map_data[[i]] %>%
    mutate(id_1 = seeds_adm1$id_1, .before = 1) %>%
    mutate(iso3c = seeds_adm1$iso3c, .before = 1)
}


# --------------------------------------------------------------------------#
# 2. Run our simulation -----------------------------------------------------
# --------------------------------------------------------------------------#

# Initialise our model
# get an adjacency matrix to figure out spread
adj_mat <- spdep::poly2nb(afr_map)
names(adj_mat) <- afr_map$id_1
spread_model <- arms:::R6_res_spread$new(afr_map, res_mod, adj_mat = adj_mat)

# Set up our initial conditions
seeds_med <- data.frame(
  seeds_adm1$crt_76T_q_50,
  seeds_adm1$mdr1_86Y_q_50,
  seeds_adm1$mdr1_184F_q_50,
  seeds_adm1$mdr1_CNV_q_50,
  seeds_adm1$k13_valid_q_50,
  seeds_adm1$pfpm23_CNV_q_50
) %>% setNames(paste0("a_", 1:6)) %>%
  mutate(id_1 = seeds_adm1$id_1) %>%
  split(.$id_1)

seeds_lci <- data.frame(
  seeds_adm1$crt_76T_q_025,
  seeds_adm1$mdr1_86Y_q_025,
  seeds_adm1$mdr1_184F_q_025,
  seeds_adm1$mdr1_CNV_q_025,
  seeds_adm1$k13_valid_q_025,
  seeds_adm1$pfpm23_CNV_q_025
) %>% setNames(paste0("a_", 1:6)) %>%
  mutate(id_1 = seeds_adm1$id_1) %>%
  split(.$id_1)

seeds_uci <- data.frame(
  seeds_adm1$crt_76T_q_975,
  seeds_adm1$mdr1_86Y_q_975,
  seeds_adm1$mdr1_184F_q_975,
  seeds_adm1$mdr1_CNV_q_975,
  seeds_adm1$k13_valid_q_975,
  seeds_adm1$pfpm23_CNV_q_975
) %>% setNames(paste0("a_", 1:6)) %>%
  mutate(id_1 = seeds_adm1$id_1) %>%
  split(.$id_1)

# function to run the uncertainty ranges
run_lci_med_uci_sim <- function(spread_model, md) {

  # first simulate central
  spread_model$set_map_data(md)
  med <- spread_model$simulate_multiallelic_spread(export_freq = 0.25, t_break = 1, s_name =  paste0("s_a_", 1:6), t_end = 40)
  med <- add_tf_to_output(med, spread_model)

  # then worst
  md_worst <- md %>% filter(id_1 %in% afr_map$id_1) %>%
    select(-matches("s_a_\\d$")) %>% setNames(gsub("_max", "", names(.)))

  spread_model$set_map_data(md_worst)
  uci <- spread_model$simulate_multiallelic_spread(export_freq = 0.25, t_break = 1, s_name =  paste0("s_a_", 1:6), t_end = 40)
  uci <- add_tf_to_output(uci, spread_model)

  # then best
  md_best <- md %>% filter(id_1 %in% afr_map$id_1) %>%
    select(-matches("s_a_\\d$")) %>% setNames(gsub("_min", "", names(.)))

  spread_model$set_map_data(md_best)
  lci <- spread_model$simulate_multiallelic_spread(export_freq = 0.25, t_break = 1, s_name =  paste0("s_a_", 1:6), t_end = 40)
  lci <- add_tf_to_output(lci, spread_model)

  lci %>% select(id_1, t, artR_lci = a_5, tf_lci = tf) %>%
    left_join(med %>% select(id_1, t, artR_med = a_5, tf_med = tf)) %>%
    left_join(uci %>% select(id_1, t, artR_uci = a_5, tf_uci = tf)) %>%
    select(id_1, t, matches("artR"), matches("tf"))

}

# Simulate Central
spread_model$set_seeds(seeds_med)
n_central <- which(apply(scenarios, 1, function(x){all(x == "med")}))
central_out <- run_lci_med_uci_sim(spread_model, map_data[[n_central]])

# Simulate Worst
# We use the central seeding so as to just compare the effects of epi parameter best/worst cases
# and not further complicated by uncertainty in resistance
spread_model$set_seeds(seeds_med)
n_worst <- scenarios %>% mutate(n = seq_len(n())) %>%
  filter(micro210 == "low" & ft == "high" & al == "high", res == "high") %>%
  pull(n)
worst_out <- run_lci_med_uci_sim(spread_model, map_data[[n_worst]])

# Simulate Best
# We use the central seeding so as to just compare the effects of epi parameter best/worst cases
# and not further complicated by uncertainty in resistance
spread_model$set_seeds(seeds_med)
n_best <- scenarios %>% mutate(n = seq_len(n())) %>%
  filter(micro210 == "high" & ft == "low" & al == "low", res == "low") %>%
  pull(n)
best_out <- run_lci_med_uci_sim(spread_model, map_data[[n_best]])

# ---------------------------------------------------- #
# 3. Creating Data Outputs  ----
# ---------------------------------------------------- #

meta_data <- left_join(scenario_maps$map, scenario_maps$map_data[[1]]) %>% sf::st_drop_geometry() %>%
  select(iso, name_0, name_1,id_1) %>%
mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1)

md_central_out <- central_out %>%
  left_join(meta_data %>% select(continent, name_0, name_1, id_1, iso)) %>%
  mutate(scenario = "Central") %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t, matches("artR|tf"), scenario)

md_worst_out <- worst_out %>%
  left_join(meta_data %>% select(continent, name_0, name_1, id_1, iso)) %>%
  mutate(scenario = "Worst") %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t, matches("artR|tf"), scenario)

md_best_out <- best_out %>%
  left_join(meta_data %>% select(continent, name_0, name_1, id_1, iso)) %>%
  mutate(scenario = "Best") %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t, matches("artR|tf"), scenario)

write.csv(md_central_out, "analysis/data-out/central_longitudinal_times_prospective.csv", row.names = FALSE)
write.csv(md_worst_out, "analysis/data-out/pessimistic_longitudinal_times_prospective.csv", row.names = FALSE)
write.csv(md_best_out, "analysis/data-out/optimistic_longitudinal_times_prospective.csv", row.names = FALSE)


# -----------------------

md_central_out <- read.csv("analysis/data-out/central_longitudinal_times_prospective.csv")
md_worst_out <- read.csv("analysis/data-out/pessimistic_longitudinal_times_prospective.csv")
md_best_out <- read.csv("analysis/data-out/optimistic_longitudinal_times_prospective.csv")

pop <- squire::population %>% group_by(iso3c) %>%
  summarise(n = sum(n)) %>%
  rename(iso = iso3c)

tf_med <- left_join(md_central_out, pop) %>% group_by(t) %>%
  summarise(tf_med = weighted.mean(tf_med, n))
tf_low <- left_join(md_best_out, pop) %>% group_by(t) %>%
  summarise(tf_low = weighted.mean(tf_lci, n))
tf_high <- left_join(md_worst_out, pop) %>% group_by(t) %>%
  summarise(tf_high = weighted.mean(tf_uci, n))

left_join(tf_low, tf_med) %>% left_join(tf_high) %>%
  ggplot(aes(x = t, y = tf_med, ymin = tf_low, ymax = tf_high), color = "A") +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  theme_bw() +
  MetBrewer::scale_color_met_d("Egypt", direction = -1, name = "Treatment Failure") +
  MetBrewer::scale_fill_met_d("Egypt", direction = -1, name = "Treatment Failure") +
  ylab("Treatment Failure (%)") +
  xlab("Year") +
  scale_y_continuous(labels = scales::percent)




