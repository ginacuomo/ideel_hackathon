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
map_0 <- malariaAtlas::getShp(ISO = isos, admin_level = c("admin0")) %>% sf::st_as_sf()

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


# loop over and generate selection and times
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

  map_data[[i]] <- res_mod$predict(al, asaq, dhappq, art_res, ppq_res, aq_res, lu_res, ft,
                                   micro210, "s_a_5", f1 = 0.01, f2 = 0.10)
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
spread_model$set_seeds(setNames(seeds_adm1$k13_valid_q_50, seeds_adm1$id_1))

# Set up our selection speed map data
spread_model$set_map_data(map_data[[which(apply(scenarios, 1, function(x){all(x == "med")}))]] %>% filter(id_1 %in% afr_map$id_1))

# Simulate our data
out <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1, s_name =  "s_a_5")

# set up parallel cluster here
n.cores <- 8
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
future::plan(future::cluster, workers = my.cluster)


# Now simulate spread for all our scenarios (15m to run)
out_all <- furrr::future_map(map_data, function(x){
  spread_model$set_map_data(x %>% filter(id_1 %in% afr_map$id_1))
  central <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1, s_name = "s_a_5") %>%
    group_by(id_1) %>%
    summarise(t = t[freq > 0.1][1]) %>%
    left_join(x %>% filter(id_1 %in% afr_map$id_1) %>% select(id_1, s = s_a_5))
  spread_model$set_map_data(x %>% filter(id_1 %in% afr_map$id_1) %>% mutate(s_a_5 = s_a_5_max))
  worst <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1) %>%
    group_by(id_1) %>%
    summarise(tmin = t[freq > 0.1][1]) %>%
    left_join(x %>% filter(id_1 %in% afr_map$id_1) %>% select(id_1, s_max = s_a_5_max))
  spread_model$set_map_data(x %>% filter(id_1 %in% afr_map$id_1) %>% mutate(s_a_5 = s_a_5_min))
  best <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1) %>%
    group_by(id_1) %>%
    summarise(tmax = t[freq > 0.1][1]) %>%
    left_join(x %>% filter(id_1 %in% afr_map$id_1) %>% select(id_1, s_min = s_a_5_min))

  out_test <- left_join(central, worst) %>% left_join(best) %>%
    select(id_1, s, s_min, s_max, t, tmin, tmax)

  return(out_test)

}, .progress = TRUE, .options = furrr_options(packages = "arms"))
parallel::stopCluster(my.cluster)
private$simulate_selection(t = t, res_pos = res_pos_list[[t]], s_name = s_name)
# --------------------------------------------------------------------------#
# 3. Plot our spread -----------------------------------------------------
# --------------------------------------------------------------------------#

isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])

# Make a map at each t interval
pl_list <- vector("list", length(unique(out$t)))
for(i in c(1, as.integer(length(unique(out$t))/4)+1, as.integer(length(unique(out$t))/2)+1)) {
  x <- unique(out$t)[i]
  pl_list[[i]] <-
    afr_map %>% filter(id_1 %in% out$id_1) %>%
    left_join(out %>% filter(t == x) %>% select(-t) %>% rename(t = freq)) %>%
    ggplot() +
    geom_sf(aes(fill = t), color = "grey", show.legend = TRUE, lwd = 0.1) +
    geom_sf(fill = NA, color = "black", show.legend = FALSE,
            data = map_0 %>% filter(iso %in% isos), lwd = 0.2) +
    coord_sf() +
    theme_void(base_size = 16) +
    ggtitle(as.integer(2020+x)) +
    theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white")) +
    scale_fill_viridis_c(name = "K13 Mutation Frequency (%)",
                         labels = scales::percent, limits = c(0,1)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.height = unit(1, "inch"),
          legend.key.width = unit(0.5, "inch"),
          legend.text = element_text(size = 14),
          legend.title.align = 0)
}

spread_gg <- cowplot::plot_grid(
  pl_list[[1]] + theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("April 2020"),
  pl_list[[as.integer(length(pl_list)/4)+1]] +
    theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("April 2030"),
  pl_list[[as.integer(length(pl_list)/2)+1]] +
    theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("April 2040"),
  cowplot::get_legend(pl_list[[1]] + theme(text = element_text("Helvetica"))),
  ncol = 4, rel_widths = c(1,1,1,0.5), label_size = 18)
save_figs("spread_africa", spread_gg, width = 25, height = 8, pdf_plot = FALSE, font_family = "Helvetica")

# --------------------------------------------------------------------------#
# 4. Plot our output video -----------------------------------------------------
# --------------------------------------------------------------------------#

# set up parallel cluster here
n.cores <- 8
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
future::plan(future::cluster, workers = my.cluster)

# Plot all to the time series directory
dir.create("analysis/plots/time_series")
out_plots <- furrr::future_map(seq_along(unique(out$t)), function(i) {

  x <- unique(out$t)[i]
  gg_x <-
    afr_map %>% filter(id_1 %in% out$id_1) %>%
    left_join(out %>% filter(t == x) %>% select(-t) %>% rename(t = freq)) %>%
    ggplot() +
    geom_sf(aes(fill = t), color = "grey", show.legend = TRUE, lwd = 0.1) +
    geom_sf(fill = NA, color = "black", show.legend = FALSE,
            data = map_0 %>% filter(iso %in% isos), lwd = 0.2) +
    coord_sf() +
    theme_void(base_size = 16) +
    ggtitle(as.integer(2020+x)) +
    theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white")) +
    scale_fill_viridis_c(name = "K13 Mutation Frequency (%)",
                         labels = scales::percent, limits = c(0,1)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.height = unit(1, "inch"),
          legend.key.width = unit(0.5, "inch"),
          legend.text = element_text(size = 14),
          legend.title.align = 0)

  arms:::save_figs(paste0("test_",sprintf("%03d",i)), fig = gg_x,
                    width = 12, height = 12, pdf_plot = FALSE, plot_dir = "analysis/plots/time_series/",
                    font_family = "Helvetica")

  return(1L)

}, .progress = TRUE, .options = furrr_options(packages = "arms"))
parallel::stopCluster(my.cluster)

# Make a movie of these
mapmate::ffmpeg(dir = "analysis/plots/time_series/",
                output_dir = "analysis/plots/time_series/",
                pattern = "test_%03d.png", output = "test_video.mp4",
                delay = 1/12, overwrite = TRUE
)

file.remove(grep("test_\\d", list.files("analysis/plots/time_series/", full.names = TRUE), value = TRUE))

# ---------------------------------------------------- #
# 5. Creating Data Outputs  ----
# ---------------------------------------------------- #

compl <- list()
for(i in seq_along(out_all)) {
  compl[[i]] <- left_join(afr_map, out_all[[i]], by = "id_1") %>%
    sf::st_drop_geometry() %>%
    select(iso, name_0, name_1, t, tmin, tmax, s, s_min, s_max) %>%
    mutate(micro210_scen = scenarios$micro210[i],
           ft_scen = scenarios$ft[i],
           al_scen = scenarios$al[i],
           asaq_scen = scenarios$asaq[i],
           dhappq_scen = scenarios$dhappq[i],
           res_scen = scenarios$res[i]
    )
}

compl_df <- do.call(rbind, compl)
saveRDS(compl_df, "analysis/data-derived/complete_prospective_times.rds")

isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])

# Times censored at 40 years
central_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(if_all(ends_with("scen"), ~ . == "med")) %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  select(continent:s_max) %>%
  mutate(t_med = round(replace(t, is.na(t) & iso %in% isos, Inf),2),
         t_uci = round(replace(tmin, is.na(t) & iso %in% isos, Inf),2),
         t_lci = round(replace(tmax, is.na(t) & iso %in% isos, Inf),2),
         s_med = round(s,4),
         s_lci = round(s_min,4),
         s_uci = round(s_max,4),
         scenario = "Central") %>%
  mutate(t_med = replace(t_med, t_med>40, "40+"),
         t_uci = replace(t_uci, t_uci>40, "40+"),
         t_lci = replace(t_lci, t_lci>40, "40+")) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

worst_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(micro210_scen == "low" & ft_scen == "high" & al_scen == "high", res_scen == "high") %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  select(continent:s_max) %>%
  mutate(t_med = round(replace(t, is.na(t) & iso %in% isos, Inf),2),
         t_uci = round(replace(tmin, is.na(t) & iso %in% isos, Inf),2),
         t_lci = round(replace(tmax, is.na(t) & iso %in% isos, Inf),2),
         s_med = round(s,4),
         s_lci = round(s_min,4),
         s_uci = round(s_max,4),
         scenario = "Pessimistic") %>%
  mutate(t_med = replace(t_med, t_med>40, "40+"),
         t_uci = replace(t_uci, t_uci>40, "40+"),
         t_lci = replace(t_lci, t_lci>40, "40+")) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

best_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(micro210_scen == "high" & ft_scen == "low" & al_scen == "low", res_scen == "low") %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  select(continent:s_max) %>%
  mutate(t_med = round(replace(t, is.na(t) & iso %in% isos, Inf),2),
         t_uci = round(replace(tmin, is.na(t) & iso %in% isos, Inf),2),
         t_lci = round(replace(tmax, is.na(t) & iso %in% isos, Inf),2),
         s_med = round(s,4),
         s_lci = round(s_min,4),
         s_uci = round(s_max,4),
         scenario = "Optimistic") %>%
  mutate(t_med = replace(t_med, t_med>40, "40+"),
         t_uci = replace(t_uci, t_uci>40, "40+"),
         t_lci = replace(t_lci, t_lci>40, "40+")) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

write.csv(central_df, "analysis/data-out/central_times_prospective.csv", row.names = FALSE)
write.csv(worst_df, "analysis/data-out/pessimistic_times_prospective.csv", row.names = FALSE)
write.csv(best_df, "analysis/data-out/optimistic_times_prospective.csv", row.names = FALSE)

# and out put the shape files
dir.create("analysis/data-out/shape_files")
afr_map_full <- readRDS(here::here("analysis/data-derived/map_for_CAR.rds"))
saveRDS(afr_map_full %>% select(iso:country_level), "analysis/data-out/shape_files/admin1.rds")
sf::write_sf(afr_map_full %>% select(iso:country_level), "analysis/data-out/shape_files/admin1.shp")
saveRDS(map_0, "analysis/data-out/shape_files/admin0.rds")
sf::write_sf(map_0, "analysis/data-out/shape_files/admin0.shp")

