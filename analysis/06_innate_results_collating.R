# ---------------------------------------------------- #
# 0. Get the scenario mapped data  ----
# ---------------------------------------------------- #

library(tidyverse)
devtools::load_all()

# scenario data
scenario_maps <- readRDS("analysis/data-derived/scenario_maps_full.rds")

# ---------------------------------------------------- #
# 1. Creating Data Outputs  ----
# ---------------------------------------------------- #

compl <- list()
for(i in seq_along(scenario_maps$map_data)) {
  compl[[i]] <- left_join(scenario_maps$map, scenario_maps$map_data[[i]]) %>% sf::st_drop_geometry() %>%
    select(iso, name_0, name_1, id_1, t_s_a_5, t_s_a_5_min, t_s_a_5_max,
           s_a_5, s_a_5_min, s_a_5_max, micro210, ft) %>%
    mutate(micro210_scen = scenario_maps$scenarios$micro210[i],
           ft_scen = scenario_maps$scenarios$ft[i],
           al_scen = scenario_maps$scenarios$al[i],
           asaq_scen = scenario_maps$scenarios$asaq[i],
           dhappq_scen = scenario_maps$scenarios$dhappq[i]
    )


}

compl_df <- do.call(rbind, compl)
compl_df <- compl_df %>%
  filter(!is.na(micro210)) %>%
  setNames(gsub("_s_a_5|", "", names(.))) %>%
  setNames(gsub("_a_5|", "", names(.))) %>%
  mutate(resistance = "Artemisinin") %>%
  as_tibble()

saveRDS(compl_df, "analysis/data-derived/complete_times.rds")

# Times censored at 40 years
tlim <- 40
tlim_text <- "40+"
central_df <- compl_df %>% group_by(iso, name_0, name_1, id_1) %>%
  filter(if_all(ends_with("scen"), ~ . == "med")) %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  na.omit %>%
  select(continent:s_max) %>%
  mutate(t_med = round(replace(t, t<0, Inf),2),
         t_lci = round(replace(t_min, t_min<0, Inf),2),
         t_uci = round(replace(t_max, t_max<0, Inf),2),
         s_med = round(s,4),
         s_lci = round(s_min,4),
         s_uci = round(s_max,4),
         scenario = "Central") %>%
  mutate(t_med = replace(t_med, t_med>tlim, tlim_text),
         t_uci = replace(t_uci, t_uci>tlim, tlim_text),
         t_lci = replace(t_lci, t_lci>tlim, tlim_text)) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

worst_df <- compl_df %>% group_by(iso, name_0, name_1, id_1) %>%
  filter(micro210_scen == "low" & ft_scen == "high" & al_scen == "high") %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  na.omit %>%
  select(continent:s_max) %>%
  mutate(t_med = round(replace(t, t<0, Inf),2),
         t_lci = round(replace(t_min, t_min<0, Inf),2),
         t_uci = round(replace(t_max, t_max<0, Inf),2),
         s_med = round(s,4),
         s_lci = round(s_min,4),
         s_uci = round(s_max,4),
         scenario = "Pessimistic") %>%
  mutate(t_med = replace(t_med, t_med>tlim, tlim_text),
         t_uci = replace(t_uci, t_uci>tlim, tlim_text),
         t_lci = replace(t_lci, t_lci>tlim, tlim_text)) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

best_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(micro210_scen == "high" & ft_scen == "low" & al_scen == "low") %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  na.omit %>%
  select(continent:s_max) %>%
  mutate(t_med = round(replace(t, t<0, Inf),2),
         t_lci = round(replace(t_min, t_min<0, Inf),2),
         t_uci = round(replace(t_max, t_max<0, Inf),2),
         s_med = round(s,4),
         s_lci = round(s_min,4),
         s_uci = round(s_max,4),
         scenario = "Optimistic") %>%
  mutate(t_med = replace(t_med, t_med>tlim, tlim_text),
         t_uci = replace(t_uci, t_uci>tlim, tlim_text),
         t_lci = replace(t_lci, t_lci>tlim, tlim_text)) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

write.csv(central_df, here::here("analysis/data-out/central_times.csv"), row.names = FALSE)
write.csv(worst_df, here::here("analysis/data-out/pessimistic_times.csv"), row.names = FALSE)
write.csv(best_df, here::here("analysis/data-out/optimistic_times.csv"), row.names = FALSE)

# ---------------------------------------------------- #
# 2. Creating Longitudinal Outputs  ----
# ---------------------------------------------------- #

# we need to grab all the selection terms for this
comp_all <- list()
for(i in seq_along(scenario_maps$map_data)) {
  comp_all[[i]] <- left_join(scenario_maps$map, scenario_maps$map_data[[i]]) %>% sf::st_drop_geometry() %>%
    select(iso, name_0, name_1,id_1, matches("^s_a_"), micro210, ft, al, asaq, dhappq) %>%
    mutate(micro210_scen = scenario_maps$scenarios$micro210[i],
           ft_scen = scenario_maps$scenarios$ft[i],
           al_scen = scenario_maps$scenarios$al[i],
           asaq_scen = scenario_maps$scenarios$asaq[i],
           dhappq_scen = scenario_maps$scenarios$dhappq[i]
    )
}

comp_all <- do.call(rbind, comp_all)



# Now do the central scenario simulation
md_central <- comp_all %>% group_by(iso, name_0, name_1, id_1) %>%
  filter(if_all(ends_with("scen"), ~ . == "med")) %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  na.omit %>% ungroup

md_central_sim <- innate_simulation(md_central, res_mod,
                                    t_seq = seq(0, 40, 1),
                                    al_lpf = magenta:::drug_create_al()$lpf,
                                    asaq_lpf = magenta:::drug_create_asaq()$lpf,
                                    dhappq_lpf = magenta:::drug_create_dhappq()$lpf)

md_central_out <- md_central_sim %>%
  left_join(md_central %>% select(continent, name_0, name_1, id_1)) %>%
  mutate(scenario = "Central") %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t, matches("artR|tf"), scenario)

# Now do the worst scenario simulation
md_worst <- comp_all %>% group_by(iso, name_0, name_1, id_1) %>%
  filter(micro210_scen == "low" & ft_scen == "high" & al_scen == "high") %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  na.omit %>% ungroup

md_worst_sim <- innate_simulation(md_worst, res_mod,
                                    t_seq = seq(0, 40, 1),
                                    al_lpf = magenta:::drug_create_al()$lpf,
                                    asaq_lpf = magenta:::drug_create_asaq()$lpf,
                                    dhappq_lpf = magenta:::drug_create_dhappq()$lpf)

md_worst_out <- md_worst_sim %>%
  left_join(md_central %>% select(continent, name_0, name_1, id_1)) %>%
  mutate(scenario = "Worst") %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t, matches("artR|tf"), scenario)

# Now do the best scenario simulation
md_best <- comp_all %>% group_by(iso, name_0, name_1, id_1) %>%
  filter(micro210_scen == "high" & ft_scen == "low" & al_scen == "low") %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  na.omit %>% ungroup

md_best_sim <- innate_simulation(md_best, res_mod,
                                  t_seq = seq(0, 40, 1),
                                  al_lpf = magenta:::drug_create_al()$lpf,
                                  asaq_lpf = magenta:::drug_create_asaq()$lpf,
                                  dhappq_lpf = magenta:::drug_create_dhappq()$lpf)

md_best_out <- md_best_sim %>%
  left_join(md_central %>% select(continent, name_0, name_1, id_1)) %>%
  mutate(scenario = "Best") %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, id_1, t, matches("artR|tf"), scenario)

# And save to file
write.csv(md_central_out, here::here("analysis/data-out/central_longitudinal_times.csv"), row.names = FALSE)
write.csv(md_worst_out, here::here("analysis/data-out/pessimistic_longitudinal_times.csv"), row.names = FALSE)
write.csv(md_best_out, here::here("analysis/data-out/optimistic_longitudinal_times.csv"), row.names = FALSE)

