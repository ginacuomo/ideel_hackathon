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
    select(iso, name_0, name_1, t_s_a_5, t_s_a_5_min, t_s_a_5_max,
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
central_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
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
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

worst_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
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
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

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
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, s_lci, s_med, s_uci, scenario)

write.csv(central_df, here::here("analysis/data-out/central_times.csv"), row.names = FALSE)
write.csv(worst_df, here::here("analysis/data-out/pessimistic_times.csv"), row.names = FALSE)
write.csv(best_df, here::here("analysis/data-out/optimistic_times.csv"), row.names = FALSE)

