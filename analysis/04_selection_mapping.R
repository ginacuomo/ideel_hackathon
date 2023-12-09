# ---------------------------------------------------- #
# 1. Get the world map and merge together with our covariate ranges
# ---------------------------------------------------- #

# Grab our covariate parameter ranges
covars <-  readRDS("analysis/data-derived/global_covariate_ranges.rds")
isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
covars <- filter(covars, iso3c %in% isos)

# Get our selection model
selection_model <- readRDS("analysis/data-derived/res_mod.rds")

# MAP world map to
world_map <- malariaAtlas::getShp(ISO = na.omit(unique(covars$iso3c)), admin_level = c("admin1")) %>% sf::st_as_sf()
available_admin <- malariaAtlas::listShp(printed = FALSE, admin_level = "admin0")
world_map_0 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin0")) %>% sf::st_as_sf()

# make scenarios for the map
scenarios <- expand.grid(
  # "al" = c("low", "med", "high"),
  # "asaq" = c("low", "med", "high"),
  # "dhappq" = c("low", "med", "high"),
  # "art_res" = c("low", "med", "high"),
  # "pd_res" = c("low", "med", "high"),
  "ft" = c("low", "med", "high"),
  "micro210" = c("low", "med", "high")
) %>%
  setNames(c("ft","micro210"))
  #setNames(c("al","asaq","dhappq","art_res","pd_res","ft","micro210"))

# ---------------------------------------------------- #
# 2. Estimate times to 5% for each scenario worldwide
# ---------------------------------------------------- #

# make scenario map data
map_data <- vector("list", nrow(scenarios))

# loop over and generate selection and times
for(i in seq_along(map_data)) {

  micro210 <- case_when(scenarios$micro210[i] == "med" ~ covars$Micro.2.10_mean,
                        scenarios$micro210[i] == "low" ~ covars$Micro.2.10_low,
                        scenarios$micro210[i] == "high" ~ covars$Micro.2.10_high)

  ft <- case_when(scenarios$ft[i] == "med" ~ covars$ft_mean,
                  scenarios$ft[i] == "low" ~ covars$ft_high,
                  scenarios$ft[i] == "high" ~ covars$ft_low)

  # al <- case_when(scenarios$al[i] == "med" ~ covars$AL,
  #                 scenarios$al[i] == "low" ~ covars$AL_min,
  #                 scenarios$al[i] == "high" ~ covars$AL_max)
  al <- covars$AL

  # asaq <- case_when(scenarios$asaq[i] == "med" ~ covars$ASAQ,
  #                 scenarios$asaq[i] == "low" ~ covars$ASAQ_min,
  #                 scenarios$asaq[i] == "high" ~ covars$ASAQ_max)
  asaq <- covars$ASAQ

  # dhappq <- case_when(scenarios$dhappq[i] == "med" ~ covars$DP,
  #                 scenarios$dhappq[i] == "low" ~ covars$DP_min,
  #                 scenarios$dhappq[i] == "high" ~ covars$DP_max)
  dhappq <- covars$DP

  # art_res <- case_when(scenarios$art_res[i] == "med" ~ 0.01,
  #                      scenarios$art_res[i] == "low" ~ 0.01,
  #                      scenarios$art_res[i] == "high" ~ 0.01)
  art_res <- 0.01

  # ppq_res <- case_when(scenarios$pd_res[i] == "med" ~ 0.001,
  #                      scenarios$pd_res[i] == "low" ~ 0.001,
  #                      scenarios$pd_res[i] == "high" ~ 0.001)
  ppq_res <- 0.001

  # aq_res <- case_when(scenarios$pd_res[i] == "med" ~ 0.001,
  #                      scenarios$pd_res[i] == "low" ~ 0.001,
  #                      scenarios$pd_res[i] == "high" ~ 0.001)
  aq_res <- 0.001

  # lu_res <- case_when(scenarios$pd_res[i] == "med" ~ 0.001,
  #                      scenarios$pd_res[i] == "low" ~ 0.001,
  #                      scenarios$pd_res[i] == "high" ~ 0.001)
  lu_res <- 0.001

  map_data[[i]] <- selection_model$predict(al, asaq, dhappq, art_res, ppq_res, aq_res, lu_res, ft, micro210, "s_a_5")
  map_data[[i]] <- map_data[[i]] %>%
    mutate(id_1 = covars$id_1, .before = 1) %>%
    mutate(iso3c = covars$iso3c, .before = 1)
}

# remove water bodies going forwards
world_map <- filter(world_map, type_1 != "Water Body") %>%
  filter(iso %in% isos)

# group together and save
scenario_maps <- list("scenarios" = scenarios,
                      "map_data" = map_data,
                      "map" = world_map)
saveRDS(scenario_maps, "analysis/data-derived/scenario_maps_full.rds")
scenario_maps <- readRDS("analysis/data-derived/scenario_maps_full.rds")


isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
world <- left_join(world_map, map_data[[5]]) %>%
  filter(iso %in% isos)

gg_map <- world %>%
  #filter(Micro.2.10 > 0.0001) %>%
  ggplot() +
  geom_sf(aes(fill = s_a_5), color = NA, show.legend = TRUE) +
  #geom_sf(aes(fill = scales::rescale(s_a_5)), color = "grey", show.legend = TRUE, lwd = 0.1) +
  scale_fill_viridis_c(name = "Risk", option = "C", direction = -1, values = c(0,0.3, 1)) +
  theme_bw()

gg_map +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = world_map_0 %>% filter(iso %in% isos), lwd = 0.2) +
  coord_sf() +
  theme_void() +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white"))
