library(tidyverse)

# ---------------------------------------------------- o
# 1. Get the world map and merge together with our covariate ranges ----
# ---------------------------------------------------- o

# Grab our covariate parameter ranges
covars <-  readRDS("analysis/data-derived/global_covariate_ranges.rds")
isos <- na.omit(unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"]))
covars <- filter(covars, iso3c %in% isos)

# Get our selection model
res_mod <- readRDS("analysis/data-derived/res_mod.rds")

# MAP world map to
world_map <- malariaAtlas::getShp(ISO = na.omit(unique(covars$iso3c)), admin_level = c("admin1")) %>% sf::st_as_sf()
available_admin <- malariaAtlas::listShp(printed = FALSE, admin_level = "admin0")
world_map_0 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin0")) %>% sf::st_as_sf()

# make scenarios for the map
drugs <- c("al", "asaq", "dhappq")
scenarios <- lapply(drugs, function(x){
  expand.grid(
    "ft" = c("low", "med", "high"),
    "micro210" = c("low", "med", "high"),
    x = c("low", "med", "high")
  ) %>%
    setNames(c("ft","micro210", x)) %>%
    mutate(!!sym(setdiff(drugs, x)[1]) := "med") %>%
    mutate(!!sym(setdiff(drugs, x)[2]) := "med")
}) %>%
  do.call(rbind, .) %>%
  unique

# ---------------------------------------------------- o
# 2. Estimate times to 10% for each scenario worldwide ----
# ---------------------------------------------------- o

# make scenario map data
map_data <- vector("list", nrow(scenarios))

# loop over and generate selection and times
for(i in seq_along(map_data)) {

  micro210 <- case_when(scenarios$micro210[i] == "med" ~ covars$Micro.2.10_mean,
                        scenarios$micro210[i] == "low" ~ covars$Micro.2.10_low,
                        scenarios$micro210[i] == "high" ~ covars$Micro.2.10_high)

  ft <- case_when(scenarios$ft[i] == "med" ~ covars$ft_mean,
                  scenarios$ft[i] == "low" ~ covars$ft_low,
                  scenarios$ft[i] == "high" ~ covars$ft_high)

  al <- case_when(scenarios$al[i] == "med" ~ covars$AL,
                  scenarios$al[i] == "low" ~ covars$AL_min,
                  scenarios$al[i] == "high" ~ covars$AL_max)

  asaq <- case_when(scenarios$asaq[i] == "med" ~ covars$ASAQ,
                  scenarios$asaq[i] == "low" ~ covars$ASAQ_min,
                  scenarios$asaq[i] == "high" ~ covars$ASAQ_max)

  dhappq <- case_when(scenarios$dhappq[i] == "med" ~ covars$DP,
                  scenarios$dhappq[i] == "low" ~ covars$DP_min,
                  scenarios$dhappq[i] == "high" ~ covars$DP_max)

  # normalise to sum to 1
  drugmat <- cbind(al, asaq, dhappq)
  drugmat <- apply(drugmat, 1, function(x){x/sum(x)}) %>% t()
  al <- drugmat[,1]
  asaq <- drugmat[,2]
  dhappq <- drugmat[,3]

  # art_res <- case_when(scenarios$art_res[i] == "med" ~ 0.01,
  #                      scenarios$art_res[i] == "low" ~ 0.01,
  #                      scenarios$art_res[i] == "high" ~ 0.01)
  art_res <- 0.01

  # ppq_res <- case_when(scenarios$pd_res[i] == "med" ~ 0.001,
  #                      scenarios$pd_res[i] == "low" ~ 0.001,
  #                      scenarios$pd_res[i] == "high" ~ 0.001)
  ppq_res <- 0.01

  # aq_res <- case_when(scenarios$pd_res[i] == "med" ~ 0.001,
  #                      scenarios$pd_res[i] == "low" ~ 0.001,
  #                      scenarios$pd_res[i] == "high" ~ 0.001)
  aq_res <- 0.01

  # lu_res <- case_when(scenarios$pd_res[i] == "med" ~ 0.001,
  #                      scenarios$pd_res[i] == "low" ~ 0.001,
  #                      scenarios$pd_res[i] == "high" ~ 0.001)
  lu_res <- 0.01

  map_data[[i]] <- res_mod$predict(al, asaq, dhappq, art_res, ppq_res, aq_res, lu_res, ft,
                                   micro210, "s_a_5", f1 = 0.01, f2 = 0.10)
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


# ---------------------------------------------------- o
# 3. Main Figures For Innate Risk Going from s to time ----
# ---------------------------------------------------- o

isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
# range of scenarios c(7,14,21)

create_comb_plot <- function(scen = 14, lablet = c("A","B"), top = 0.51, t_top = 80){

  world <- left_join(world_map, map_data[[scen]]) %>%
    filter(iso %in% isos) %>%
    filter(!is.na(micro210))

# s map
gg_map_s <- world %>% filter(!is.na(s_a_5)) %>%
  ggplot() +
  geom_sf(aes(fill = s_a_5), color = NA, show.legend = TRUE) +
  scale_fill_viridis_c(name = "Selection \nCoefficient\n", option = "A", direction = 1, values = c(0,0.3, 1), limits = c(0,top),
                       breaks = seq(0, 0.5, 0.1), ) +
  theme_bw()

gg_map_s <- gg_map_s +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = world_map_0 %>% filter(iso %in% isos), lwd = 0.2) +
  coord_sf() +
  theme_void() +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white"),
        text = element_text("Helvetica"))


# t map
gg_map_t <- world %>% filter(!is.na(s_a_5)) %>%
  ggplot() +
  geom_sf(aes(fill = t_s_a_5), color = NA, show.legend = TRUE) +
  scale_fill_viridis_c(name = "Years for ArtR\n1% → 10%\n", option = "C", direction = -1,
                       trans = "log", breaks = c(5, 10,20,40,t_top), limits = c(5,t_top)) +
  theme_bw()

gg_map_t <- gg_map_t +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = world_map_0 %>% filter(iso %in% isos), lwd = 0.2) +
  coord_sf() +
  theme_void() +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white"),
        text = element_text("Helvetica"))

# combine plots
arrow <- ggplot() +
  geom_line(arrow = grid::arrow(),lwd = 1,
            data = data.frame(
              x = c(0, 1),
              y = c(0.5, 0.5)
            ),
            aes(x = x, y=  y)
  ) +
  cowplot:::theme_nothing()

comb_gg <- cowplot::plot_grid(gg_map_s,NA,
                              arrow,NA,
                              gg_map_t,
                              labels = c(lablet[1], "","","", lablet[2]),
                              rel_widths = c(1, 0.05,0.2,-0.01, 1),
                              ncol = 5, scale = 0.95) +
  theme(plot.background = element_rect("white","white"),
        text = element_text("Helvetica"))

return(comb_gg)
}

comb_low <- create_comb_plot(7, lablet = c("E","F"))
comb_mid <- create_comb_plot(14, lablet = c("C","D"))
comb_high <- create_comb_plot(21, lablet = c("A","B"))

comb_all <- cowplot::plot_grid(comb_high, comb_mid, comb_low, ncol = 1)
comb_mid_solo <- create_comb_plot(14, lablet = c("A","B"))

# save the comparison conversion figure out
save_figs("selcoef_to_time_all", comb_all, width = 12, height = 12, pdf_plot = FALSE)
save_figs("selcoef_to_time", comb_mid_solo, width = 12, height = 4, pdf_plot = FALSE)


