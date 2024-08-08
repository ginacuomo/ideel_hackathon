library(tidyverse)

# ---------------------------------------------------- o
# 1. Get our resistance modelled data and maps together ----
# ---------------------------------------------------- o

# read in our res from the CAR modelling
res_car <- readRDS(here::here("analysis/data-derived/CAR_modelled_res.rds"))

# read in our scenario map
scenario_maps <- readRDS(here::here("analysis/data-derived/scenario_maps_full.rds"))
map_1 <- sf::st_make_valid(scenario_maps$map)
isos <- na.omit(unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"]))
map_0 <- readRDS("analysis/data-derived/admin0_sf.rds") %>%
  filter(iso %in% isos) %>% sf::st_as_sf()

# ---------------------------------------------------- o
# 2. Function to plot median and upper and lower for a marker ----
# ---------------------------------------------------- o


create_comb_plot <- function(map_1, map_0, res_car, scenario = 14, marker = "crt_76T",
                             scale_lims = c(0,1),
                             lablet = rep("",6), legend_title = "\npfcrt 76T \nFrequency (%)\n"){

  medgg <- map_1  %>%
    left_join(scenario_maps$map_data[[scenario]]) %>%
    filter(!is.na(micro210)) %>%
    left_join(res_car %>% filter(mut == marker), by = "id_1") %>%
    ggplot() +
    geom_sf(aes(fill = q_50))  +
    scale_fill_gradientn(colours = c("#FEECE3", 'red'),labels = scales::percent,
                         values = c(0, 1), name = legend_title, limits = scale_lims) +
    geom_sf(fill = NA, color = "black", show.legend = FALSE,
            data = map_0, lwd = 0.2) +
    coord_sf() +
    theme_bw() +
    theme_void() +
    theme(plot.caption = element_text(face = "italic"),
          plot.background = element_rect(fill = "white", color = "white"),
          text = element_text("Helvetica"),
          plot.title = element_text(hjust = 0.5))

  lowgg <- map_1  %>%
    left_join(scenario_maps$map_data[[scenario]]) %>%
    filter(!is.na(micro210)) %>%
    left_join(res_car %>% filter(mut == marker), by = "id_1") %>%
    ggplot() +
    geom_sf(aes(fill = q_025))  +
    scale_fill_gradientn(colours = c("#FEECE3", 'red'),labels = scales::percent,
                         values = c(0, 1), name = legend_title, limits = scale_lims) +
    geom_sf(fill = NA, color = "black", show.legend = FALSE,
            data = map_0, lwd = 0.2) +
    coord_sf() +
    theme_bw() +
    theme_void() +
    theme(plot.caption = element_text(face = "italic"),
          plot.background = element_rect(fill = "white", color = "white"),
          text = element_text("Helvetica"),
          plot.title = element_text(hjust = 0.5))

  highgg <- map_1  %>%
    left_join(scenario_maps$map_data[[scenario]]) %>%
    filter(!is.na(micro210)) %>%
    left_join(res_car %>% filter(mut == marker), by = "id_1") %>%
    ggplot() +
    geom_sf(aes(fill = q_975))  +
    scale_fill_gradientn(colours = c("#FEECE3", 'red'),labels = scales::percent,
                         values = c(0, 1), name = legend_title, limits = scale_lims) +
    geom_sf(fill = NA, color = "black", show.legend = FALSE,
            data = map_0, lwd = 0.2) +
    coord_sf() +
    theme_bw() +
    theme_void() +
    theme(plot.caption = element_text(face = "italic"),
          plot.background = element_rect(fill = "white", color = "white"),
          text = element_text("Helvetica"),
          plot.title = element_text(hjust = 0.5))




  comb_gg <- cowplot::plot_grid(
    lowgg + ggtitle(label = "2.5% Percentile") + theme(legend.position = "none"), NA,
    medgg + ggtitle(label = "50% Percentile") + theme(legend.position = "none"),NA,
    highgg + ggtitle(label = "97.5% Percentile") + theme(legend.position = "none"),
    cowplot::get_legend(highgg),
    rel_widths = c(1, 0.01,1,0.01, 1,0.2),
    labels = c(lablet[1], "",lablet[2],"", lablet[3],""),
    ncol = 6, scale = 0.95) +
    theme(plot.background = element_rect("white","white"),
          text = element_text("Helvetica"))

  return(comb_gg)
}

# Create our plots
crt_76T_gg <- create_comb_plot(map_1, map_0, res_car, marker = "crt_76T", legend_title = "\npfcrt 76T \nFrequency (%) \n")
mdr1_86Y_gg  <- create_comb_plot(map_1, map_0, res_car, marker = "mdr1_86Y", legend_title = "\npfmdr1 86Y \nFrequency (%) \n")
mdr1_184F_gg  <- create_comb_plot(map_1, map_0, res_car, marker = "mdr1_184F", legend_title = "\npfmdr1 184F \nFrequency (%) \n")
mdr1_CNV_gg  <- create_comb_plot(map_1, map_0, res_car, marker = "mdr1_CNV", legend_title = "\npfmdr1 CNV \nFrequency (%) \n", scale_lims = c(0,0.5))
k13_valid_gg  <- create_comb_plot(map_1, map_0, res_car, marker = "k13_valid", legend_title = "\npfk13 Validated \nMutation \nFrequency (%) \n", scale_lims = c(0,0.5))
pfpm23_CNV_gg  <- create_comb_plot(map_1, map_0, res_car, marker = "pfpm23_CNV", legend_title = "\npfpm23 CNV \nFrequency (%) \n")

# Save our plots
save_figs("crt_76T_modelled_maps", crt_76T_gg, width = 16, height = 5, pdf_plot = FALSE, res=600)
save_figs("mdr1_86Y_modelled_maps", mdr1_86Y_gg, width = 16, height = 6, pdf_plot = FALSE, res=600)
save_figs("mdr1_184F_modelled_maps", mdr1_184F_gg, width = 16, height = 6, pdf_plot = FALSE, res=600)
save_figs("mdr1_CNV_modelled_maps", mdr1_CNV_gg, width = 16, height = 6, pdf_plot = FALSE, res=600)
save_figs("k13_valid_modelled_maps", k13_valid_gg, width = 16, height = 6, pdf_plot = FALSE, res=600)
save_figs("pfpm23_CNV_modelled_maps", pfpm23_CNV_gg, width = 16, height = 6, pdf_plot = FALSE, res=600)

