library(tidyverse)

# ---------------------------------------------------- o
# 1. Collate into one nice object ----
# ---------------------------------------------------- o

# read in each database
who_res_df <- readRDS(here::here("analysis/data-derived/who_res_df.rds"))
pf7k_res_df <- readRDS(here::here("analysis/data-derived/pf7k_res_df.rds"))
wwarn_res_df <- readRDS(here::here("analysis/data-derived/wwarn_res_df.rds"))

# bring together
all_res <- rbind(who_res_df, pf7k_res_df, wwarn_res_df)
all_res <- all_res %>% mutate(lat = as.numeric(lat), long = as.numeric(long))
all_res <- all_res %>% mutate(continent = countrycode::countrycode(iso3c, "iso3c", "continent") )
all_res <- all_res %>% mutate(syear = year) %>%
  mutate(syear = replace(syear, syear == 0, NA)) %>%
  mutate(syear = replace(syear, is.na(syear), study_start_year[is.na(syear)])) %>%
  mutate(syear = as.integer(syear))

crt_df <- all_res %>% filter(mut == "crt_76T")
mdr1_86Y_df <- all_res %>% filter(mut == "mdr1_86Y")
mdr1_184F_df <- all_res %>% filter(mut == "mdr1_184F")
mdr1_CNV_df <- all_res %>% filter(mut == "mdr1_CNV")
art_df <- all_res %>% filter(mut == "k13_valid")
ppq_df <- all_res %>% filter(mut == "pfpm23_CNV")

resis <- list(
  "crt_76T" = crt_df,
  "mdr1_86Y" = mdr1_86Y_df,
  "mdr1_184F" = mdr1_184F_df,
  "mdr1_CNV" = mdr1_CNV_df,
  "k13_valid" = art_df,
  "pfpm2_CNV" = ppq_df
)
saveRDS(resis, here::here("analysis/data-derived/resistance_dfs.rds"))

# 1. Get map info--------

world_map_0 <- readRDS("analysis/data-derived/admin0_sf.rds") %>%
  filter(iso %in% afr_isos) %>% sf::st_as_sf()
isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])

# ---------------------------------------------------- o
# 2. crt Mapping ----
# ---------------------------------------------------- o

crt_76T_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = crt_df %>% filter(x == 0) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = crt_df %>% filter(x>0) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\ncrt 76T\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,3000)), limits = c(0,3000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
crt_76T_map
save_figs("crt_76T_map", crt_76T_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 3. mdr1 86Y Mapping ----
# ---------------------------------------------------- o

mdr1_86Y_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = mdr1_86Y_df %>% filter(x == 0) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = mdr1_86Y_df %>% filter(x>0) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\nmdr1 86Y\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,2000)), limits = c(0,2000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
mdr1_86Y_map
save_figs("mdr1_86Y_map", mdr1_86Y_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 4. mdr1 184F Mapping ----
# ---------------------------------------------------- o

mdr1_184F_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = mdr1_184F_df %>% filter(x == 0) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = mdr1_184F_df %>% filter(x>0) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\nmdr1 184F\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,1000)), limits = c(0,1000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
mdr1_184F_map
save_figs("mdr1_184F_map", mdr1_184F_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 5. mdr1 CNV Mapping ----
# ---------------------------------------------------- o

mdr1_CNV_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = mdr1_CNV_df %>% filter(x == 0) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = mdr1_CNV_df %>% filter(x>0) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\nmdr1 CNV\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,1000)), limits = c(0,1000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
mdr1_CNV_map
save_figs("mdr1_CNV_map", mdr1_CNV_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 6. ARTR Mapping ----
# ---------------------------------------------------- o

art_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = art_df %>% filter(x == 0) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = art_df %>% filter(x>0) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                       values = c(0, .Machine$double.eps, 1), name = "\nArtR\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,3000)), limits = c(0,3000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
art_map
save_figs("art_map", art_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 7. PPQ Mapping ----
# ---------------------------------------------------- o

pm2_CNV_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = ppq_df %>% filter(x == 0) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = ppq_df %>% filter(x>0) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\npm23 CNV\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,200)), limits = c(0,200)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
pm2_CNV_map
save_figs("pm23_CNV_map", pm2_CNV_map, width = 6, height = 4.5, pdf_plot = FALSE)


# 2016 onwards only --------------------------------


# ---------------------------------------------------- o
# 2. crt Mapping ----
# ---------------------------------------------------- o

crt_76T_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = crt_df %>% filter(x == 0 & syear >= 2016) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = crt_df %>% filter(x>0 & syear >= 2016) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\ncrt 76T\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,3000)), limits = c(0,3000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
crt_76T_map
save_figs("crt_76T_map_2016", crt_76T_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 3. mdr1 86Y Mapping ----
# ---------------------------------------------------- o

mdr1_86Y_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = mdr1_86Y_df %>% filter(x == 0 & syear >= 2016) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = mdr1_86Y_df %>% filter(x>0 & syear >= 2016) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\nmdr1 86Y\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,2000)), limits = c(0,2000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
mdr1_86Y_map
save_figs("mdr1_86Y_map_2016", mdr1_86Y_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 4. mdr1 184F Mapping ----
# ---------------------------------------------------- o

mdr1_184F_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = mdr1_184F_df %>% filter(x == 0 & syear >= 2016) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = mdr1_184F_df %>% filter(x>0 & syear >= 2016) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\nmdr1 184F\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,1000)), limits = c(0,1000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
mdr1_184F_map
save_figs("mdr1_184F_map_2016", mdr1_184F_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 5. mdr1 CNV Mapping ----
# ---------------------------------------------------- o

mdr1_CNV_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = mdr1_CNV_df %>% filter(x == 0 & syear >= 2016) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = mdr1_CNV_df %>% filter(x>0 & syear >= 2016) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\nmdr1 CNV\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,1000)), limits = c(0,1000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
mdr1_CNV_map
save_figs("mdr1_CNV_map_2016", mdr1_CNV_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 6. ARTR Mapping ----
# ---------------------------------------------------- o

art_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = art_df %>% filter(x == 0 & syear >= 2016) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = art_df %>% filter(x>0 & syear >= 2016) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\nArtR\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,500,3000)), limits = c(0,3000)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
art_map
save_figs("art_map_2016", art_map, width = 6, height = 4.5, pdf_plot = FALSE)

# ---------------------------------------------------- o
# 7. PPQ Mapping ----
# ---------------------------------------------------- o

pm2_CNV_map <- world_map_0 %>%
  filter(iso %in% isos) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE, lwd = 0.2) +
  theme_bw() +
  geom_point(aes(long,lat,size=n), color = "black", alpha = 0.8,data = ppq_df %>% filter(x == 0 & syear >= 2016) %>% filter(continent == "Africa")) +
  geom_point(aes(long,lat,color = prev,size=n), data = ppq_df %>% filter(x>0 & syear >= 2016) %>% filter(continent == "Africa"),alpha=0.8) +
  coord_sf() +
  theme_void() +
  scale_color_gradientn(colours = c('black', "#FEECE3", 'red'),labels = scales::percent,
                        values = c(0, .Machine$double.eps, 1), name = "\npm23 CNV\n", limits = c(0,1)) +
  scale_size_binned(name = "\n\nSample Size \n", range = c(0.01,8), breaks = c(c(0,10,25,50,100,200)), limits = c(0,200)) +
  theme(plot.background = element_rect(fill = "white", color = "white")) +
  guides(colour = guide_colourbar(order = 1))
pm2_CNV_map
save_figs("pm23_CNV_map_2016", pm2_CNV_map, width = 6, height = 4.5, pdf_plot = FALSE)
