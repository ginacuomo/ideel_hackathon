library(tidyverse)

central_times <- read_csv("analysis/data-out/central_times.csv")
best_times <- read_csv("analysis/data-out/optimistic_times.csv")
worst_times <- read_csv("analysis/data-out/pessimistic_times.csv")

uganda_comp <- rbind(
central_times %>% filter(iso == "UGA") %>% mutate(scen = "Central"),
worst_times %>% filter(iso == "UGA") %>% mutate(scen = "Worst"),
best_times %>% filter(iso == "UGA") %>% mutate(scen = "Best")
) %>%
  filter(admin_1 %in% c("Kaabong", "Lamwo", "Koboko", "Arua", "Kole",
                        "Agago", "Katakwi", "Kapchorwa", "Tororo", "Jinja",
                        "Hoima", "Mubende", "Kasese", "Kanungu")) %>%
  group_by(scen) %>%
  summarise(med = median(s_med),
            low = quantile(s_med, 0.025),
            high = quantile(s_med, 0.975)
  ) %>%
  ggplot(aes(y = med, ymin = low, ymax = high, x = scen)) +
  geom_rect(ymin = 0.528, ymax = 0.247,
            xmin = -Inf, xmax = Inf, fill = '#00008B', alpha = 0.03) +
  geomtextpath::geom_texthline(yintercept = c(0.383), color = "#00008B", lwd = 1, family = "Helvetica",
                               label = "Median Uganda K13 Selection Coefficient \n(Meier-Scherling et al. 2024)") +
  geomtextpath::geom_texthline(yintercept = c(0.247), linetype = "dashed", color = "#00008B",  family = "Helvetica",
                               label = "2.5% CrI Uganda K13 Selection Coefficient \n(Meier-Scherling et al. 2024)") +
  geomtextpath::geom_texthline(yintercept = c(0.528), linetype = "dashed", color = "#00008B",  family = "Helvetica",
                               label = "97.5% CrI Uganda K13 Selection Coefficient \n(Meier-Scherling et al. 2024)") +
  geom_pointrange() +
  theme_bw() +
  ylab("Modelled K13 Selection Coeffients") +
  xlab("Scenarios")


save_figs("uganda_comp", uganda_comp, width = 8, height = 6, pdf_plot = FALSE, font_family = "Helvetica")
