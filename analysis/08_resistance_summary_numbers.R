resdfs <- readRDS("analysis/data-derived/resistance_dfs.rds")

gg1 <- resdfs$crt_76T %>% filter(continent == "Africa") %>%
  group_by(syear) %>% summarise(n = n()) %>%
  filter(syear > 2000 & syear < 2020) %>% ggplot(aes(syear, y = n)) +
  geom_bar(stat = "identity") + theme_bw() + scale_y_continuous(expand = c(0,0)) +
  ylab("Data Points") + xlab("Year") + ggtitle("pfcrt 76T") +
  scale_x_continuous(limits = c(2000,2020)) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", color = "white"),
                     text = ggplot2::element_text(family = "Helvetica"))

gg2 <- resdfs$mdr1_86Y %>% filter(continent == "Africa") %>%
  group_by(syear) %>% summarise(n = n()) %>%
  filter(syear > 2000 & syear < 2020) %>% ggplot(aes(syear, y = n)) +
  geom_bar(stat = "identity") + theme_bw() + scale_y_continuous(expand = c(0,0)) +
  ylab("Data Points") + xlab("Year") + ggtitle("mdr1 86Y")+
  scale_x_continuous(limits = c(2000,2020)) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", color = "white"),
                     text = ggplot2::element_text(family = "Helvetica"))

gg3 <- resdfs$mdr1_184F %>% filter(continent == "Africa") %>%
  group_by(syear) %>% summarise(n = n()) %>%
  filter(syear > 2000 & syear < 2020) %>% ggplot(aes(syear, y = n)) +
  geom_bar(stat = "identity") + theme_bw() + scale_y_continuous(expand = c(0,0)) +
  ylab("Data Points") + xlab("Year") + ggtitle("mdr1 184F")+
  scale_x_continuous(limits = c(2000,2020)) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", color = "white"),
                     text = ggplot2::element_text(family = "Helvetica"))

gg4 <- resdfs$mdr1_CNV %>% filter(continent == "Africa") %>%
  group_by(syear) %>% summarise(n = n()) %>%
  filter(syear > 2000 & syear < 2020) %>% ggplot(aes(syear, y = n)) +
  geom_bar(stat = "identity") + theme_bw() + scale_y_continuous(expand = c(0,0)) +
  ylab("Data Points") + xlab("Year") + ggtitle("mdr1 CNV")+
  scale_x_continuous(limits = c(2000,2020)) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", color = "white"),
                     text = ggplot2::element_text(family = "Helvetica"))

gg5 <- resdfs$k13_valid %>% filter(continent == "Africa") %>%
  group_by(syear) %>% summarise(n = n()) %>%
  filter(syear > 2000 & syear < 2020) %>% ggplot(aes(syear, y = n)) +
  geom_bar(stat = "identity") + theme_bw() + scale_y_continuous(expand = c(0,0)) +
  ylab("Data Points") + xlab("Year") + ggtitle("pfk13 validated markers")+
  scale_x_continuous(limits = c(2000,2020)) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", color = "white"),
                     text = ggplot2::element_text(family = "Helvetica"))

gg6 <- resdfs$pfpm2_CNV %>% filter(continent == "Africa") %>%
  group_by(syear) %>% summarise(n = n()) %>%
  filter(syear > 2000 & syear < 2020) %>% ggplot(aes(syear, y = n)) +
  geom_bar(stat = "identity") + theme_bw() + scale_y_continuous(expand = c(0,0)) +
  ylab("Data Points") + xlab("Year") + ggtitle("pfpm2-3 CNV")+
  scale_x_continuous(limits = c(2000,2020)) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", color = "white"),
                     text = ggplot2::element_text(family = "Helvetica"))

samples_gg <- cowplot::plot_grid(gg1,gg2,gg3,gg4,gg5,gg6, ncol = 6)
save_figs("res_sample_n", samples_gg, width = 16, height = 4, pdf_plot = FALSE)
