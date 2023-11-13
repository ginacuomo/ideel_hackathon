library(tidyverse)

## ----------------------------------------------------o
## 1. Sourcing previous covariates --------------
## ----------------------------------------------------o

# source covariates from hrp2 analysis into data-derived
download.file("https://github.com/OJWatson/hrpup/raw/v0.2.0/analysis/data_derived/global_covariate_ranges.rds",
              here::here("analysis/data-derived/global_covariate_ranges.rds"))

dat <- readRDS(here::here("analysis/data-derived/global_covariate_ranges.rds"))

## ----------------------------------------------------o
## 2. Checking drug ratios observed across Africa --------------
## ----------------------------------------------------o

# check GF for drug ratios
# https://insights.theglobalfund.org/t/Public/views/PriceQualityReportingTransactionSummary/TransactionSummary?iframeSizedToWindow=true&%3Aembed=y&%3AshowAppBanner=false&%3Adisplay_count=no&%3AshowVizHome=no
gf <- readxl::read_xlsx(here::here("analysis/data-raw/gf_pqr_antimalarial.xlsx"))
gf <- gf %>% fill(1:22, .direction = "down")
gf <- gf %>%
  rename(product = Product) %>%
  mutate(product = replace(product, grep("Lumefantrine", product, fixed = TRUE), "AL")) %>%
  mutate(product = replace(product, grep("Amodiaquine+[Sulfadoxine+Pyrimethamine]", product, fixed = TRUE), "AQSP")) %>%
  mutate(product = replace(product, grep("Chloroquine", product, fixed = TRUE), "CQ")) %>%
  mutate(product = replace(product, grep("Artesunate + Mefloquine", product, fixed = TRUE), "ASMQ")) %>%
  mutate(product = replace(product, grep("Artesunate + [Sulfadoxine+Pyrimethamine]", product, fixed = TRUE), "ASSP")) %>%
  mutate(product = replace(product, grep("Artesunate + Amodiaquine", product, fixed = TRUE), "ASAQ")) %>%
  mutate(product = replace(product, grep("Artesunate+Pyronaridine", product, fixed = TRUE), "ASPY")) %>%
  mutate(product = replace(product, grep("Dihydroartemisinin+Piperaquine", product, fixed = TRUE), "DHAPPQ")) %>%
  mutate(product = replace(product, grep("Quinine", product, fixed = TRUE), "QU")) %>%
  mutate(product = replace(product, grep("Sulfadoxine+Pyrimethamine", product, fixed = TRUE), "SP")) %>%
  mutate(product = replace(product, grep("Primaquine", product, fixed = TRUE), "PQ")) %>%
  mutate(product = replace(product, grep("Mefloquine", product, fixed = TRUE), "MQ")) %>%
  mutate(product = replace(product, nchar(product) > 6, "ART"))

# get quants for each item
quants <- read.csv("analysis/data-raw/gf_product_quants.csv")
gf$quant <- quants$Number.of.Drugs[match(gf$`Product Pack`, quants$Product.Description)]
gf$vol <- gf$quant * as.numeric(gsub(",", "", gf$`Pack quantity`, fixed = TRUE))

# convert into relative proportion of drugs in each country
drug_vols <- gf %>%
  rename(country = `Country/Teritorry`) %>%
  mutate(year = lubridate::year(strptime(gf$`Actual Delivery Date`, "%d-%b-%y"))) %>%
  group_by(country, product, year) %>%
  summarise(vol = sum(vol))

drug_vols <- drug_vols %>% group_by(country, year) %>%
  mutate(vol_prop = vol/sum(vol)) %>%
  mutate(continent = countrycode::countrycode(country, "country.name.en", "continent"))

# quick plot
drug_vols %>%
  filter(continent %in% "Africa") %>%
  ungroup %>%
  complete(country, product, year, fill = list(vol = 0, vol_prop = 0)) %>%
  ggplot(aes(year, vol_prop, fill = product)) +
  geom_area(na.rm = TRUE) +
  facet_wrap(~country) +
  theme_bw() +
  xlab("Year") +
  ylab("Proportion of GF Volumes")

# confirmed we need to span across 0 - 1 for all 3 front line ACTs
drug_vols %>% filter(year > 2017) %>%
  group_by(country, product, continent) %>%
  summarise(p = sum(vol)) %>%
  group_by(country, continent) %>%
  mutate(p = p/sum(p)) %>%
  filter(continent == "Africa") %>%
  group_by(product) %>%
  summarise(min = min(p),
            max = max(p))

## ----------------------------------------------------o
## 3. Define ranges for each parameter --------------
## ----------------------------------------------------o

# Define our ranges
al_range <- c(0, 1)
asaq_range <- c(0, 1)
dhappq_range <- c(0, 1)

# Do no fitness for the time being
fitness_range <- 1

# Mutation fixed at the 7 years
mu_range <- c(0.00084)

# Starting resistance
art_range <- c(0, 1)
ppq_range <- c(0, 1)
aq_range <- c(0, 1)
lu_range <- c(0, 1)

# Define ranges for prev as well to span known prev
prev_range <- c(min(dat$Micro.2.10_low, na.rm = TRUE),
                max(dat$Micro.2.10_high, na.rm = TRUE))

# the prev range corresponds to EIRs from 0.2 to 200.
eir_range <- log(c(0.2, 200))

# Define ranges for ft as well
ft_range <- c(min(dat$ft_low, na.rm = TRUE),
              max(dat$ft_high, na.rm = TRUE))
ft_range <- c(0.01, 0.8)

## ----------------------------------------------------o
## 4 Latin Hypercube sampling --------------
## ----------------------------------------------------o

n <- 2908 # 908 for corner cases that have more than 1 drug coverage or all resistance fixed

ranges <- list(
  al = al_range,
  asaq = asaq_range,
  dhappq = dhappq_range,
  art_res = art_range,
  ppq_res = ppq_range,
  aq_res = aq_range,
  lu_res = lu_range,
  eir = eir_range,
  ft = ft_range,
  linked = c(0,1) # linked flag
)

low_res_corners <- list(
  al = al_range,
  asaq = asaq_range,
  dhappq = dhappq_range,
  art_res = c(0, 0.25),
  ppq_res = c(0, 0.25),
  aq_res = c(0, 0.25),
  lu_res = c(0, 0.25),
  eir = eir_range,
  ft = ft_range,
  linked = c(0,1) # linked flag
)

mid_res_corners <- list(
  al = al_range,
  asaq = asaq_range,
  dhappq = dhappq_range,
  art_res = c(0, 0.5),
  ppq_res = c(0, 0.5),
  aq_res = c(0, 0.5),
  lu_res = c(0, 0.5),
  eir = eir_range,
  ft = ft_range,
  linked = c(0,1) # linked flag
)

# Generate corner samples
corner_samples <- expand.grid(lapply(ranges, function(r) r))
corner_samples_low <- expand.grid(lapply(low_res_corners, function(r) r))
corner_samples_mid <- expand.grid(lapply(mid_res_corners, function(r) r))
corner_samples <- rbind(corner_samples, corner_samples_low, corner_samples_mid)

# Generate LHS samples
n_lhs <- n - nrow(corner_samples)/2
lhs_samples <- lhs::maximinLHS(n_lhs, length(ranges))
colnames(lhs_samples) <- names(ranges)

# Scale the LHS samples to the desired ranges
scaled_lhs_samples <- sapply(names(ranges), function(var) {
  lhs_samples[, var] <- scales::rescale(lhs_samples[, var], to = ranges[[var]])
}, simplify = "data.frame")

# for AL and DHAPPQ and ASAQ these must internally add to 1
# this will give an overly clusterd sample so let's use a
# dirichlet with suitable shape to give good coverage with a bit
# more weight at the corners of the 3d simplex
rdirichlet <- function(shape = rep(1, 3)) {

  x <- shape + 1.1

  while(sum(x) != 1) {
  # Draw from distribution
  x <- rgamma(length(shape), shape = shape, rate = max(shape))
  x <- x / sum(x)
  pos <- sample(3, 1)
  x[pos] <- 1.0 - sum(x[-pos])
  }

  x
}

# make our drug ratio draws
drugs <- t(replicate(nrow(scaled_lhs_samples), rdirichlet(rep(0.75,3))))
drugs %>% as.data.frame() %>% ggplot(aes(V1,V2,color=V3)) + geom_point()

# assign these
scaled_lhs_samples[,1:3] <- drugs

# Combine corner samples and LHS samples
combined_samples <- rbind(corner_samples, scaled_lhs_samples)

# Sort the rows by decreasing eir
combined_samples <- combined_samples %>% arrange(desc(eir))

# convert prevalence back from log(EIR) scale
combined_samples$eir <- exp(combined_samples$eir)

# convert linked back to binary
combined_samples$linked <- ifelse(combined_samples$linked < 0.5, 0, 1)

# add on extra fixed vars
combined_samples$mu <- mu_range
combined_samples$fitness <- fitness_range

# remove instances with all mutations fixed
combined_samples <- combined_samples %>% filter(!(art_res == 1 & ppq_res == 1 & aq_res == 1 & lu_res == 1))

# remove instances where drug usage is greater than 1
combined_samples <- combined_samples[which((combined_samples %>% select(1:3) %>% rowSums()) == 1),]

saveRDS(combined_samples, "analysis/data-derived/lhs_sample.rds")

