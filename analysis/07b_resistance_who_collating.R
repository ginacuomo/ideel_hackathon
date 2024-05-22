library(tidyverse)

# ---------------------------------------------------- o
# 1. WHO database compile ----
# ---------------------------------------------------- o

# read in WHO study info
whodf <- readxl::read_xlsx("analysis/data-raw/WHO_res_database_02-01-2024.xlsx", sheet = 2)

# read in res WHO info
whores <- readxl::read_xlsx("analysis/data-raw/WHO_res_database_02-01-2024.xlsx", sheet = 3)

# combine by the res
who <- left_join(whores, whodf, by = "ID")

# sort names as wanted
who <- who %>%
  rename(admin_0 = COUNTRY_NAME) %>%
  rename(admin_1 = ADMIN2) %>%
  mutate(iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c")) %>%
  rename(lat = LATITUDE) %>%
  rename(long = LONGITUDE) %>%
  rename(study_start_year = YEAR_START) %>%
  mutate(year = NA) %>%
  mutate(study_end_year = NA) %>%
  rename(n = SAMPLE_SIZE) %>%
  mutate(prev = as.numeric(PROPORTION)/100) %>%
  mutate(x = round(n*prev)) %>%
  rename(gene = MM_TYPE) %>%
  rename(mut = GENOTYPE) %>%
  rename(url = CITATION_URL) %>%
  # happen to know that this is the Conrad paper early
  mutate(url = replace(url, url == "99999", "http://www.ncbi.nlm.nih.gov/pubmed/37611122")) %>%
  mutate(pmid = gsub(".*pubmed/(\\d+).*|.*gov/(\\d+).*","\\1",url)) %>%
  mutate(pmid = replace(pmid, pmid == "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7029593/", "32070355")) %>%
  mutate(pmid = replace(pmid, pmid == "", NA)) %>%
  mutate(database = "WHO") %>%
  rename(source = DATA_SOURCE) %>%
  rename(site = SITE_NAME) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

# standardise gene names
who <- who %>%
  mutate(
    gene = replace(gene, gene == "Pfkelch13", "k13"),
    gene = replace(gene, gene == "Pfcrt", "crt"),
    gene = replace(gene, gene == "Pfmdr1", "mdr1"),
    gene = replace(gene, gene == "Pfplasmepsin 2-3", "pfpm23")
  )

# split these out into each marker type
whospl <- who %>% split(who$gene)

# 1. Handle each separately for ease - crt

# sanity check that data is correctly formatted from WHO
# i.e. that all samples for a uuid grouping sum to the
# number of samples tested,
whospl$crt %>%
  mutate(mut = replace(mut, mut == "Pfcrt", "crt_76T")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup() %>%
  group_by(uuid) %>%
  summarise(p = sum(prev)) %>%
  pull(p) %>%
  all()

whospl$crt <- whospl$crt %>%
  mutate(mut = replace(mut, mut == "Pfcrt", "crt_76T")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  summarise(n = sum(unique(n)), x = n - sum(x[mut == "WT"]), prev = x/n) %>%
  mutate(mut = "crt_76T") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source) %>%
  ungroup

# 2. Handle each separately for ease - mdr1

# sanity check that data is correctly formatted from WHO
# i.e. that all samples for a uuid grouping sum to the
# number of samples tested,
whospl$mdr1 %>%
  mutate(mut = replace(mut, mut == "MC", "mdr1_CNV")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup() %>%
  group_by(uuid) %>%
  summarise(p = sum(prev)) %>%
  pull(p) %>%
  all()

whospl$mdr1 <- whospl$mdr1 %>%
  mutate(mut = replace(mut, mut == "MC", "mdr1_CNV")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  summarise(n = sum(unique(n)), x = n - sum(x[mut == "WT"]), prev = x/n) %>%
  mutate(mut = "mdr1_CNV") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source) %>%
  ungroup

# 3. Handle each separately for ease - pfpm23

# sanity check that data is correctly formatted from WHO
# i.e. that all samples for a uuid grouping sum to the
# number of samples tested,
whospl$pfpm23 %>%
  mutate(mut = replace(mut, mut == "MC", "pfpm23_CNV")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup() %>%
  group_by(uuid) %>%
  summarise(p = sum(prev)) %>%
  pull(p) %>%
  all()

whospl$pfpm23 <- whospl$pfpm23 %>%
  mutate(mut = replace(mut, mut == "MC", "pfpm23_CNV")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  summarise(n = sum(unique(n)), x = n - sum(x[mut == "WT"]), prev = x/n) %>%
  mutate(mut = "pfpm23_CNV") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source) %>%
  ungroup

# 4. Handle each separately for ease - k13_valid

# A. ART
# as of 2 Jan 2024 ART-R markers (valid and candidate)
validated <- c("F446I", "N458Y", "C469Y", "M476I", "Y493H", "R539T", "I543T",
               "P553L", "R561H", "P574L", "C580Y", "R622I", "A675V",
               "P441L", "G449A", "C469F", "A481V", "R515K", "P527H",
               "N537D", "N537I", "G538V", "V568G")
validated <- paste0(validated, collapse = "|")

# sanity check that data is correctly formatted from WHO
# i.e. that all samples for a uuid grouping sum to the
# number of samples tested,
whospl$k13 %>%
  mutate(mut = replace(mut, grepl(validated, mut, ignore.case = TRUE), "k13_valid")) %>%
  mutate(mut = replace(mut, mut != "k13_valid", "WT")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup() %>%
  group_by(uuid) %>%
  summarise(p = sum(prev)) %>%
  pull(p) %>%
  all()


# need to account for mixed infections
marker <- whospl$k13$mut %>% unique()
spls <- strsplit(marker, "&", fixed = TRUE)
prevw <- lapply(spls, function(x){
  #message(x)
  if (length(x) == 3) {
    return(sum(grepl(validated, x, ignore.case = TRUE))/3)
  } else if (length(x) == 2) {
    return(sum(grepl(validated, x, ignore.case = TRUE))/2)
  } else {
    if (is.na(x)){
      return(NA)
    } else if (x == "WT") {
      return(0)
    } else {
      return(as.integer(grepl(validated, x, ignore.case = TRUE)))
    }
  }
}) %>% unlist
mixed_res_loci <- marker[which(prevw == 0.5)]

whospl$k13 <- whospl$k13 %>%
  # Adjust for mixed res and WT infections to try and derive frequency
  mutate(x = replace(x, mut %in% mixed_res_loci, x[mut %in% mixed_res_loci])) %>%
  mutate(mut = replace(mut, grepl(validated, mut, ignore.case = TRUE), "k13_valid")) %>%
  mutate(mut = replace(mut, mut != "k13_valid", "WT")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  summarise(n = sum(unique(n)), x = n - sum(x[mut == "WT"]), prev = x/n) %>%
  mutate(mut = "k13_valid") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source) %>%
  ungroup


# bring it all back together
who_res_df <- do.call(rbind, whospl)
saveRDS(who_res_df, here::here("analysis/data-derived/who_res_df.rds"))

