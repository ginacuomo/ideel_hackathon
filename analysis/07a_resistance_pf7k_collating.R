library(tidyverse)

# ---------------------------------------------------- o
# 1. Malariagen wrangle
# ---------------------------------------------------- o

mdf <- data.table::fread("https://www.malariagen.net/wp-content/uploads/2023/11/Pf7_drug_resistance_marker_genotypes.txt") %>%
  rename(sample = V1)
# writing in case they remove it later
write.csv(mdf, "analysis/data-raw/pf7k_raw.csv")

# get the vcf grabbed GTs
mdrex <- read.csv(here::here("analysis/data-raw/mdr.csv")) %>%
  rename(sample = Sample)

# wrangle these to be similar to mdf styling
# For note, there is one sample with the second allele for mdr1_86
# but I can't find what this could be. Sample SPT35516. Setting to NA for now.
aa_annot <- function(ref, alt, ref_aa = "N", alt_aa = "Y") {
  hom_ref <- which(ref == 0 & alt == 0)
  het <- which(ref == 0 & alt == 1)
  hom_alt <- which(ref == 1 & alt == 1)
  aa <- rep(NA, length(ref))
  aa[hom_ref] <- ref_aa
  aa[het] <- paste(ref_aa, alt_aa, sep = ",")
  aa[hom_alt] <- alt_aa
  return(aa)
}
mdrex_df <- mdrex %>%
  mutate(`mdr1_86[Y]` = aa_annot(mdr1_N86Y, mdr1_N86Y_oth, "N", "Y")) %>%
  mutate(`mdr1_184[F]` = aa_annot(mdr1_Y184F, mdr1_Y184F_oth, "Y", "F")) %>%
  select(sample, `mdr1_86[Y]`, `mdr1_184[F]`)

# bind together
mdf <- left_join(mdf, mdrex_df, by = "sample")

# grab the meta information for getting lat long and sample year
meta <- read.csv("analysis/data-raw/Pf7_samples.txt", sep = "\t") # weblink no longer works
mdf <- left_join(mdf, meta %>% rename(sample = Sample), by = "sample")

# sort this into a consistent format that is the same as the WHO/WWWARN information

# 1. Rename annoying markers
mdf <- mdf %>%
  rename(crt_K76T = `crt_76[K]`,
         mdr1_N86Y = `mdr1_86[Y]`,
         mdr1_Y184F = `mdr1_184[F]`,
         k13_markers = `kelch13_349-726_ns_changes`
  )

# select all that are needed
mdf <- mdf %>% select(
  sample,
  crt_K76T, mdr1_N86Y, mdr1_Y184F, mdr1_dup_call, k13_markers, pm2_dup_call,
  study = Study,
  admin_0 = Country,
  admin_1 = Admin.level.1,
  lat = Admin.level.1.latitude,
  long = Admin.level.1.longitude,
  year = Year
)

# now start collapsing into helpful formats

# A. ART
# as of 2 Jan 2024 ART-R markers (valid and candidate)
validated <- c("F446I", "N458Y", "C469Y", "M476I", "Y493H", "R539T", "I543T",
               "P553L", "R561H", "P574L", "C580Y", "R622I", "A675V",
               "P441L", "G449A", "C469F", "A481V", "R515K", "P527H",
               "N537D", "N537I", "G538V", "V568G")
validated <- paste0(validated, collapse = "|")

# First clean the misings
mdf <- mdf %>%
# 1. All blanks in k13 are WT
mutate(k13_markers = replace(k13_markers, k13_markers == "", "WT")) %>%
  # 2. All "-" in k13 are missing (N = 843)
  mutate(k13_markers = replace(k13_markers, k13_markers == "-", NA)) %>%
  # 3. All "!" in k13 are frame-shift in the haplotype, consider missing (N = 4)
  mutate(k13_markers = replace(k13_markers, k13_markers == "!", NA)) %>%
  # 4. All "!*" in k13 are frameshift and unphased het followed by het. consider missing (N = 1)
  mutate(k13_markers = replace(k13_markers, k13_markers == "!*", NA)) %>%
  # 5. All "*" in k13 are unphased het followed by het. consider missing (N = 8)
  mutate(k13_markers = replace(k13_markers, k13_markers == "*", NA))

# FOR UNDERSTANDING
# Upper case = homozygous mutations
# lower case = heterozygous
# lower case without second mutation noted is heterozygous at SNP level but same AA
# , dictates separate and distinct haplotypes (different clones)
# / dictates additional NS changes in the same clone
# * in the context of mutations, these indicate that the sample could not be phased

# From this, anything that is homozygous alt and the alt is a validated marker is resistant (=1)
# anything that is heterozygous alt and both NS are validated is full resistant (=1)
# anything that is heterozygous alt but not both NS are validated, then weighted resistant is full resistant (=0.5)
# This function works with the data as is here, which has only ether two haplotypes reported (i.e. 1 comma (","))
k13_conv <- function(marker) {

  spls <- strsplit(marker, ",", fixed = TRUE)

  lapply(spls, function(x){
    #message(x)
    if (length(x) == 2) {
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

}

# 6. convert our k13 markers into numeric for k13 valid markers
mdf <- mdf %>%
  mutate(k13_valid = k13_conv(k13_markers))

# B. crt_K76T mdr1 SNPs

# Function to convert biallelix res markers
bires_conv <- function(marker, res_AA, sens_AA) {

  spls <- strsplit(marker, ",", fixed = TRUE)

  lapply(spls, function(x){
    #message(x)
    if (length(x) == 2) {
      return(sum(grepl(res_AA, x, ignore.case = TRUE))/2)
    } else {
      if (is.na(x)){
        return(NA)
      } else if (x == sens_AA) {
        return(0)
      } else if (x == res_AA) {
        return(1)
      } else {
        return("ERROR")
      }
    }
  }) %>% unlist

}

# First clean the misings
mdf <- mdf %>%
  # 1. All "-" in crt are missing (N = 26)
  mutate(crt_K76T = replace(crt_K76T, crt_K76T == "-", NA)) %>%
  # 2. convert our crt markers into numeric for crt76T
  mutate(crt_76T = bires_conv(crt_K76T, "T", "K")) %>%
  # 3. convert our mdr markers into numeric for mdr1_86Y
  mutate(mdr1_86Y = bires_conv(mdr1_N86Y, "Y", "N")) %>%
  # 3. convert our mdr markers into numeric for mdr1_184F
  mutate(mdr1_184F = bires_conv(mdr1_Y184F, "F", "Y"))

# C. Duplications - convert the dup calls to resi status and NAs
mdf <- mdf %>%
  mutate(mdr1_CNV = mdr1_dup_call) %>%
  mutate(pfpm23_CNV = pm2_dup_call) %>%
  mutate(mdr1_CNV = replace(mdr1_CNV, mdr1_CNV == -1, NA)) %>%
  mutate(pfpm23_CNV = replace(pfpm23_CNV, pfpm23_CNV == -1, NA))

# FINALLY:
# and work out totals etc
pf7k_res_df <- mdf %>%
  relocate(k13_valid, .before = last_col()) %>%
  pivot_longer(crt_76T:pfpm23_CNV) %>%
  mutate(gene = gsub("(.*)_(.*)", "\\1", name)) %>%
  rename(mut = name) %>%
  mutate(url = "https://pubmed.ncbi.nlm.nih.gov/36864926/",
         pmid = "36864926",
         database = "Pf7k",
         site = paste("Pf7k Study:", study),
         source = NA) %>%
  group_by(
    site, source, admin_0, admin_1, lat, long, year,
    gene, mut, pmid, url, database
  ) %>%
  summarise(n = sum(!is.na(value)),
            x = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c")) %>%
  mutate(study_start_year = NA, study_end_year = NA) %>%
  mutate(prev = x/n) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

saveRDS(wwarn_res_df, here::here("analysis/data-derived/pf7k_res_df.rds"))

