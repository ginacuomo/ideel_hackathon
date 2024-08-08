library(magenta)
library(tidyverse)

## ----------------------------------------------------o
## 1. Setting up a cluster configuration --------------
## ----------------------------------------------------o

# Setting Up Cluster From New

# Log in to didehpc
credentials = "C:/Users/ow813/.smbcredentials"
options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "ow813")

# not if T is not mapped then map network drive
didehpc::didehpc_config_global(temp=didehpc::path_mapping("tmp",
                                                          "T:",
                                                          "//fi--didef3.dide.ic.ac.uk/tmp",
                                                          "T:"),
                               home=didehpc::path_mapping("OJ",
                                                          "Z:",
                                                          "//wpia-hn/Malaria/nas5_malaria",
                                                          "Z:"),
                               credentials=credentials,
                               cluster = "wpia-hn")
#cluster = "fi--didemrchnb")

# Creating a Context
context_name <- "analysis/context-dide2"

ctx <- context::context_save(
  path = context_name,
  package_sources = conan::conan_sources(
    packages= c(
      "local::scripts/binaries/magenta_1.3.5.zip",
      "local::scripts/binaries/dde_1.0.2.zip",
      "rrq"),
    repos = "https://mrc-ide.github.io/drat/"
  )
)

# set up a specific config for here as we need to specify the large RAM nodes
config <- didehpc::didehpc_config(use_workers = TRUE, parallel = FALSE, cores = 6)

# Configure the Queue
obj <- didehpc::queue_didehpc(ctx, config = config)

## ----------------------------------------------------o
## 1/2. Create an alternate drug efficacy table --------------
## ----------------------------------------------------o

drug_table <- readr::read_csv("https://raw.githubusercontent.com/OJWatson/art_resistance_consensus_modelling/main/analysis/data/uncertainty_drug_tables/EfficacyTable_equal.csv",
                              locale = readr::locale(encoding = "UTF-16"))

dt <- drug_table %>%
  mutate(tbl = 1) %>%
  mutate(
    Genotype = gsub("--", "0", Genotype),
    Genotype = gsub("x", "", Genotype),
    Genotype = gsub("X", "", Genotype),
    Genotype = gsub("NFNF", "NF1", Genotype),
    Genotype = gsub("YYYY", "YY1", Genotype),
    Genotype = gsub("YFYF", "YF1", Genotype),
    Genotype = gsub("NYNY", "NY1", Genotype)
  ) %>%
  select(Genotype, ASAQ, AL, DHAPPQ, tbl) %>%
  unique

# copy our table
new_drug_table <- dt

# what are the gene options
options <- list(c("K", "T"),
                c("N", "Y"),
                c("Y", "F"),
                c("0", "1"),
                c("C", "Y"),
                c("1", "2"))

# set up these as binary bitsets
bit_options <- sapply(0:63,intToBits)[1:6,]
get_gen <- function(x){

  paste0(unlist(lapply(1:6, function(y) {
    options[[y]][as.numeric(x[y])+1]
  })), collapse = "")

}

# where do these occur and then match them
ordered_options <- apply(bit_options, 2, get_gen)
new_drug_table <- new_drug_table[match(ordered_options, new_drug_table$Genotype),]

# assign as drug table
alt_drug_table <- new_drug_table
alt_drug_table$Genotype <- gsub("2$", "1", gsub("1$","0",alt_drug_table$Genotype))


## ----------------------------------------------------o
## 2. Create cluster param submissions --------------
## ----------------------------------------------------o

# read in parms from lhs
parms <- readRDS(here::here("analysis/data-derived/lhs_sample.rds"))
param_list <- vector("list", length = nrow(parms))

# main parameter set up
N <- 100000
tl <- 40
hd <- 20
nl <- 6

# Create our parameter lists to be run
for(i in seq_len(nrow(parms))) {

  drug_list <-
    magenta:::drug_list_create(
      resistance_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
      mft_flag = TRUE,
      artemisinin_loci = 5,
      absolute_fitness_cost_flag = FALSE,
      partner_drug_ratios = c(parms$al[i],
                              parms$asaq[i],
                              parms$dhappq[i]),
      drugs = list(magenta:::drug_create_al(),
                   magenta:::drug_create_asaq(),
                   magenta:::drug_create_dhappq()),
      cost_of_resistance = rep(parms$fitness[i], 6),
      sequential_cycling = -1,
      number_of_drugs = 3,
      sequential_update = -1,
      number_of_resistance_loci = 6
    )

  # reassign the lpfs to be the new drug efficacies
  drug_list$drugs[[1]]$lpf <- alt_drug_table$AL
  drug_list$drugs[[2]]$lpf <- alt_drug_table$ASAQ
  drug_list$drugs[[3]]$lpf <- alt_drug_table$DHAPPQ

  plaf <- matrix(
    c(rep(c(parms$aq_res[i],parms$aq_res[i],
            parms$lu_res[i],parms$lu_res[i],
            parms$art_res[i],
            parms$ppq_res[i]),hd),
      rep(0,nl*tl)),
    ncol=nl,
    byrow=TRUE
  )

  # list to pass to pipeline
  param_list[[i]] <- list(
    EIR = parms$eir[i],
    N = N,
    years = tl+hd,
    itn_cov = 0,
    save_lineages = TRUE,
    irs_cov = 0,
    ft = parms$ft[i],
    num_loci = nl,
    sample_reps = 1,
    mutation_flag=c(rep(FALSE,hd),rep(TRUE,tl)),
    mutation_treated_modifier = 100,
    mutation_rate=rep(parms$mu[i], nl),
    survival_percentage = 0.22,
    genetics_df_without_summarising = TRUE,
    spatial_incidence_matrix = c(rep(1,hd),rep(0,tl)),
    spatial_mosquitoFOI_matrix = c(rep(1,hd),rep(0,tl)),
    human_only_full_save=FALSE,
    spatial_type = "island",
    use_historic_interventions = TRUE,
    seed = as.integer(runif(1, 1, 1000000000)),
    sample_size = c(100,1000),
    sample_states = c(1,2,4),
    ibd_length = 1,
    update_length = 30,
    plaf=plaf,
    update_save = TRUE,
    human_update_save = TRUE,
    summary_saves_only = TRUE,
    housekeeping_list = magenta:::housekeeping_list_create(quiet = TRUE,cluster = TRUE),
    drug_list = drug_list,
    nmf_list = magenta:::nmf_list_create(nmf_flag = TRUE, prob_of_testing_nmf = parms$ft[i]),
    island_imports_plaf_linked_flag = as.logical(parms$linked[i])
  )

}

## ----------------------------------------------------o
## 3. Make cluster submissions --------------
## ----------------------------------------------------


## safe submission
try_fail_catch <- function(expr, attempts = 3){
  r <- NULL
  attempt <- 1
  while( is.null(r) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      r <- eval(expr)
    )
  }

}

# 1. First do the non linked, i.e. independent loci
sub <- which(parms$linked == 0)

for(rep in 1:10){
  message(rep)
  try_fail_catch(
    grp <- obj$lapply(
      X = param_list[sub],
      timeout=0,
      FUN = function(x){
        return(magenta::pipeline(
          EIR = x$EIR,
          seed = as.integer(runif(1, 1, 1000000000)),
          save_lineages = x$save_lineages,
          N = x$N,
          mutation_rate = x$mutation_rate,
          mutation_flag = x$mutation_flag,
          sample_size = x$sample_size,
          sample_reps = x$sample_reps,
          years = x$years,
          survival_percentage = x$survival_percentage,
          itn_cov = x$itn_cov,
          irs_cov = x$irs_cov,
          ft = x$ft,
          genetics_df_without_summarising = x$genetics_df_without_summarising,
          spatial_incidence_matrix = x$spatial_incidence_matrix,
          spatial_mosquitoFOI_matrix = x$spatial_mosquitoFOI_matrix,
          spatial_type = x$spatial_type,
          use_historic_interventions = x$use_historic_interventions,
          human_only_full_save = x$human_only_full_save,
          ibd_length = x$ibd_length,
          num_loci = x$num_loci,
          sample_states = x$sample_states,
          update_length = x$update_length,
          update_save = x$update_save,
          human_update_save = x$human_update_save,
          summary_saves_only = x$summary_saves_only,
          housekeeping_list = x$housekeeping_list,
          drug_list = x$drug_list,
          nmf_list = x$nmf_list,
          plaf = x$plaf,
          island_imports_plaf_linked_flag = x$island_imports_plaf_linked_flag
        ))
      },
      name = paste0("worstcasebankmu_linked_0_nmf_TRUE_rep_", rep), overwrite = TRUE)
  )

}

# clear those from queue that need extra cores
rrq <- obj$rrq_controller()
ql <- rrq$queue_list()
tid <- unlist(lapply(ql, function(x){(rrq$task_data(x)$expr %>% as.list)$task_id}))
eirs <- vapply(tid, function(x){as.list(obj$task_get(x)$expr())[[2]]$EIR}, numeric(1))
rrq$queue_remove(ql[which(eirs>150)])

# now submit workers
didehpc::web_login()
workers <- obj$submit_workers(150)

# and get the grp list here
grps <- lapply(grep("worstcasebankmu",obj$task_bundle_list(), value = TRUE), function(x){obj$task_bundle_get(x)})
st <- lapply(grps, function(x){x$status()})
running <- unlist(lapply(st[1:10], function(x){names(which(x == "RUNNING"))}))
obj$reconcile(running)


# resubmit those that errored or were not submitted in the first place
for(i in seq_along(grps)) {
  st <- grps[[i]]$status()
  if(any(st == "ERROR" )) {
    obj$submit(grps[[i]]$ids[as.integer(which(st == "ERROR"))])
  }
  if(any(st == "PENDING" )) {
    obj$submit(grps[[i]]$ids[as.integer(which(st == "PENDING"))])
  }
}

## ----------------------------------------------------o
## 4. Fetch Objects --------------
## ----------------------------------------------------o

# let's put our sim outputs here
dir.create(here::here("analysis/data-derived/sims_nmf_higher_res"))

# first save the parms used to generate (X)
saveRDS(grps[[1]]$X, here::here("analysis/data-derived/sims_nmf_higher_res/X.rds"))
X <- grps[[1]]$X

# Save our simulations across the parameter space
# 103 114 122 147 149 158 163 169 still to grab
for(i in c(103, 114, 122, 147, 149, 158, 163, 169)) {
  #for(i in seq_along(X)) {

  message(i)

  sti <- map_chr(seq_along(grps), function(x) {
    grps[[x]]$tasks[[i]]$status()
  })

  if(all(sti == "COMPLETE")) {

    # create our results
    r_i <- map(seq_along(grps), function(x) {
      grps[[x]]$db$get_value(grps[[x]]$db$get_hash(grps[[x]]$tasks[[i]]$id, "task_results"), FALSE)
    })

    if(any(unlist(lapply(r_i, object.size))==0)){

    } else {

      # remove the final loggers that are very large and unneeded
      for(x in seq_along(r_i)) {
        to_rm <- which(names(r_i[[x]][[length(r_i[[x]])]]$Loggers) %in%
                         c("InfectionStates", "Ages", "IB", "ICA", "ICM", "ID") )
        r_i[[x]][[length(r_i[[x]])]]$Loggers[to_rm] <- NULL
      }

      # and save to file
      fn_i <- paste0("rnmf_", i, ".rds")
      saveRDS(r_i, here::here("analysis/data-derived/sims_nmf_higher_res", fn_i))
      gc()

    }
  }

}

# These sims have not been pushed to Github due to Github memory/size constraints

# Analysis in next script uses these sims and the output of this is then saved
# in data-derived

