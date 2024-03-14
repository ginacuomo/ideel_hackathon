#' @title R6 Class for simulating spread of resistance on a given map
#'
#' @description A simulating model
#'
#' @importFrom R6 R6Class
#' @importFrom spdep poly2nb
#' @importFrom purrr list_rbind
#' @importFrom dplyr filter
#' @importFrom stats setNames weighted.mean
R6_res_spread <- R6::R6Class(
  classname = "res_spread",
  cloneable = FALSE,

  # PUBLIC METHODS
  public = list(

    # INITIALISATION
    #' @description
    #' Create a new res spread model object
    #' @param map Map sf object
    #' @param res_mod res_mod selection model
    #' @param adj_mat Adjacency matrix for `map`. Default = NULL and will calculate internally
    initialize = function(map, res_mod, adj_mat = NULL) {

      private$map <- map
      private$res_mod <- res_mod

      # check adj_mat has the same names
      if(is.null(adj_mat)) {
        adj_mat <- spdep::poly2nb(map)
        names(adj_mat) <- map$id_1
      } else {
        order_names(adj_mat, as.character(map$id_1))
      }
      private$adj_mat <- adj_mat
    },

    #' @description
    #' Set up simulation seeds
    #' @param seeds Named vector or list of regions seeding spread and iniial frequency
    #'   E.g. list("region_1" = 0.2, "region_2" = 0.3, "region_3" = 0.002)
    #'
    set_seeds = function(seeds){
      if(!all(names(seeds) %in% private$map$id_1)) {
        stop("All seeds names must be in spread model map")
      }
      private$seeds <- seeds
      invisible(private$seeds)
    },

    #' @description
    #' Set up map_data
    #' @param map_data Data.frame with region names ("id_1") and selection ("s")
    #'
    set_map_data = function(map_data){
      if(!all(private$map$id_1 %in% map_data$id_1)) {
        stop("map_data must include values for all spread model map regions")
      }
      private$map_data <- map_data
      invisible(private$map_data)
    },

    #' @description
    #' Simulated spread
    #' @param import_freq What frequency does importation result in. Default = 0.01
    #' @param export_freq At what frequency does exportation occur at. Default = 0.25
    #' @param t_end What year does simulation end. Default = 40
    #' @param t_break Gap between time breaks. Default = 1
    #' @param import_gap Number of years for importation to occur over. Default = 1
    #' @param s_name Name of selection coefficient to be simulated
    #'
    simulate_spread = function(import_freq = 0.01,
                               export_freq = 0.25,
                               t_end = 40,
                               t_break = 1,
                               import_gap = 1,
                               s_name = "s_a_5") {

      # set up our results object
      private$set_res_list(t_end = t_end, t_break = t_break)

      # Where are we looking for resistance first
      next_pos <- which(names(private$res_list) %in% names(private$seeds))
      res_pos <- private$find_resistant_regions(t = 1, pos = next_pos)

      # Running vector of regions that have been simulated from
      simulated <- c()

      # Vector of times
      t_s <- seq(0, t_end - t_break, t_break)
      res_pos_list <- vector("list", length(t_s) + as.integer(1/t_break))
      res_pos_list[[1]] <- res_pos

      # Simulate spread process
      for(t in seq_along(t_s)) {

        # 1. Simulate Selection at resistant regions
        private$simulate_selection(t = t, res_pos = res_pos_list[[t]], s_name = s_name)

        # 2. Update which regions have been simulated
        simulated <- c(simulated, res_pos_list[[t]])

        # 3. Which regions are exporting after selections
        export_pos <- private$find_resistant_regions(t = t, res_freq = export_freq, pos = simulated)

        # 4. If a simulated region has exported then we remove it from simulated
        simulated <- setdiff(simulated, export_pos)

        # 5. Simulate Importation only if less than import gap before end
        if (float_leq(t_s[t], (t_end - import_gap))) {
          res_pos_list[[t + as.integer(1/t_break)]] <- private$simulate_importation(
            t = t, export_pos = export_pos,
            import_freq = import_freq,
            t_break = t_break, import_gap = import_gap
          )
        }
      }

      # And return our simulation
      return(purrr::list_rbind(private$res_list))

    },

    #' @description
    #' Simulated spread of multiallelic traits
    #' @param import_freq What frequency does importation result in. Default = 0.01
    #' @param export_freq At what frequency does exportation occur at. Default = 0.25
    #' @param t_end What year does simulation end. Default = 40
    #' @param t_break Gap between time breaks. Default = 1
    #' @param import_gap Number of years for importation to occur over. Default = 1
    #' @param s_names Names of selection coefficient to be simulated
    #'
    simulate_multiallelic_spread = function(import_freq = 0.01,
                                            export_freq = 0.25,
                                            t_end = 40,
                                            t_break = 1,
                                            import_gap = 1,
                                            s_names = paste0("s_a_", 1:6)) {

      # Catch for if just one mutation
      if(length(s_names) == 1) {
        return(simulate_spread(mport_freq,
                               export_freq,
                               t_end,
                               t_break,
                               import_gap,
                               s_names))
      }

      # set up our results object
      private$set_multiallelic_res_list(t_end = t_end, t_break = t_break)

      # Where are we looking for resistance first
      next_pos <- which(names(private$res_list) %in% names(private$seeds))

      res_pos <- map(s_names, function(s){
        private$find_multiallelic_resistant_regions(t = 1, pos = next_pos, s_name = s)
      })
      names(res_pos) <- s_names

      # Running vector of regions that have been simulated from
      simulated <- vector("list", length(s_names))
      names(simulated) <- s_names

      # Vector of times
      t_s <- seq(0, t_end - t_break, t_break)
      res_pos_list <- lapply(seq_along(s_names), function(x) {
        vector("list", length(t_s) + as.integer(1/t_break))
      })
      names(res_pos_list) <- s_names

      # and fill in which are to be simulated in t = 1
      for(s in s_names) {
        res_pos_list[[s]][[1]] <- res_pos[[s]]
      }

      # Simulate spread process
      for(t in seq_along(t_s)) {

        for(s in seq_along(s_names)){

          # what s name is it
          s_name <- s_names[s]

        # 1. Simulate Selection at resistant regions
        private$simulate_multiallelic_selection(t = t, res_pos = res_pos_list[[s_name]][[t]], s_name = s_name)

        # 2. Update which regions have been simulated
        simulated[[s_name]] <- c(simulated[[s_name]], res_pos_list[[s_name]][[t]])

        # 3. Which regions are exporting after selections
        export_pos <- private$find_multiallelic_resistant_regions(t = t, res_freq = export_freq, pos = simulated[[s_name]], s_name = s_name)

        # 4. If a simulated region has exported then we remove it from simulated
        simulated[[s_name]] <- setdiff(simulated[[s_name]], export_pos)

        # 5. Simulate Importation only if less than import gap before end
        if (float_leq(t_s[t], (t_end - import_gap))) {
          res_pos_list[[s_name]][[t + as.integer(1/t_break)]] <- private$simulate_multiallelic_importation(
            t = t, export_pos = export_pos,
            import_freq = import_freq,
            t_break = t_break, import_gap = import_gap,
            s_name = s_name
          )
        }
        }
      }

      # And return our simulation
      return(purrr::list_rbind(private$res_list))

    }


  ),

  private = list(

    # Private Member Variables
    map = NULL,
    adj_mat = NULL,
    res_mod = NULL,
    seeds = NULL,
    map_data = NULL,
    res_list = NULL,

    # Private Member Functions

    # Set Up
    set_res_list = function(t_end, t_break = 1){

      res <- expand.grid("id_1" = private$map$id_1, "t" = c(seq(0, t_end - t_break, t_break), t_end), "freq" = 0)
      res$t_pos <- match(res$t, c(seq(0, t_end - t_break, t_break), t_end))
      res_list <- split(res, res$id_1)

      # make it have the same ordering as our map
      res_list <- res_list[match(private$map$id_1, names(res_list))]

      # set up initial regions
      for(i in seq_along(private$seeds)){
        res_list[[names(private$seeds)[i]]]$freq[1] <- private$seeds[i]
      }

      private$res_list <- res_list

    },

    # Set Multiallelic list
    set_multiallelic_res_list = function(t_end, t_break = 1){

      res <- expand.grid("id_1" = private$map$id_1, "t" = c(seq(0, t_end - t_break, t_break), t_end))
      res$t_pos <- match(res$t, c(seq(0, t_end - t_break, t_break), t_end))
      alleles <- private$get_multiallelic_alleles()
      for(i in seq_along(alleles)) {
        res[[alleles[i]]] <- 0
      }

      res_list <- split(res, res$id_1)

      # make it have the same ordering as our map
      res_list <- res_list[match(private$map$id_1, names(res_list))]

      # set up initial regions
      for(i in seq_along(private$seeds)){
        for(j in seq_along(alleles)) {
          res_list[[names(private$seeds)[i]]][[alleles[j]]][1] <- private$seeds[[i]][[alleles[j]]]
        }
      }

      private$res_list <- res_list

    },

    # Simulation functions
    # function to find the resistant regions
    find_resistant_regions = function(t, res_freq = NULL, pos = NULL) {

      # which positions are we finding
      if(is.null(pos)) {
        pos <- seq_along(private$res_list)
      }

      # what is our comparison criteria
      if(is.null(res_freq)) {
        res_freq_func <- function(x) {x > 0}
      } else {
        res_freq_func <- function(x) {x >= res_freq}
      }

      # Find regions if positions
      if (length(pos) > 0) {
        res <- map_lgl(private$res_list[pos], function(x){
          res_freq_func(x$freq[x$t_pos == t])
        })
        pos[which(res)]
      } else {
        integer(0L)
      }
    },

    # get alleles from seeds
    get_multiallelic_alleles = function() {
      return(setdiff(names(private$seeds[[1]]), "id_1"))
    },

    # function to find the resistant regions across multiple alleles
    find_multiallelic_resistant_regions = function(t, res_freq = NULL, pos = NULL, s_name) {

      # which positions are we finding
      if(is.null(pos)) {
        pos <- seq_along(private$res_list)
      }

      # what is our comparison criteria
      if(is.null(res_freq)) {
        res_freq_func <- function(x) {x > 0}
      } else {
        res_freq_func <- function(x) {x >= res_freq}
      }

      a_name <- gsub("s_", "", s_name)

      # Find regions if positions
      if (length(pos) > 0) {

        res <- map_lgl(pos, function(i){
          x <- private$res_list[[i]]
              res_freq_func(x[[a_name]][x$t_pos == t])
            })
        pos[which(res)]
      } else {
        integer(0L)
      }
    },

    # simulate selection
    simulate_selection = function(t, res_pos, s_name){

      if(length(res_pos) > 0) {

        # Remaining t for this time step
        t_right_pos <- which(private$res_list[[1]]$t_pos >= t)
        t_right <- private$res_list[[1]]$t[t_right_pos]
        t_forward <- t_right - t_right[1]

        # loop over the regions that need updating
        for(i in res_pos) {

          # s for our region
          s_pos <- match(
            names(private$res_list)[i],
            as.character(private$map_data$id_1)
          )
          s <- private$map_data[[s_name]][s_pos]

          # and our update positions
          dat <- list()
          dat[[s_name]] <- s
          private$res_list[[i]]$freq[t_right_pos] <-
            private$res_mod$predict_f2(
              dat = dat,
              f1 = private$res_list[[i]]$freq[t_right_pos[1]],
              t = t_forward,
              s_name = s_name
            )

        }
      }
    },

    # simulate multiallelic_selection
    simulate_multiallelic_selection = function(t, res_pos, s_name){

      if(length(res_pos) > 0) {

        # Remaining t for this time step
        t_right_pos <- which(private$res_list[[1]]$t_pos >= t)
        t_right <- private$res_list[[1]]$t[t_right_pos]
        t_forward <- t_right - t_right[1]

        # get the allele name
        a_name <- gsub("s_", "", s_name)

        # loop over the regions that need updating
        for(i in res_pos) {

          # s for our region
          s_pos <- match(
            names(private$res_list)[i],
            as.character(private$map_data$id_1)
          )
          s <- private$map_data[[s_name]][s_pos]

          # and our update positions
          dat <- list()
          dat[[s_name]] <- s
          private$res_list[[i]][[a_name]][t_right_pos] <-
            private$res_mod$predict_f2(
              dat = dat,
              f1 = private$res_list[[i]][[a_name]][t_right_pos[1]],
              t = t_forward,
              s_name = s_name
            )

        }
      }
    },

    # find those that reach export freq
    simulate_importation = function(t, export_pos, import_freq, t_break = 1, import_gap = 1) {

      # imported regions
      import_pos <- unique(unlist(private$adj_mat[export_pos]))
      import_pos <- import_pos[import_pos != 0]

      # are we importing somewhere
      if(length(import_pos) > 0) {

        # Create a vector of where we are next going to be simulating
        next_res_pos <- c()

        # loop through imported regions
        for(j in import_pos) {

          # time plus one import_gap position
          tp1 <- which(private$res_list[[j]]$t_pos == t + as.integer(import_gap/t_break))

          # what is the frequency after the importation has occurred
          freq <- private$res_list[[j]]$freq[tp1]

          # if less than import then record and import
          if(freq < import_freq) {
            next_res_pos <- c(next_res_pos, j)
            private$res_list[[j]]$freq[tp1] <- import_freq
          }

        }

        # If all the regions being imported into are already higher than import_freq
        if(is.null(next_res_pos)) {
          next_res_pos <- integer(0L)
        }

      } else {

        next_res_pos <- integer(0L)

      }

      # return where we are simulating next time step
      return(next_res_pos)

    },

    # find those that reach export freq for multiallelic simulations
    simulate_multiallelic_importation = function(t, export_pos, import_freq, t_break = 1, import_gap = 1, s_name) {

      # imported regions
      import_pos <- unique(unlist(private$adj_mat[export_pos]))
      import_pos <- import_pos[import_pos != 0]

      # get the a name
      a_name <- gsub("s_", "", s_name)

      # are we importing somewhere
      if(length(import_pos) > 0) {

        # Create a vector of where we are next going to be simulating
        next_res_pos <- c()

        # loop through imported regions
        for(j in import_pos) {

          # time plus one import_gap position
          tp1 <- which(private$res_list[[j]]$t_pos == t + as.integer(import_gap/t_break))

          # what is the frequency after the importation has occurred
          freq <- private$res_list[[j]][[a_name]][tp1]

          # if less than import then record and import
          if(freq < import_freq) {
            next_res_pos <- c(next_res_pos, j)
            private$res_list[[j]][[a_name]][tp1] <- import_freq
          }

        }

        # If all the regions being imported into are already higher than import_freq
        if(is.null(next_res_pos)) {
          next_res_pos <- integer(0L)
        }

      } else {

        next_res_pos <- integer(0L)

      }

      # return where we are simulating next time step
      return(next_res_pos)

    }

  )
)


# Floating less than or equals
#' @noRd
float_leq <- function(x, y, tolerance = 1e-9) {
  return(x < y | abs(x - y) <= tolerance)
}




# Create haplotype frequencies from allele frequencies
#' @noRd
create_haplotype_freq <- function(a) {

  # Calculate Independent Haplotype Frequencies
  calc_haplotype_freq <- function(bo, a) {

    flippos <- which(!as.logical(bo))
    a[flippos] <- 1-a[flippos]
    prod(a)

  }

  # create the binary combinations
  bit_options <- sapply(0:63, intToBits)[1:6, ]

  # calc the freqs
  res <- apply(bit_options, 2, calc_haplotype_freq, a)

  # nomalise just in case
  res/sum(res)
}

# Calculate tf
#' @noRd
calculate_tf <- function(a, al, asaq, dhappq,
                         al_lpf, asaq_lpf, dhappq_lpf) {

  # get hap frequency
  hapf <- create_haplotype_freq(a)

  al_tf <- sum(hapf*al_lpf)
  asaq_tf <- sum(hapf*asaq_lpf)
  dhappq_tf <- sum(hapf*dhappq_lpf)

  # get drug weighted average tf
  1 - (weighted.mean(c(al_tf, asaq_tf, dhappq_tf),
                     c(al, asaq, dhappq)))
}

# Add tf to putput generated by spread model
#' @noRd
add_tf_to_output <- function(out, spread_model,
                             al_lpf = magenta:::drug_create_al()$lpf,
                             asaq_lpf = magenta:::drug_create_asaq()$lpf,
                             dhappq_lpf = magenta:::drug_create_dhappq()$lpf) {

  reout <- map(split(out, out$id_1), function(x){

    meta <- spread_model$.__enclos_env__$private$map_data %>%
      filter(id_1 == x$id_1[1])

    tf <- apply(
      x %>% select(starts_with("a_")),
      1,
      calculate_tf,
      al = meta$al,
      asaq = meta$asaq,
      dhappq = meta$dhappq,
      al_lpf = al_lpf,
      asaq_lpf = asaq_lpf,
      dhappq_lpf = dhappq_lpf
    )

    x$tf <- tf

    return(x)

  }, .progress = TRUE)

  reout <- do.call(rbind, reout)

  return(reout)
}

