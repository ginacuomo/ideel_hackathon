
# Use selection coefficients to simulate af and tf forward
#' @noRd
af_tf_simlation <- function(s,
                             res_mod,
                             t_seq = seq(0, 80, 0.1),
                             al_lpf = magenta:::drug_create_al()$lpf,
                             asaq_lpf = magenta:::drug_create_asaq()$lpf,
                             dhappq_lpf = magenta:::drug_create_dhappq()$lpf){

  af <- lapply(names(s),
               function(x) {
                 res_mod$predict_f2(
                   dat = s,
                   f1 = 0.01,
                   t = t_seq,
                   s_name = x
                 )}
  ) %>% do.call(cbind, .) %>%
    as.data.frame() %>%
    setNames(gsub("s_", "", names(s)))


  tf <- apply(
    af %>% select(starts_with("a_")),
    1,
    calculate_tf,
    al = mdli$al,
    asaq = mdli$asaq,
    dhappq = mdli$dhappq,
    al_lpf = al_lpf,
    asaq_lpf = asaq_lpf,
    dhappq_lpf = dhappq_lpf
  )

  af$tf <- tf
  return(af)

}


# Use map data to conduct innate simulation for ArtR and TF
#' @noRd
innate_simulation <- function(md, res_mod, iso = "iso",
                            t_seq = seq(0, 80, 1),
                             al_lpf = magenta:::drug_create_al()$lpf,
                             asaq_lpf = magenta:::drug_create_asaq()$lpf,
                             dhappq_lpf = magenta:::drug_create_dhappq()$lpf) {

  # split our map data into admins
  mdl <- split(md, md$id_1)

  # create out innate simulations
  md2 <- map(mdl, function(mdli) {

    # get our selection coefficients
    s <- mdli %>% select(matches("^s_a_\\d$"))
    s_min <- mdli %>% select(matches("^s_a_\\d_min"))
    s_max <- mdli %>% select(matches("^s_a_\\d_max"))

    inn <- af_tf_simlation(s, res_mod, t_seq, al_lpf, asaq_lpf, dhappq_lpf)
    inn_min <- af_tf_simlation(s_max, res_mod, t_seq, al_lpf, asaq_lpf, dhappq_lpf)
    inn_max <- af_tf_simlation(s_min, res_mod, t_seq, al_lpf, asaq_lpf, dhappq_lpf)

    res <- data.frame(
      iso = mdli[[iso]],
      id_1 = mdli$id_1,
      t = t_seq,
      artR_lci = inn_min$a_5,
      artR_med = inn$a_5,
      artR_uci = inn_max$a_5,
      tf_lci = inn_min$tf,
      tf_med = inn$tf,
      tf_uci = inn_max$tf
    )
    return(res)

  }, .progress = TRUE) %>% do.call(rbind, .)

  return(md2)

}

