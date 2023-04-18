chain_sim_susc_new <- function(offspring = c("pois", "nbinom"),
                           mn_offspring,
                           disp_offspring = NULL,
                           serial,
                           t0 = 0,
                           tf = Inf,
                           pop,
                           initial_immune = 0) {
  # check offspring arguments
  if (!offspring %in% c("pois", "nbinom")) {
    # stop(sprintf(
    #   "%s %s %s",
    #   "Unknown offspring distribution;",
    #   "`chain_sim_susc()` can only generate",
    #   "pois and nbinom offspring"
    # ))
    cli::cli_abort(c(
      "{.arg offspring} must be {.val pois} or {.val nbinom}",
      "x" = "You have supplied {offspring}",
      "i" = "{.fn chain_sim_susc} can only generate poisson and negative
      binomial offspring. Supply either `pois` or `nbinom`"
    ))
  }

  # check that "nbinom" has appropriate accompanying arguments
  if (offspring == "nbinom" && (is.null(disp_offspring) || disp_offspring <= 1)
  ) {
    # stop(sprintf(
    #   "%s %s %s",
    #   "Offspring distribution 'nbinom' requires",
    #   "argument 'disp_offspring' > 1. ",
    #   "Use 'pois' if there is no overdispersion."
    # ))

    cli::cli_abort(c(
      "{.arg disp_offspring} must be greater than 1.",
      "x" = "You have currently supplied it as {disp_offspring}.",
      "i" = "Offspring distribution {.var nbinom} requires argument
      {.arg disp_offspring} > 1. Use {.var pois} if there is no
      overdispersion."
    ))
  }

  # determine offspring distribution

  offspring_fun <- function(n, susc) {
    if (offspring == "pois") {
      #' using a right truncated poisson distribution to prevent sampling
      #' more cases than susceptibles
      pois_mean <- mn_offspring * susc / pop


        truncdist::rtrunc(n,
        spec = "pois",
        lambda = pois_mean,
        b = susc
      )
    } else if (offspring == "nbinom") {
      #' using a right truncated nbinom distribution to prevent sampling
      #' more cases than susceptibles. We are using the parameterization
      #' with mean, 'mu' and dispersion, 'size' (See ?rnbinom).
      new_mn <- mn_offspring * susc / pop ## apply susceptibility

      size <- new_mn / (disp_offspring - 1)
        truncdist::rtrunc(
          n,
          spec = "nbinom",
          b = susc,
          mu = new_mn,
          size = size
        )
      }
    }
  #' initialize simulation data frame and mark as unsimulated
  #' (using offspring_generated)
  tdf <- data.frame(
    id = 1L,
    ancestor = NA_integer_,
    generation = 1L,
    time = t0,
    offspring_generated = FALSE
  )

  # Initial susceptible individuals
  susc <- pop - initial_immune - 1L

  # run simulation until we reach tf or run out of susceptibles
  while (any(tdf$time <= tf & !tdf$offspring_generated) && susc > 0) {
    # get earliest unsimulated time
    earliest_unsim_time <- min(tdf$time[!tdf$offspring_generated])

    # get index of earliest unsimulated case
    idx <- which(tdf$time == earliest_unsim_time & !tdf$offspring_generated)[1]

    # generate offspring
    n_offspring <- offspring_fun(1, susc)
    if (n_offspring %% 1 > 0) {
      # stop("Offspring distribution must return integers")
      cli::cli_abort(c(
        "Sampled offspring must be integers.",
        "x" = "Current specified {.arg offpspring} and associated parameters
        are generating {.cls {class(n_offspring)}} offspring.",
        "i" = "Check the parameter values for {offspring}."
      ))
    }

    # mark offspring as simulated
    tdf$offspring_generated[idx] <- TRUE

    # sample serial intervals for the new offspring
    if (n_offspring > 0) {
      # generate times
      new_serials <- serial(n_offspring)
      if (any(new_serials < 0)) {
        # stop("Serial interval must be >= 0.")
        cli::cli_abort(c(
          "Serial interval must be >= 0.",
          "x" = "Current specified {.arg serial} is generating
          serial interval as {length(new_serials)}.",
        "i" = "Check the specification of {.arg serial}. For more information
        on specifying {.arg serial}, see {.code ?chain_sim}"
        ))
      }

      # create data frame for new cases
      new_ids <- seq_len(n_offspring) + max(tdf$id)

      new_df <- data.frame(
        id = new_ids,
        time = new_serials + earliest_unsim_time,
        ancestor = tdf$id[idx],
        generation = tdf$generation[idx] + 1L,
        offspring_generated = FALSE
      )

      # add new cases to data frame
      tdf <- rbind(tdf, new_df)
    }

    # update susceptible count
    susc <- susc - n_offspring
  }

  ## remove cases with time > tf that could
  ## have been generated in the last generation
  tdf <- tdf[tdf$time <= tf, ]

  ## sort output and remove columns not needed
  tdf <- tdf[order(tdf$time, tdf$id), ]
  tdf$offspring_generated <- NULL

  return(tdf)
}
