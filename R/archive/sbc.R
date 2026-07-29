sbc <- function(model, generator, N_steps, iter = 2000,
                cores = getOption("mc.cores", parallel::detectCores()),
                cache_dir = NULL,
                ...) {
  true_list <- list()
  observed_list <- list()
  args_per_fit <- list()
  for(i in 1:N_steps) {
    data <- generator()
    true_list[[i]] <- data$true
    args_per_fit[[i]] <- list(data = data$observed)
    observed_list[[i]] <- data$observed
  }

  convert_to_rstan <- function(x) {
    if(inherits(x, "stanfit")) {
      x
    } else if(inherits(x, "CmdStanMCMC")) {
      rstan::read_stan_csv(x$output_files())
    } else {
      stop("Invalid object")
    }
  }


  fits <- sampling_parallel(args_shared = list(model = model, iter = iter, ...), args_per_fit = args_per_fit,
                         cores = cores, cache_dir = cache_dir, cache_fits = FALSE,
                         cache_summaries = !is.null(cache_dir),
                         summarise_fun = convert_to_rstan)

  param_stats <-
    fits %>% imap_dfr(function(fit, data_id) {
      samples <- rstan::extract(fit) #rstan::extract(fit, pars = names(true_list[[data_id]]))
      eval <- evaluate_all_params(samples, true_list[[data_id]])
      eval %>% mutate(run = data_id)
    })
  diagnostics <-
    fits %>% imap_dfr(function(fit, run_id) {
      data.frame(run = run_id,
                 n_divergent = rstan::get_num_divergent(fit),
                 n_treedepth = rstan::get_num_max_treedepth(fit),
                 n_chains_low_bfmi = length(rstan::get_low_bfmi_chains(fit)),
                 max_rhat = max(rstan::summary(fit)$summary[,"Rhat"], na.rm = TRUE),
                 min_rel_neff = min(rstan::summary(fit)$summary[,"n_eff"], na.rm = TRUE) / iter,
                 total_time = sum(rstan::get_elapsed_time(fit))
                 )
    })

  return(list(params = param_stats, diagnostics = diagnostics, data = observed_list, true_values = true_list))
}

plot_sbc_params <- function(params, binwidth = 10, caption = NULL, plot_stat = "mean", x_axis_trans = "identity", y_axis_trans = "identity") {
  if(!plot_stat %in% c("median","mean") ){
    stop("Invalid plot_stat")
  }

  #CI - taken from https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R

  if(100 %% binwidth != 0) {
    stop("binwidth has to divide 100")
  }
  n_simulations <- length(unique(params$run))
  CI = qbinom(c(0.005,0.5,0.995), size=n_simulations,prob  =  binwidth / 100)
  ci_lower = CI[1]
  ci_mean = CI[2]
  ci_upper = CI[3]

  #The visualisation style taken as well from   https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R
  print(params %>%
          ggplot(aes(x = order_within)) +
          geom_segment(aes(x=0,y=ci_mean,xend=100,yend=ci_mean),colour="grey25") +
          geom_polygon(data=data.frame(x=c(-10,0,-10,110,100,110,-10),y=c(ci_lower,ci_mean,ci_upper,ci_upper,ci_mean,ci_lower,ci_lower)),aes(x=x,y=y),fill="grey45",color="grey25",alpha=0.5) +
          geom_histogram(breaks =  seq(1, 101, by = binwidth), closed = "left" ,fill="#A25050",colour="black") +
          facet_wrap(~param_name, scales = "free_y") +
          ggtitle("Posterior order within 100 samples")
  )

  point_alpha <- 1 / ((n_simulations * 0.03) + 1)
  print(params %>%
          ggplot(aes_string(x = "true_value", y = plot_stat)) + geom_point(alpha = point_alpha) +
          geom_abline(slope = 1, intercept = 0, color = "blue") +
          facet_wrap(~param_name, scales = "free")  +
          scale_x_continuous(trans = x_axis_trans) +
          scale_y_continuous(trans = y_axis_trans) +
          ggtitle(paste0(plot_stat, " of marginal posteriors vs. true value - ", caption))
  )
}

summarise_sbc_diagnostics <- function(sbc_results) {
  sbc_results$diagnostics %>%
    summarise(
      has_divergence = mean(n_divergent > 0),
      has_treedepth = mean(n_treedepth > 0),
      has_low_bfmi = mean(n_chains_low_bfmi > 0),
      median_total_time = median(total_time),
      has_high_rhat = mean(max_rhat > 1.01),
      has_low_neff = mean(min_rel_neff < 0.04)
      )

}
