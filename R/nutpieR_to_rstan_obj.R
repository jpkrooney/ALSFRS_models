
# function not needed anymore ad now built into nutpieR
#trim_pars <- function(nutfit, pars, exclude = TRUE){
#    allpars <- attributes(nutfit)$dimnames[[3]]
#    dims <- attributes(nutfit)$dim
#
#    idx <- lapply(1:length(pars), function(i) {
#        pat <- paste("^", pars[i], "(\\[(\\d+,)*\\d+\\])*$", sep = "")
#        grepl(pat, allpars) }) |>
#        bind_cols() |>
#        rowMeans() > 0
#
#    # If exclude = TRUE invert index
#    if (exclude == TRUE) { idx <- !idx}
#
#    # Apply the exclusions to nutfit
#    trimfit <- nutfit[ , , idx ]
#
#    attributes(trimfit)$diagnostics <- attributes(nutfit)$diagnostics
#    attributes(trimfit)$num_chains <- attributes(nutfit)$num_chains
#    attributes(trimfit)$num_warmup <- attributes(nutfit)$num_warmup
#    attributes(trimfit)$num_draws <- attributes(nutfit)$num_draws
#
#    return(trimfit)
#}








make_sim_object <- function(nutfit) {

    num_warmup <- attributes(nutfit)$num_warmup
    num_draws <- attributes(nutfit)$num_draws

    # Get diagnostics from nutfit
    diagnostics <- nutpie_diagnostics(nutfit)

    # add lp__ and repackage samples
    lp__ <- array(diagnostics$logp, dim = c(dim(nutfit)[1], dim(nutfit)[2], 1),
                  dimnames = c(dimnames(nutfit)[1:2],
                               list("lp__")))
    samples <- posterior::as_draws_list(abind::abind(nutfit, lp__, along = 3))
    names(samples) <- NULL
    class(samples) <- NULL

    # rest of code adapted from rstan::read_stan_csv()
    n <- length(diagnostics$diverging)
    chains <- attr(diagnostics, "num_chains") %||% 1L
    n_per_chain <- n %/% chains

    fnames <- names(samples[[1]])
    fnames <- rstan:::sqrfnames_to_dotfnames(fnames)
    n_save <- nrow(samples[[1]])
    paridx <- rstan:::paridx_fun(fnames)
    lp__idx <- attr(paridx, "meta")["lp__"]
    par_fnames <- c(fnames[paridx], "lp__")
    pars_oi <- rstan:::unique_par(par_fnames)
    dims_oi <- lapply(pars_oi, function(i) {
        pat <- paste("^", i, "(\\.\\d+)*$", sep = "")
        i_fnames <- par_fnames[grepl(pat, par_fnames)]
        rstan:::get_dims_from_fnames(i_fnames, i)
    })
    names(dims_oi) <- pars_oi

    sim <- list(samples = samples,
               iter = num_warmup + num_draws,
               thin = 1,
               warmup = num_warmup,
               chains = chains,
               n_save = rep(num_draws, chains),
               #warmup2 = rep(attributes(nutfit)[[1]][[1]], chains), # not sure about this line
               warmup2 = rep(0, chains), # not sure about this line
               permutation = lapply(1:chains, function(id) 1:num_draws),
               pars_oi = pars_oi,
               dims_oi = dims_oi,
               fnames_oi = rstan:::dotfnames_to_sqrfnames(par_fnames),
               n_flatnames = length(par_fnames))

    return(sim)
}



nutpieR_to_rstan_obj <- function(nutfit, stanmodel){
    # make sim object
    sim <- make_sim_object(nutfit)

    # make a stan_args list
    stan_args <- lapply(1 : sim$chains, function(i) {
        list(chain_id = i, iter = sim$iter, warmup = sim$warmup,
        thin = 1, seed = NULL, algorithm = "NUTS", engine = "nutpieR") } )

    # bundle into stanfit object
    stanfit <- new("stanfit",
                model_name = "",
                model_pars = sim$pars_oi,
                par_dims = sim$dims_oi,
                mode = as.integer(0),
                sim = sim,
                inits = list(),
                stan_args = stan_args,
                stanmodel = stanmodel,
                date = "")
    ## add divergences etc ??

    return(stanfit)
}



nutpie_to_brms <- function(form, dat, model, fit, family, stanvars = NULL){
    emptybrms <- brm(bf(form, family = family) +
                         set_rescor(FALSE), stanvars = stanvars, adapt_delta = 0.95, init = 0.1,
                     data = dat, empty = TRUE)
    #
    rstanobj <- nutpieR_to_rstan_obj(fit, model)
    brmsfit <- emptybrms
    brmsfit$fit <- rstanobj
    brmsfit <- rename_pars(brmsfit)
    brmsfit
}



# fragile hacky function to make a nutpieR fit object into a stanfit object

#nutpieR_to_rstan_obj <- function(nutfit, emptyrstan, num_draws, num_warmup){
#
#    attribs <- attributes(nutfit)
#    diagnostics <- nutpie_diagnostics(nutfit)
#
#    n <- length(diagnostics$diverging)
#    chains <- attr(diagnostics, "num_chains") %||% 1L
#    n_per_chain <- n %/% chains
#    n_save <- nrow(nutfit)
#
#
#    lp__ <- array(diagnostics$logp, dim = c(dim(nutfit)[1], dim(nutfit)[2], 1),
#                  dimnames = c(dimnames(nutfit)[1:2],
#                               list("lp__")))
#    samples <- posterior::as_draws_list(abind::abind(nutfit, lp__, along = 3))
#    names(samples) <- NULL
#    class(samples) <- NULL
#    #for(i in 1:chains){
#    #    attributes(samples[[i]])  <- attributes(emptyrstan@sim$samples[[i]])
#    #}
#
#    # Extract structure
#    fit <- emptyrstan
#
#    # Replace draws
#    fit@sim$samples <- samples
#    fit@sim$iter <- num_draws + num_warmup
#    fit@sim$n_save <- rep(num_draws, chains)
#    fit@sim$warmup <- num_warmup
#    fit@sim$warmup2 <- rep(num_warmup, chains)
#
#    fit@sim$n_save <- rep(num_draws, chains)
#    fit@sim$warmup2 <- attribs[[1]][[1]] - rep(num_warmup, chains)
#
#    # set attributes
#    attributes(fit@sim$samples[[1]])$args$warmup <- as.integer(num_warmup)
#    attributes(fit@sim$samples[[1]])$args$iter <- as.integer(num_draws + num_warmup)
#    #attr(fit@sim$samples[[1]], 'sampler_params')$accept_stat__ <- diagnostics$mean_tree_accept
#    #attr(fit@sim$samples[[1]], 'sampler_params')$stepsize__  <- diagnostics$step_size_bar
#    return(fit)
#}
