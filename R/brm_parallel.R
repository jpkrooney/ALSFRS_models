sampling_multi2 <- function(args_shared, args_per_fit,
                         cores = getOption("mc.cores", 1),
                         summarise_fun = NULL,
                         summarise_fun_dependencies = c(),
                         cache_dir = NULL,
                         R_session_init_expr = NULL) {

    if(!is.list(args_shared)) {
        stop("args_shared must be a list")
    }
    if(!is.list(args_per_fit) || length(args_per_fit) <= 0) {
        stop("args_per_fit must be a non-empty list")
    }

    if(!is.null(cache_dir) && !dir.exists(cache_dir)) {
        stop(paste0("Cache dir '", cache_dir,"'  does not exist"))
    }

    n_fits <- length(args_per_fit)

    if("cores" %in% names(args_shared) || "num_cores" %in% names(args_shared)) {
        stop("args_shared must not specify cores or num_cores")
    }

    uses_rstan <- FALSE
    uses_cmdstan <- FALSE
    model_in_shared_args <- FALSE
    data_in_shared_args <- "data" %in% names(args_shared)

    if("model" %in% names(args_shared)) {
        if(inherits(args_shared$model, "stanmodel")) {
            uses_rstan <- TRUE
        } else if(inherits(args_shared$model, "CmdStanModel")) {
            uses_cmdstan <- TRUE
        } else {
            stop("Model in shared args is not of class 'stanmodel' or 'CmdStanModel'")
        }
        model_in_shared_args <- TRUE
    }

    for(i in 1:n_fits) {
        if(!is.list(args_per_fit[[i]])) {
            stop("All elements of args_per_fit have to be lists")
        }

        if(length(intersect(names(args_shared), names(args_per_fit[[i]]))) > 0) {
            stop(paste0("No parameters provided in args_per_fit can be given in args_shared.\n
                 Found intersection at index ", i, "."))
        }

        if("model" %in% names(args_per_fit[[i]])) {
            if(inherits(args_per_fit[[i]]$model, "stanmodel")) {
                uses_rstan <- TRUE
            } else if(inherits(args_per_fit[[i]]$model, "CmdStanModel")) {
                uses_cmdstan <- TRUE
            } else {
                stop(paste0("Model for fit id ", i," is not of class 'stanmodel' or 'CmdStanModel'"))
            }
        } else if(!model_in_shared_args) {
            stop(paste0("No model argument in shared_args and fit id ", i, " does not provide model"))
        }

        if(!data_in_shared_args && !("data" %in% names(args_per_fit[[i]]))) {
            stop(paste0("No data argument in shared_args and fit id ", i, " does not provide data"))
        }


        if("cores" %in% names(args_per_fit[[i]]) || "num_cores" %in% names(args_per_fit[[i]])) {
            stop(paste0("args_per_fit[[", i, "]] must not specify cores or num_cores"))
        }
    }

    if(2 * n_fits <= cores) {
        cores_per_fit <- floor(cores / n_fits)
        cluster_cores <- n_fits
    } else {
        cores_per_fit <- 1
        cluster_cores <- min(c(cores, n_fits))
    }

    cl <- parallel::makeCluster(cluster_cores, useXDR = FALSE)
    on.exit(parallel::stopCluster(cl))


    fit_fun <- function(args, args_shared, summarise_fun, cores_per_fit,
                        cache_dir, cmdstan_fit_dir) {
        all_args <- c(args_shared, args)
        all_args$cores <- cores_per_fit

        model <- all_args$model

        if(inherits(model, "stanmodel")) {
            model_code <- model@model_code
            is_rstan <- TRUE
        } else if(inherits(model, "CmdStanModel")) {
            model_code <- model$code()
            is_rstan <- FALSE
        } else {
            stop("Invalid model")
        }


        cached <-  FALSE
        not_cached_msg <- ""
        if(!is.null(cache_dir)) {
            data <- all_args$data

            code_hash <- rlang::hash(model_code)
            data_hash <- rlang::hash(data)


            cache_file <- paste0(cache_dir, "/", code_hash, "_", data_hash, ".rds")
            if(file.exists(cache_file)) {
                fit_from_file <- readRDS(cache_file)
                if((is_rstan && inherits(fit_from_file, "stanmodel"))
                   || (!is_rstan && inherits(fit_from_file, "CmdStanModel"))
                   ) {
                    fit <- fit_from_file
                    cached <- TRUE
                }
            }
        }

        if(!cached) {
            if(inherits(model,"stanmodel")) {
                fit <- do.call(rstan::sampling, args = all_args)
                if(!is.null(cache_dir)) {
                    saveRDS(fit, cache_file)
                }
            } else {
                translated_args <- list()
                for(old in names(all_args)) {
                    if(old == "chains") {
                        translated_args$num_chains = all_args$chains
                    } else if(old == "cores") {
                        translated_args$num_cores = all_args$cores
                    } else if(old == "control") {
                        if(!is.null(all_args$control$adapt_delta)) {
                            translated_args$adapt_delta = all_args$control$adapt_delta
                        }
                        if(!is.null(all_args$control$max_treedepth)) {
                            translated_args$max_depth = all_args$control$max_treedepth
                        }
                    } else if(old == "iter") {
                        translated_args$num_warmup = all_args$iter / 2
                        translated_args$num_samples = all_args$iter/ 2
                    } else {
                        translated_args[[old]] = all_args[[old]]
                    }
                }
                translated_args$model <- NULL
                fit <- do.call(model$sample, args = translated_args)
                fit$save_output_files(cmdstan_fit_dir)
                if(!is.null(cache_dir)) {
                    fit$save_object(cache_file)
                }
            }
        }

        if(!is.null(summarise_fun)) {
            summarise_fun(fit)
        } else {
            fit
        }
    }


    dependencies <- c()
    if(uses_rstan) {
        dependencies <- c(dependencies, "rstan", "Rcpp")
    }
    if(uses_cmdstan) {
        dependencies <- c(dependencies, "cmdstanr")
    }
    dependencies <- c(dependencies, summarise_fun_dependencies)
    .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
        dirname(system.file(package = d))
    })))
    .paths <- .paths[.paths != ""]
    parallel::clusterExport(cl, varlist = c(".paths","dependencies"), envir = environment())
    parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
    parallel::clusterEvalQ(cl, expr =
                               for(dep in dependencies) {
                                   suppressPackageStartupMessages(require(dep, quietly = TRUE, character.only = TRUE))
                               }
    )

    #    parallel::clusterExport(cl, varlist = "args_shared", envir = environment())

    parallel::clusterExport(cl, varlist =
                                c("R_session_init_expr"),
                            envir = environment())
    parallel::clusterEvalQ(cl, expr = R_session_init_expr)

    results <- parallel::parLapplyLB(cl, X = args_per_fit,
                                     fun = fit_fun,
                                     chunk.size = 1,
                                     args_shared = args_shared,
                                     summarise_fun = summarise_fun,
                                     cores_per_fit = cores_per_fit,
                                     cache_dir = cache_dir,
                                     cmdstan_fit_dir = tempdir())

    results
}




brm_parallel <- function(args_shared, args_per_fit,
                         cores = options("mc.cores"),
                         summarise_fun = NULL,
                         summarise_fun_dependencies = c(),
                         backend = options("brms.backend"),
                         pre_compile = backend == "cmdstanr",
                         R_session_init_expr = NULL) {

    if(!is.list(args_shared)) {
        stop("args_shared must be a list")
        formulas <- list(formulas)
    }
    if(!is.list(args_per_fit) || length(args_per_fit) <= 0) {
        stop("args_per_fit must be a non-empty list")
    }

    n_fits <- length(args_per_fit)

    if("cores" %in% names(args_shared)) {
        stop("args_shared must not specify cores")
    }

    for(i in 1:n_fits) {
        if(!is.list(args_per_fit[[i]])) {
            stop("All elements of args_per_fit have to be lists")
        }

        if(length(intersect(names(args_shared), names(args_per_fit[[i]]))) > 0) {
            stop(paste0("No parameters provided in args_per_fit can be given in args_shared.\n
                 Found intersection at index ", i, "."))
        }

        if("cores" %in% names(args_per_fit[[i]])) {
            stop(paste0("args_per_fit[[", i, "]] must not specify cores"))
        }

        if(pre_compile) {
            all_args <- c(args_shared, args_per_fit[[i]])
            all_args$fit <- FALSE
        }
    }

    if(2 * n_fits <= cores) {
        cores_per_fit <- floor(cores / n_fits)
        cluster_cores <- n_fits
    } else {
        cores_per_fit <- 1
        cluster_cores <- cores
    }

    cl <- parallel::makeCluster(cluster_cores, useXDR = FALSE)
    on.exit(parallel::stopCluster(cl))


    fit_fun <- function(args, args_shared, summarise_fun, cores_per_fit) {
        all_args <- c(args_shared, args)
        all_args$cores <- cores_per_fit
        fit <- do.call(brm, args = all_args)
        if(!is.null(summarise_fun)) {
            summarise_fun(fit)
        } else {
            fit
        }
    }



    dependencies <- c("brms", summarise_fun_dependencies)
    .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
        dirname(system.file(package = d))
    })))
    .paths <- .paths[.paths != ""]
    parallel::clusterExport(cl, varlist = c(".paths","dependencies"), envir = environment())
    parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
    parallel::clusterEvalQ(cl, expr =
                               for(dep in dependencies) {
                                   suppressPackageStartupMessages(require(dep, quietly = TRUE, character.only = TRUE))
                               }
    )

#    parallel::clusterExport(cl, varlist = "args_shared", envir = environment())

    parallel::clusterExport(cl, varlist =
                                c("R_session_init_expr"),
                            envir = environment())
    parallel::clusterEvalQ(cl, expr = R_session_init_expr)

    results <- parallel::parLapplyLB(cl, X = args_per_fit,
                                     args_shared = args_shared,
                                     summarise_fun = summarise_fun,
                                     cores_per_fit = cores_per_fit,
                                     fun = fit_fun,
                                     chunk.size = 1)

    results
}
