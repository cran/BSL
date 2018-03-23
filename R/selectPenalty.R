#' Selecting BSLasso Penalty
#' 
#' @description This is the main function for selecting the penalty for BSLasso based on a point estimate of the parameters. 
#' Parallel computing is supported with the R package \code{foreach}.
#' 
#' @param ssy                   A summary statistic vector for the observed data.
#' @param n                     A vector of possible values of \code{n}, the number of simulations from the model per MCMC iteration, to test.
#' @param lambda_all            A list, with each entry containing the vector of penalty values to test for the corresponding choice of \code{n}.
#' @param theta                 A point estimate of the parameter value which all of the simulations will be based on.
#' @param M                     The number of repeats to use in estimating the standard deviation of the estimated log synthetic likelihood.
#' @param sigma                 The standard deviation of the log synthetic likelihood to aim for.
#' @param fn_sim                A function that simulates data for a given parameter value. The function can have at most two arguments.
#' The first is the vector of tuning parameters and the second (optional) argument is a list of other necessary arguments. See \code{sim_options}.
#' @param fn_sum                A function for computing summary statistics of data. The function can have at most two arguments.
#' The first is the vector of tuning parameters and the second (optional) argument is a list of other necessary arguments. See \code{sum_options}.
#' @param sim_options           A list of additional arguments to pass into the simulation function. Only use when the input \code{fn_sim} requires additional arguments. The default is \code{NULL}.
#' @param sum_options           A list of additional arguments to pass into the summary statistic function. Only use when the input \code{fn_sum} requires additional arguments. The default is \code{NULL}.
#' @param standardise           A logical argument that determines whether to standardise the summary statistics before applying the graphical lasso. This is only valid if penalty is not \code{NULL}.
#' The diagonal elements will not be penalised if the penalty is not \code{NULL}. The default is \code{FALSE}.
#' @param parallel_sim          A logical value indicating whether parallel computing should be used for simulation and summary statistic evaluation. Default is \code{FALSE}.
#' @param parallel_sim_packages A character vector of package names to pass into the \code{foreach} function as argument '.package'. Only used when \code{parallel_sim} is \code{TRUE}, default is \code{NULL}.
#' @param parallel_main         A logical value indicating whether parallel computing should be used to computing the graphical lasso function. Default is \code{FALSE}.
#' @param verbose               A logical argument indicating whether the iteration numbers (\code{1:M}) should be printed to track progress. The default is \code{FALSE}.
#'
#' @return 				An object of class \code{penbsl} is returned, containing the following components:
#' \itemize{
#' \item \code{resultsDF}: A data frame containing the following:
#'    \itemize{
#'    \item \code{n}: The choices of \code{n} that were specified.
#'    \item \code{penalty}: The choices of the penalty that were specified.
#'    \item \code{sigma}: The standard deviation of the log synthetic likelihood under the above choices.
#'    \item \code{sigmaOpt}: An indicator of whether it was the closest \code{sigma} to the desired one for each choice of \code{n}.
#'    }
#' \item \code{call}: The original code that was used to call the method.
#' }
#' The functions print() and plot() are both available for types of class \code{penbsl}.
#'
#' @references 
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018). Accelerating
#' Bayesian synthetic likelihood with the graphical lasso. https://eprints.qut.edu.au/102263/
#' 
#' @author 								Ziwen An, Christopher C. Drovandi and Leah F. South
#' @seealso 							\code{\link{bsl}} for a function to run BSLasso after selecting the tuning parameter
#' and \code{\link{plot}} for functions related to visualisation.
#' @export
selectPenalty <- function(ssy, n, lambda_all, theta, M, sigma, fn_sim, fn_sum, sim_options = NULL, 
                          sum_options = NULL, standardise = FALSE, parallel_sim = FALSE, parallel_sim_packages = NULL, 
						  parallel_main = FALSE, verbose = TRUE) {
    n <- as.vector(n)
    lambda_all <- as.list(lambda_all)
    N <- length(n)
    if (length(lambda_all) != N) {
        stop('lambda_all must be a list with the same length as n')
    }
    ns <- length(ssy)
    cl <- match.call()
	
    n_max <- max(n)
    K <- max(sapply(lambda_all, length))
    logSL <- array(NA, c(M, N, K))

    for (m in 1 : M) {

        flush.console()
        if (verbose) cat('m = ', m, '\n')

        # simulate with theta_prop and calculate summaries
        if (!parallel_sim) {
            ssx <- array(0, c(n_max, ns))
            for (j in 1 : n_max) {
                if (is.null(sim_options)) {
                    x <- fn_sim(theta)
                } else {
                    x <- fn_sim(theta, sim_options)
                }
                if (is.null(sum_options)) {
                    ssx[j, ] <- fn_sum(x)
                } else {
                    ssx[j, ] <- fn_sum(x, sum_options)
                }
            }
        } else { # use foreach for parallel computing
            ssx <- foreach (j = 1 : n_max, .combine = rbind, .packages = parallel_sim_packages) %dopar% {
                if (is.null(sim_options)) {
                    x <- fn_sim(theta)
                } else {
                    x <- fn_sim(theta, sim_options)
                }
                if (is.null(sum_options)) {
                    fn_sum(x)
                } else {
                    fn_sum(x, sum_options)
                }
            }
        }

        for (i in 1 : N) {
            n_curr <- n[i]
            ssx_curr <- ssx[sample(n_max, n_curr), ]
            mu <- colMeans(ssx_curr)
            S <- cov(ssx_curr)
            if (standardise) {
                std <- apply(ssx_curr, MARGIN = 2, FUN = sd)
                ssx_curr_std <- (ssx_curr - matrix(mu, n_curr, ns, byrow = TRUE)) / matrix(std, n_curr, ns, byrow = TRUE)
				cov_ssx_curr_std <- cov(ssx_curr_std)
            }
            
            if (!parallel_main) {
                for (k in 1 : K) {
				    if (is.na(lambda_all[[i]][k])) {
				        next
				    }
                    if (standardise) {
                        gl <- glasso(cov_ssx_curr_std, rho = lambda_all[[i]][k])
                        corr <- gl$w
                        Sigma <- std %*% t(std) * corr
                        Omega <- solve(Sigma)
                    } else {
                        gl <- glasso(S, rho = lambda_all[[i]][k])
                        Sigma <- gl$w
                        Omega <- gl$wi
                    }
                    logSL[m, i, k] <- as.numeric(-0.5 * log(det(Sigma)) - 0.5 * t(ssy-mu) %*% Omega %*% (ssy-mu))
                }
            } else {
                logSL[m, i, ] <- foreach (k = 1 : K, .combine = rbind, .packages = c('glasso', 'cvTools'), .export='glasso_cv') %dopar% {
				    if (is.na(lambda_all[[i]][k])) {
				        return(NA)
				    }
                    if (standardise) {
                        gl <- glasso(cov_ssx_curr_std, rho = lambda_all[[i]][k])
                        corr <- gl$w
                        Sigma <- std %*% t(std) * corr
                        Omega <- solve(Sigma)
                    } else {
                        gl <- glasso(S, rho = lambda_all[[i]][k])
                        Sigma <- gl$w
                        Omega <- gl$wi
                    }
                    as.numeric(-0.5 * log(det(Sigma)) - 0.5 * t(ssy-mu) %*% Omega %*% (ssy-mu))
                }
            }
            
        }
    }

    resultsDF <- array(list(), N)
    for (i in 1 : N) {
        temp <- apply(logSL[, i, ], MARGIN = 2, FUN = sd)
        resultsDF[[i]] <- data.frame(n = n[i], penalty = lambda_all[[i]], 
                              sigma = temp[!is.na(temp)], sigmaOpt = FALSE)
        resultsDF[[i]]$sigmaOpt[which.min(abs(resultsDF[[i]]$sigma - sigma))] <- TRUE
    }
    resultsDF <- do.call(rbind, resultsDF)
	
    results <- list(resultsDF = resultsDF, call = cl)
    class(results) <- 'penbsl'
    
    return(results)
}
