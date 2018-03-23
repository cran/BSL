#' Performing BSL and BSLasso
#'
#' @description This is the main function for performing MCMC BSL and MCMC BSLasso.
#' Parallel computing is supported with the R package \code{foreach}.
#'
#' @param y				The observed data - note this should be the raw dataset NOT the set of summary statistics.
#' @param n				The number of simulations from the model per MCMC iteration.
#' @param M				The number of MCMC iterations.
#' @param start			Initial guess of the parameter value.
#' @param cov_rw		A covariance matrix to be used in multivariate normal random walk proposals.
#' @param fn_sim        A function that simulates data for a given parameter value. The function can have at most two arguments.
#' The first should be the vector of tuning parameters and the second (optional) argument should be a list of other necessary arguments. See \code{sim_options}.
#' @param fn_sum        A function for computing summary statistics of data. The function can have at most two arguments.
#' The first should be the vector of tuning parameters and the second (optional) argument should be a list of other necessary arguments. See \code{sum_options}.
#' @param penalty		The penalty value of the graphical lasso. BSLasso is performed if a numeric value is specified. The default is \code{NULL}, which means that BSL is performed.
#' @param fn_prior		A function that computes the prior density for a parameter. The default is \code{NULL}, which is an improper flat prior over the real line for each parameter. The function must have a single input: the vector of proposed parameters.
#' @param sim_options	A list of additional arguments to pass into the simulation function. Only use when the input \code{fn_sim} requires additional arguments. The default is \code{NULL}.
#' @param sum_options	A list of additional arguments to pass into the summary statistics function. Only use when the input \code{fn_sum} requires additional arguments. The default is \code{NULL}.
#' @param standardise	A logical argument that determines whether to standardise the summary statistics before applying the graphical lasso. This is only valid if penalty is not \code{NULL}.
#' The diagonal elements will not be penalised if the penalty is not \code{NULL}. The default is \code{FALSE}.
#' @param parallel		A logical value indicating whether parallel computing should be used for simulation and summary statistic evaluation. Default is \code{FALSE}.
#' @param parallel_packages	A character vector of package names to pass into the \code{foreach} function as argument '.package'. Only used when parallel computing is enabled, default is \code{NULL}.
#' @param theta_names	A character vector of parameter names, which must has the same length as the parameter vector. The default is \code{NULL}.
#' @param verbose       A logical argument indicating whether the iteration numbers (\code{1:M}) and accepted proposal flags should be printed to track progress. The default is \code{FALSE}.

#' @return 				An object of class \code{bsl} is returned, containing the following components:
#' \itemize{
#' \item \code{theta}: MCMC samples from the joint approximate posterior distribution of the parameters.
#' \item \code{loglike}: MCMC samples of the estimated log-likelihood values.
#' \item \code{acceptanceRate}: The acceptance rate of the MCMC algorithm.
#' \item \code{earlyRejectionRate}: The early rejection rate of the algorithm (early rejection may occur when using bounded prior distributions).
#' \item \code{call}: The original code that was used to call the method.
#' \item \code{theta_names}: A character vector of parameter names.
#' }
#' The functions print(), summary() and plot() are all available for types of class \code{bsl}. Multiple results
#' can be plotted with overlapping densities using \code{\link{combinePlotsBSL}}.
#'
#' @references
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018). Accelerating
#' Bayesian synthetic likelihood with the graphical lasso. \url{https://eprints.qut.edu.au/102263/}
#' 
#' Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018).
#' Bayesian synthetic likelihood. To appear in Journal of Computational and Graphical Statistics.
#' \url{https://eprints.qut.edu.au/92795/}
#'
#' @author 								Ziwen An, Christopher C. Drovandi and Leah F. South
#' @seealso 							\code{\link{selectPenalty}} for a function to tune the BSLasso tuning parameter
#' and \code{\link{plot}} for functions related to visualisation.
#' @export
bsl <- function(y, n, M, start, cov_rw, fn_sim, fn_sum, penalty = NULL, 
    fn_prior = NULL, sim_options = NULL, sum_options = NULL, standardise = FALSE, 
    parallel = FALSE, parallel_packages = NULL, theta_names = NULL, verbose = TRUE) {
    if (!parallel & !is.null(parallel_packages)) {
        warning('parallel_packages is not used in serial computing')
    }
    if (is.null(penalty) & standardise) {
        warning('standardisation is only supported in BSLasso')
    }
    
    cl <- match.call()
    k <- length(start)
    if (is.null(sum_options)) {
        ssy <- fn_sum(y)
    } else {
        ssy <- fn_sum(y, sum_options)
    }
    ns <- length(ssy)
    theta_curr <- start
    theta <- array(0, c(M, k), dimnames = list(NULL, theta_names))
    loglike <- numeric(M)
    countAcc <- 0
    countEr <- 0

    # simulate with theta_prop and calculate the summary statistics
    if (!parallel) {
        ssx <- array(0, c(n, ns))
        for (j in 1 : n) {
            if (is.null(sim_options)) {
                x <- fn_sim(theta_curr)
            } else {
                x <- fn_sim(theta_curr, sim_options)
            }
            if (is.null(sum_options)) {
                ssx[j, ] <- fn_sum(x)
            } else {
                ssx[j, ] <- fn_sum(x, sum_options)
            }
        }
    } else { # use foreach for parallel computing
        ssx <- foreach (j = 1 : n, .combine = rbind, .packages = parallel_packages) %dopar% {
            if (is.null(sim_options)) {
                x <- fn_sim(theta_curr)
            } else {
                x <- fn_sim(theta_curr, sim_options)
            }
            if (is.null(sum_options)) {
                fn_sum(x)
            } else {
                fn_sum(x, sum_options)
            }
        }
    }

    # compute the loglikelihood
    if (is.null(penalty)) { # BSL if penalty is not provided
        mu <- colMeans(ssx)
        Sigma <- cov(ssx)
        Omega <- solve(Sigma)
    } else { # BSLasso if a penalty value is detected
        mu <- colMeans(ssx)
        S <- cov(ssx)
        if (!standardise) { # use graphical lasso without standardisation
            gl <- glasso(S, rho = penalty)
            Sigma <- gl$w
            Omega <- gl$wi
        } else { # standardise the summary statistics before passing into the graphical lasso function
            std <- apply(ssx, MARGIN = 2, FUN = sd)
            ssx_std <- (ssx - matrix(mu, n, ns, byrow = TRUE)) / matrix(std, n, ns, byrow = TRUE)
            gl <- glasso(cov(ssx_std), rho = penalty, penalize.diagonal = FALSE) # do not penalise the diagonal entries since we want the correlation matrix
            corr <- gl$w
            Sigma <- std %*% t(std) * corr
            Omega <- solve(Sigma)
        }
    }
    loglike_curr <- -0.5 * log(det(Sigma)) - 0.5 * t(ssy-mu) %*% Omega %*% (ssy-mu)

    for (i in 1:M) {

        flush.console()
        if (verbose) cat('i = ', i, '\n')
		
		# multivariate normal random walk to the proposed value of theta
        theta_prop <- mvrnorm(1, theta_curr, cov_rw)
        
        # early rejection if the proposed theta falls outside of prior coverage / feasible region
        if (!is.null(fn_prior)) {
            p = fn_prior(theta_prop) / fn_prior(theta_curr)
            if (p == 0) {
                if (verbose) cat('*** early rejection ***\n')
                theta[i, ] <- theta_curr
                loglike[i] <- loglike_curr
                countEr <- countEr + 1
                next
            }
        } else {
            p <- 1
        }

        # simulate with theta_prop and calculate the summary statistics
        if (!parallel) {
            ssx <- array(0, c(n, ns))
            for (j in 1 : n) {
                if (is.null(sim_options)) {
                    x <- fn_sim(theta_prop)
                } else {
                    x <- fn_sim(theta_prop, sim_options)
                }
                if (is.null(sum_options)) {
                    ssx[j, ] <- fn_sum(x)
                } else {
                    ssx[j, ] <- fn_sum(x, sum_options)
                }
            }
        } else { # use foreach for parallel computing
            ssx <- foreach (j = 1 : n, .combine = rbind, .packages = parallel_packages) %dopar% {
                if (is.null(sim_options)) {
                    x <- fn_sim(theta_prop)
                } else {
                    x <- fn_sim(theta_prop, sim_options)
                }
                if (is.null(sum_options)) {
                    fn_sum(x)
                } else {
                    fn_sum(x, sum_options)
                }
            }
        }

        # compute the loglikelihood
        if (is.null(penalty)) { # BSL if penalty is not provided
            mu <- colMeans(ssx)
            Sigma <- cov(ssx)
            temp <- try(solve(Sigma))
            if (class(temp) == 'try-error') {
                theta[i, ] <- theta_curr
                loglike[i] <- loglike_curr
                next
            } else {
                Omega <- temp
            }
        } else { # BSLasso if a penalty value is detected
            mu <- colMeans(ssx)
            S <- cov(ssx)
            if (!standardise) { # use graphical lasso without standardisation
                gl <- glasso(S, rho = penalty)
                Sigma <- gl$w
                Omega <- gl$wi
            } else { # standardise the summary statistics before passing into the graphical lasso function
                std <- apply(ssx, MARGIN = 2, FUN = sd)
                ssx_std <- (ssx - matrix(mu, n, ns, byrow = TRUE)) / matrix(std, n, ns, byrow = TRUE)
                gl <- glasso(cov(ssx_std), rho = penalty, penalize.diagonal = FALSE) # do not penalise the diagonal entries since we want the correlation matrix
                corr <- gl$w
                Sigma <- std %*% t(std) * corr
                Omega <- solve(Sigma)
            }
        }
        loglike_prop <- -0.5 * log(det(Sigma)) - 0.5 * t(ssy-mu) %*% Omega %*% (ssy-mu)

		# accept the proposed theta with a probability
        if (runif(1) < p * exp(loglike_prop - loglike_curr)) {
            if (verbose) cat('*** accept ***\n')
            theta_curr <- theta_prop
            loglike_curr <- loglike_prop
            countAcc <- countAcc + 1
        }

        theta[i, ] <- theta_curr
        loglike[i] <- loglike_curr
    }

    accRate <- countAcc / M
    erRate <- countEr / M

    results <- list(theta = theta, loglike = loglike, acceptanceRate = accRate, earlyRejectionRate = erRate, 
                    call = cl, theta_names = theta_names)
    class(results) <- 'bsl'
    return(results)
}
