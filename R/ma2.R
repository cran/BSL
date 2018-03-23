#' An MA(2) model
#'
#' @description In this example we wish to estimate the parameters of a simple MA(2) time series model.
#' We provide the data and tuning parameters required to reproduce the results in An et al. (2018).
#'
#' @param theta     A vector of proposed model parameters, \eqn{\theta_1} and \eqn{\theta_2}
#' @param options	A list of options for simulating data from the model. For this example, the list contains only \eqn{T}, the number of observations.
#' @param x			Observed or simulated data in the format of a vector of length \eqn{T}.
#'
#' @details
#' This example is based on estimating the parameters of a basic MA(2) time series model
#' of the form
#' 
#' \deqn{y_t = z_t + \theta_1 z_{t-1} + \theta_2 z_{t-2},}
#' 

#' where \eqn{t=1,\ldots,T} and \eqn{z_t ~ N(0,1)} for \eqn{t=-1,0,\ldots,T}.
#' A uniform prior is used for this example, subject to the restrictions that
#' \eqn{-2<\theta_1<2}, \eqn{\theta_1+\theta_2>-1} and \eqn{\theta_1-\theta_2<1}
#' so that invertibility of the time series is satisfied. The summary statistics
#' are simply the full data.
#'
#' @section A simulated dataset:
#'
#' An example 'observed' dataset and the tuning parameters relevant to that example
#' can be obtained using \code{data(ma2)}. This 'observed' data is a simulated dataset
#' with \eqn{\theta_1 = 0.6}, \eqn{\theta_2=0.2} and \eqn{T=50}. Further information
#' about this model and the specific choices of tuning parameters used in BSL and
#' BSLasso can be found in An et al. (2018).
#'
#' \itemize{
#'  \item \code{data}: A time series dataset, in the form of a vector of length \eqn{T}
#'  \item \code{sim_options}: A list containing \eqn{T=50}
#'  \item \code{start}: A vector of suitable initial values of the parameters for MCMC
#'  \item \code{cov}: Covariance matrix of the multivariate normal random walk, in the form of a \eqn{2 x 2} matrix
#' }
#'
#' @examples
#' \dontshow{ 
#' # Loading the data for this example
#' data(ma2)
#' true_ma2 <- c(0.6,0.2)
#' 
#' # Performing BSL
#' resultMa2BSL <- bsl(y = ma2$data, n = 300, M = 5, start = ma2$start, cov_rw = ma2$cov,
#'                  fn_sim = ma2_sim, fn_sum = ma2_sum, fn_prior = ma2_prior,
#'                  sim_options = ma2$sim_options, theta_names = c('theta1', 'theta2'),
#'                  verbose = FALSE)
#' summary(resultMa2BSL)
#' plot(resultMa2BSL, true_value = true_ma2, thin = 1)
#' 
#' # Performing tuning for BSLasso 
#' lambda_all <- list(exp(seq(-3,0.5,length.out=2)), exp(seq(-4,-0.5,length.out=2)))
#' 
#' sp_ma2 <- selectPenalty(ssy = ma2_sum(ma2$data), n = c(50, 150), lambda_all,
#'                  theta = true_ma2, M = 5, sigma = 1.5, fn_sim = ma2_sim,
#'                  fn_sum = ma2_sum, sim_options = ma2$sim_options,
#'                  verbose = FALSE)
#' sp_ma2
#' plot(sp_ma2)
#' 
#' # Performing BSLasso with a fixed penalty
#' resultMa2BSLasso <- bsl(y = ma2$data, n = 150, M = 5, start = ma2$start, cov_rw = ma2$cov,
#'                  fn_sim = ma2_sim, fn_sum = ma2_sum, penalty = 0.027, fn_prior = ma2_prior,
#'                  sim_options = ma2$sim_options, theta_names = c('theta1', 'theta2'),
#'                  verbose = FALSE)
#' summary(resultMa2BSLasso)
#' plot(resultMa2BSLasso, true_value = true_ma2, thin = 1)
#' 
#' # Plotting the results together for comparison
#' combinePlotsBSL(resultMa2BSL, resultMa2BSLasso, true_value = true_ma2, thin = 1)
#' }
#' \dontrun{ 
#' # Loading the data for this example
#' data(ma2)
#' true_ma2 <- c(0.6,0.2)
#' 
#' # Performing BSL
#' resultMa2BSL <- bsl(y = ma2$data, n = 500, M = 300000, start = ma2$start, cov_rw = ma2$cov,
#'                  fn_sim = ma2_sim, fn_sum = ma2_sum, fn_prior = ma2_prior,
#'                  sim_options = ma2$sim_options, theta_names = c('theta1', 'theta2'))
#' summary(resultMa2BSL)
#' plot(resultMa2BSL, true_value = true_ma2, thin = 20)
#' 
#' # Performing tuning for BSLasso 
#' lambda_all <- list(exp(seq(-3,0.5,length.out=20)), exp(seq(-4,-0.5,length.out=20)), 
#'                  exp(seq(-5.5,-1.5,length.out=20)), exp(seq(-7,-2,length.out=20)))
#' 
#' sp_ma2 <- selectPenalty(ssy = ma2_sum(ma2$data), n = c(50, 150, 300, 500), lambda_all,
#'                  theta = true_ma2, M = 100, sigma = 1.5, fn_sim = ma2_sim,
#'                  fn_sum = ma2_sum, sim_options = ma2$sim_options)
#' sp_ma2
#' plot(sp_ma2)
#' 
#' # Performing BSLasso with a fixed penalty
#' resultMa2BSLasso <- bsl(y = ma2$data, n = 300, M = 250000, start = ma2$start, cov_rw = ma2$cov,
#'                  fn_sim = ma2_sim, fn_sum = ma2_sum, penalty = 0.027, fn_prior = ma2_prior,
#'                  sim_options = ma2$sim_options, theta_names = c('theta1', 'theta2'))
#' summary(resultMa2BSLasso)
#' plot(resultMa2BSLasso, true_value = true_ma2, thin = 20)
#' 
#' # Plotting the results together for comparison
#' combinePlotsBSL(resultMa2BSL, resultMa2BSLasso, true_value = true_ma2, thin = 20)
#' }
#' 
#' @references
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018). Accelerating
#' Bayesian synthetic likelihood with the graphical lasso. \url{https://eprints.qut.edu.au/102263/}
#' 
#' @author Ziwen An, Christopher C. Drovandi and Leah F. South
#' 
#' @name ma2
#' @usage data(ma2)
NULL


#' The function \code{ma2_sim(theta,options)} simulates an MA(2) time series.
#' @rdname ma2
ma2_sim <- function(theta, options) {
    T <- options$T
    rand <- rnorm(T + 2)
    y <- rand[3 : (T+2)] + theta[1] * rand[2 : (T+1)] + theta[2] * rand[1 : T]
    return(y)
}

#' The function \code{ma2_sum(x)} returns the summary statistics for a given data set. Since the summary statistics are the data, this function simply returns the data.
#' @rdname ma2
ma2_sum <- function(x) {
    return(x)
}

#' \code{ma2_prior(theta)} evaluates the (unnormalised) prior, which is uniform subject to several restrictions related to invertibility of the time series.
#' @rdname ma2
ma2_prior <- function(theta) {
    theta[2] < 1 & sum(theta) > -1 & diff(theta) > -1
}
