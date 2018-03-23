#' Cell biology example
#'
#' @description This example estimates the probabilities of cell motility and
#' cell proliferation for a discrete-time stochastic model of
#' cell spreading. We provide the data and tuning parameters required to
#' reproduce the results in An et al. (2018).
#'
#' @param theta		    A vector of proposed model parameters, \eqn{Pm} and \eqn{Pp}.
#' @param Y          A \code{rows} \eqn{x} \code{cols} \eqn{x} \code{num_obs} array of the cell presences at times \code{1:num_obs} (not time 0).
#' @param sim_options	A list of options for simulating data from the model. For this example, the list contains
#' \itemize{
#' \item \code{Yinit}: The initial matrix of cell presences of size \code{rows} \eqn{x} \code{cols}.
#' \item \code{rows}: The number of rows in the lattice (rows in the cell location matrix).
#' \item \code{cols}: The number of columns in the lattice (columns in the cell location matrix).
#' \item \code{sim_iters}: The number of discretisation steps to get to when an observation is
#' actually taken. For example, if observations are taken every 5 minutes but the discretisation level is 2.5 minutes,
#' then \code{sim_iters} would be 2. Larger values of \code{sim_iters} lead to more "accurate" simulations from the model, but they also increase the simulation time.
#' \item \code{num_obs}: The total number of images taken after initialisation.
#' }
#' @param sum_options A list of options for simulating data from the model.
#' For this example, the list just contains \code{Yinit}, the same \code{rows} \eqn{x} \code{cols} matrix as in \code{sim_options}.
#
#' @details
#' Cell motility (movement) and proliferation (reproduction) cause
#' tumors to spread and wounds to heal. If we can measure cell proliferation
#' and cell motility under different situations, then we may be able to use
#' this information to determine the efficacy of different medical treatments.
#'
#' A common method for measuring in vitro cell movement and proliferation is
#' the scratch assay. Cells form a layer on an assay and, once
#' they are completely covering the assay, a scratch is
#' made to separate the cells. Images of the cells are taken until the
#' scratch has closed up and the cells are in contact again.
#' Each image can be converted to a binary matrix by forming a lattice
#' and recording the binary matrix (of size \code{rows} \eqn{x} \code{cols}) of cell presences.
#'
#' The model that we consider is a random walk model with parameters for the probability
#' of cell movement (\eqn{Pm}) and the probability of cell proliferation (\code{Pp})
#' and it has no tractable likelihood function. We use the uninformative priors
#' \eqn{\theta_1 ~ U(0,1)} and \eqn{\theta_2 ~ U(0,1)}.
#'
#' We have a total of 145 summary statistics, which are made up of the Hamming distances
#' between the binary matrices for each time point and the total number of cells at the final time.
#' 
#' Details about the types of cells that this model is suitable for
#' and other information can be found in Price et al. (2018) and An et al. (2018). Johnston et al. (2014)
#' use a different ABC method and different summary statistics for a similar example.
#'
#' @section A simulated dataset:
#' 
#' An example 'observed' dataset and the tuning parameters relevant to that example
#' can be obtained using \code{data(cell)}. This 'observed' data is a simulated dataset
#' with \eqn{Pm = 0.35} and \eqn{Pp = 0.001}. The lattice has 27 \code{rows} and 36 \code{cols}
#' and there are \code{num_obs = 144} observations after time 0
#' (to mimic images being taken every 5 minutes for 12 hours).
#' The simulation is based on there initially being 110 cells in the assay.
#' 
#' Further information about the specific choices of tuning parameters
#' used in BSL and BSLasso can be found in An et al. (2018).
#' 
#' \itemize{
#'  \item \code{data}:  The \code{rows} \eqn{x} \code{cols} \eqn{x} \code{num_obs} array of the cell presences at times 1:144.
#'  \item \code{sim_options}: Values of \code{sim_options} relevant to this example.
#'  \item \code{sum_options}: Values of \code{sim_options} relevant to this example, i.e. just the value of \code{Yinit}.
#'  \item \code{start}: A vector of suitable initial values of the parameters for MCMC.
#'  \item \code{cov}: Covariance matrix of the multivariate normal random walk, in the form of a \eqn{2 x 2} matrix.
#' }
#'
#' @examples
#' \dontrun{
#' require(doParallel) # You can use a different package to set up the parallel backend
#' 
#' # Loading the data for this example
#' data(cell)
#' true_cell <- c(0.35, 0.001)
#' 
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores())
#' registerDoParallel(cl)
#' 
#' # Performing BSL
#' resultCellBSL <- bsl(cell$data, n = 5000, M = 10000, start = cell$start, cov_rw = cell$cov,
#'                  fn_sim = cell_sim, fn_sum = cell_sum, fn_prior = cell_prior,
#'                  sim_options = cell$sim_options, sum_options = cell$sum_options,
#'                  parallel = TRUE, parallel_packages = 'BSL',
#'                  theta_names = c('Pm', 'Pp'))
#' summary(resultCellBSL)
#' plot(resultCellBSL, true_value = true_cell, thin = 20)
#' 
#' # Performing tuning for BSLasso
#' lambda_all <- list(exp(seq(0.5,2.5,length.out=20)), exp(seq(0,2,length.out=20)), 
#'                    exp(seq(-1,1,length.out=20)), exp(seq(-1,1,length.out=20)))
#'
#' sp_cell <- selectPenalty(ssy = cell_sum(cell$data, cell$sum_options),
#'                  n = c(500, 1000, 1500, 2000), lambda_all, theta = true_cell,
#'                  M = 100, sigma = 1.5, fn_sim = cell_sim,
#'                  fn_sum = cell_sum, sim_options = cell$sim_options,
#'                  sum_options = cell$sum_options, parallel_sim = TRUE,
#'                  parallel_sim_packages = 'BSL', parallel_main = TRUE)
#' 
#' sp_cell
#' plot(sp_cell)
#' 
#' # Performing BSLasso with a fixed penalty
#' resultCellBSLasso <- bsl(cell$data, n = 1500, M = 10000, start = cell$start,
#'                  cov_rw = cell$cov, fn_sim = cell_sim, fn_sum = cell_sum,
#'                  penalty = 1.35, fn_prior = cell_prior,
#'                  sim_options = cell$sim_options, sum_options = cell$sum_options,
#'                  parallel = TRUE, parallel_packages = 'BSL',
#'                  theta_names = c('Pm', 'Pp'))
#' summary(resultCellBSLasso)
#' plot(resultCellBSLasso, true_value = true_cell, thin = 20)
#' 
#' # Plotting the results together for comparison
#' combinePlotsBSL(resultCellBSL, resultCellBSLasso, true_value = true_cell, thin = 20)
#' 
#' # Closing the parallel pools
#' stopCluster(cl)
#' }
#'
#' @references
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018). Accelerating
#' Bayesian synthetic likelihood with the graphical lasso. \url{https://eprints.qut.edu.au/102263/}
#'
#' Johnston, S., Simpson, M. J., McElwain, D. L. S., Binder, B. J. &
#' Ross, J. V. (2014). Interpreting Scratch Assays Using Pair Density
#' Dynamic and Approximate Bayesian Computation. Open Biology, 4, 1-11.
#'
#' Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018).
#' Bayesian synthetic likelihood. To appear in Journal of Computational and Graphical Statistics. \url{https://eprints.qut.edu.au/92795/}
#' 
#' @author 								Ziwen An, Christopher C. Drovandi and Leah F. South
#' @name cell
NULL



#' The function \code{cell_sim(theta, sim_options)} simulates data from the model, using C++ in the backend.
#' @rdname cell
cell_sim <-function(theta, sim_options) {
    Pm <- theta[1]
    Pp <- theta[2]
    Y <- simulate_cell(sim_options$Yinit, sim_options$rows, sim_options$cols, Pm, Pp, sim_options$sim_iters, sim_options$num_obs)
    return(Y)
}

#' The function \code{cell_sum(Y,sum_options)} calculates the summary statistics for this example.
#' @rdname cell
cell_sum <- function(Y, sum_options) {
    num_obs = dim(Y)[3]
    summ_stat = numeric(num_obs+1)
  
    # Hamming distances between cell locations across time
    summ_stat[1] = sum(abs(sum_options$Yinit-Y[, , 1]))
    for (i in 2:num_obs) {
        summ_stat[i] = sum(abs(Y[, , i-1]-Y[, , i]))
    }

    # Total number of cells in the final time period
    summ_stat[num_obs + 1] = sum(Y[, , num_obs])
	
	return(summ_stat)
}

#' The function \code{cell_prior(theta)} evaluates the prior at the chosen parameters. 
#' @rdname cell
cell_prior <- function(theta) {
    theta[1] > 0 & theta[1] < 1 & theta[2] > 0 & theta[2] < 1
}
