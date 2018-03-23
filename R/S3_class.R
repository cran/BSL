# print a bsl class result
print.bsl <- function(x, digits = max(3L, getOption("digits") - 3L),...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (nrow(x$theta)) {
        cat("Summary of theta:\n")
        summ <- summary(x$theta)
        attr(summ, 'dimnames') = list(NULL, x$theta_names)
        print.default(format(summ, digits = digits), print.gap = 2L, 
                      quote = FALSE)
    }
    else cat("No theta\n")
    if (length(x$loglike)) {
        cat("Summary of loglikelihood:\n")
        summ <- summary(x$loglike)
        print.default(format(summ, digits = digits), print.gap = 2L, 
                      quote = FALSE)
    }
    else cat("No loglikelihood\n")
    if (length(x$acceptanceRate)) {
        cat("Acceptance Rate:\n")
        print.default(format(x$acceptanceRate, digits = digits), print.gap = 2L, 
                      quote = FALSE)
    }
    else cat("No acceptance rate\n")
    if (length(x$earlyRejectionRate)) {
        cat("Early Rejection Rate:\n")
        print.default(format(x$earlyRejectionRate, digits = digits), print.gap = 2L, 
                      quote = FALSE)
    }
    else cat("No early rejection rate\n")
    cat("\n")
}

# print a penbsl class result
print.penbsl <- function(x, digits = max(3L, getOption("digits") - 4L), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(x$resultsDF)) {
        r1 <- x$resultsDF[which(x$resultsDF$sigmaOpt), c('n', 'penalty', 'sigma')]
        cat('Penalty selected based on the standard deviation of the loglikelihood:\n')
        print(format(r1, digits = digits))
    } else {
        cat("No result to show\n")
    }
    return(invisible(r1))
}

# summarise a bsl class result
summary.bsl <- function(..., n = NULL, digits = c(3, 0, 0), scale = 1000000) {
    if (length(list(...)) == 0) {
        stop('no expressions to be evaluated')
    }
    arguments <- list(...)
    na <- length(arguments)
    if (!is.null(n) & length(n) != na) {
        stop('length of n and arguments dismatch')
    }
    parameters <- names(match.call()[2 : (1+na)])
    if (is.null(parameters)) {
        parameters <- as.character(match.call()[2:(1+na)])
    }
    p <- sapply(lapply(arguments, '[[', 'theta'), FUN = ncol)
    if (length(unique(p)) != 1L) {
        stop(parameters, 'have different length of parameters')
    } else {
        p <- unique(p)
    }
    summ <- data.frame(matrix(nrow = na, ncol = 3L + p))
    dimnames(summ) <- list(parameters, c('n', 'penalty', 'acc. rate (%)', paste('scaled ESS', arguments[[1]]$theta_names)))
    for (i in 1 : na) {
        if (is.null(n)) {
            summ$n[i] <- as.numeric(as.character(arguments[[i]]$call['n']))
        } else {
            summ$n[i] <- n[i]
        }
        M <- nrow(arguments[[i]]$theta)
        suppressWarnings(summ$penalty[i] <- round(as.numeric(as.character(arguments[[i]]$call['penalty'])), digits[1]))
        summ[i, 3] <- round(arguments[[i]]$acceptanceRate * 100,  digits[2])
        suppressWarnings(summ[i, 4 : (3+p)] <- round(effectiveSize(arguments[[i]]$theta) / summ$n[i] / M * scale, digits[3]))
    }
    summ
}



#' Plotting BSL and BSLasso Results
#'
#' @description There are several functions included in this package to help visualise the results of the MCMC run and of the tuning process.
#'
#' @param x                         An object to be plotted.
#' @param ...						The first inputs to the \code{combinePlotsBSL} function should be two or more results of type \code{bsl}, separated by commas.
#' @param true_value                A set of values to be included on the plots as a reference line. The default is \code{NULL}.
#' @param thin                      The gap between samples to be taken when thinning the MCMC draws. The default is 1 (no thinning).
#' @param logscale                  An indicator for whether the penalty values should be shown on the log scale in plots of \code{penbsl} objects. The default is \code{TRUE}.
#'
#' @author 								Ziwen An, Christopher C. Drovandi and Leah F. South
#' 
#' @name plot
NULL


#' The function \code{combinePlotsBSL} can be used to plot multiple BSL and BSLasso densities together, optionally with the true values for the parameters.
#' @rdname plot
combinePlotsBSL <- function(..., true_value = NULL, thin = 1) {
    if (length(list(...)) == 0) {
        stop('no expressions to be evaluated')
    }
    arguments <- list(...)
    na <- length(arguments)
    parameters <- names(match.call()[2:(1+na)])
    if (is.null(parameters)) {
        parameters <- as.character(match.call()[2:(1+na)])
    }
    thetaList <- lapply(arguments, '[[', 'theta')
    n <- sapply(thetaList, FUN = nrow)
    p <- ncol(thetaList[[1]])
    a <- floor(sqrt(p))
    b <- ceiling(p / a)
    theta_names <- colnames(thetaList[[1]])
    for (i in 2 : na) {
        colnames(thetaList[[i]]) <- theta_names
    }
    x <- list()
    for (i in 1 : na) {
        x[[i]] <- data.frame(thetaList[[i]][seq(1, n[i], by = thin), ], label = parameters[i])
    }
    x <- do.call('rbind', x)
    plist <- array(list(), p)
    for (i in 1 : p) {
        plist[[i]] <- ggplot(x, aes_string(x = theta_names[i])) + 
            geom_density(aes_string(color = 'label', linetype = 'label'), size = 1) + # using aes_string to satisfy R CMD check
            geom_hline(yintercept = 0, color = 'black', alpha = 0.6, size = 0.7) + {
                if (!is.null(true_value)) {
                    geom_vline(xintercept = true_value[i], color = 'forestgreen', linetype = 'dashed', size = 0.5)
                }
            } + {
                labs(y = NULL)
            } + {
			    theme(axis.title=element_text(size = 12), legend.position = 'none', legend.text = element_text(size = 8, angle = 0), legend.direction = 'horizontal', 
				        legend.title = element_blank())
            }
    }
    g <- ggplotGrob(plist[[1]] + {
	    theme(legend.position = 'bottom')
	})$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lwidth <- sum(legend$width)
	lheight <- sum(legend$height)
    combined <- arrangeGrob(do.call(arrangeGrob, c(plist, nrow = a, ncol = b, top = 'Posteriors')), legend, nrow = 2, heights = unit.c(unit(1, "npc") - lheight, lheight))
	grid.newpage()
    grid.draw(combined)

}


#' The function \code{plot.BSL} can be used to plot BSL and BSLasso densities, optionally with the true values for the parameters.
#' @rdname plot
plot.bsl <- function(x, true_value = NULL, thin = 1, ...) {
    n <- nrow(x$theta)
    p <- ncol(x$theta)
    a <- floor(sqrt(p))
    b <- ceiling(p / a)
    
    if (!is.null(true_value) & length(true_value) != p) {
        stop('length of true_value dismatches number of parameters')
    }
    
    objData <- data.frame(x$theta[seq(1, n, by = thin), ])
    theta_names <- colnames(objData)
    plist <- list()
    for (i in 1 : p) {
        plist[[i]] <- ggplot(objData, aes_string(x = theta_names[i])) + 
            geom_density(color = 'darkblue', linetype = 'solid', size = 1) + 
            geom_hline(yintercept = 0, colour = "grey90", size = 0.75) + {
                if (!is.null(true_value)) {
                    geom_vline(xintercept = true_value[i], color = 'forestgreen', linetype = 'dashed', size = 0.5)
                }
            } + {
                labs(y = NULL)
            }
            
    }
    do.call('grid.arrange', c(plist, nrow = a, ncol = b, top = 'Posteriors'))
}

#' The function \code{plot.penbsl} can be used to plot the results from tuning to select the optimal penalty for BSLasso.
#' @rdname plot
plot.penbsl <- function(x, logscale = TRUE, ...) {
    sigma <- sigmaOpt <- penalty <- logPenalty <- NULL # to satisfy R CMD check
    x <- x$resultsDF
	a <- floor(sqrt(length(unique(x$n))))
	b <- ceiling(length(unique(x$n)) / a)
    nRepeats <- sapply(unique(x$n), FUN = function(xx) sum(x['n'] == xx))
    yPosSigma <- sapply(unique(x$n), FUN = function(xx) mean(range(x[which(x$n == xx), 'sigma'])))
    textYSigma <- c(unlist(mapply(yPosSigma, nRepeats, FUN = rep)))
    sigmaOpt.NA <- x$sigmaOpt
    sigmaOpt.NA[which(sigmaOpt.NA == FALSE)] <- NA
    sigmaOpt.NA[which(sigmaOpt.NA == TRUE)] <- 1
    logPenalty <- log(x$penalty)
    if (logscale) {
		ggplot(data = x, aes(x = logPenalty, y = sigma)) + 
		    geom_line(color = 'darkblue', linetype = 'dashed', size = 1) + 
		    facet_wrap( ~ n, scales = 'free', nrow = a, ncol = b, labeller = label_both) + 
		    geom_vline(aes(xintercept = logPenalty*sigmaOpt.NA), na.rm = TRUE, color = 'forestgreen', linetype = 4) + 
		    geom_label(aes(x = logPenalty*sigmaOpt.NA, y = textYSigma, label = paste('lambda == ', round(penalty, 3))), 
		               hjust = 0.5, vjust = "inward", parse = TRUE, color = 'white', fill = '#FE66A9', size = 2.7, 
					   alpha = 0.8, na.rm = TRUE) + 
		    labs(x = 'log(penalty)', title = 'Penalty Selection') + 
			theme(plot.title = element_text(size = 14, hjust = 0.5)) + 
			theme(strip.text.x = element_text(size = 12, face = 'bold'), axis.title = element_text(size = 12))
	} else {
		ggplot(data = x, aes(x = penalty, y = sigma)) + 
		    geom_line(color = 'darkblue', linetype = 'dashed', size = 1) + 
		    facet_wrap( ~ n, scales = 'free', nrow = a, ncol = b, labeller = label_both) + 
		    geom_vline(aes(xintercept = penalty*sigmaOpt.NA), na.rm = TRUE, color = 'forestgreen', linetype = 4) + 
		    geom_label(aes(x = penalty*sigmaOpt.NA, y = textYSigma, label = paste('lambda == ', round(penalty, 3))), 
		               hjust = 0.5, vjust = "inward", parse = TRUE, color = 'white', fill = '#FE66A9', size = 2.7, 
					   alpha = 0.8, na.rm = TRUE) + 
		    labs(x = 'penalty', title = 'Penalty Selection') + 
			theme(plot.title = element_text(size = 14, hjust = 0.5)) + 
			theme(strip.text.x = element_text(size = 12, face = 'bold'), axis.title = element_text(size = 12))
    }
}
