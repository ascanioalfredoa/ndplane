#' Modified Maxent Model Builder Using glmnet
#'
#' This function is a modified version of the original maxnet function from the
#' \pkg{maxnet} package, adapted for specific modeling needs.
#' For details on the original function, see \code{\link[maxnet]{maxnet}}.
#'
#' @inheritParams maxnet::maxnet
#' @param ... Additional parameters passed to \code{glmnet}.
#'
#' @return An object of class \code{maxnet} containing the fitted model and metadata.
#'
#' @details
#' This version modifies the fitting process by adding unique presence samples to the background,
#' setting controlled regularization adjustments, and calculating model entropy and penalty details.
#'
#' @importFrom glmnet glmnet glmnet.control
#' @export
maxnet2 <- function (p, data, f = maxnet.formula(p, data), regmult = 1, 
    regfun = maxnet.default.regularization, addsamplestobackground = TRUE, 
    ...) 
{
    if (anyNA(data)) 
        stop("NA values in data table. Please remove them and rerun.")
    if (addsamplestobackground) {
        pdata <- data[p == 1, , drop = FALSE]
        ndata <- data[p == 0, , drop = FALSE]
        toadd <- apply(pdata, 1, function(rr) !any(apply(ndata, 
            1, function(r) identical(r, rr))))
        p <- c(p, rep(0, sum(toadd)))
        data <- rbind(data, pdata[toadd, , drop = FALSE])
    }
    mm <- model.matrix(f, data)
    reg <- regfun(p, mm) * regmult
    weights <- p + (1 - p) * 100
    glmnet::glmnet.control(pmin = 1e-08, fdev = 0)
    model <- glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial", 
        standardize = FALSE, penalty.factor = reg, lambda = 10^(seq(4, 
            0, length.out = 200)) * sum(reg)/length(reg) * sum(p)/sum(weights), 
        weights = weights, ...)
    class(model) <- c("maxnet", class(model))
    if (length(model$lambda) < 200) {
        msg <- "Error: glmnet failed to complete regularization path.  Model may be infeasible."
        if (!addsamplestobackground) 
            msg <- paste(msg, " Try re-running with addsamplestobackground=T.")
        stop(msg)
    }
    bb <- model$beta[, 200]
    model$betas <- bb[bb != 0]
    model$alpha <- 0
    rr <- predict.maxnet(model, data[p == 0, , drop = FALSE], 
        type = "exponent", clamp = FALSE)
    raw <- rr/sum(rr)
    model$entropy <- -sum(raw * log(raw))
    model$alpha <- -log(sum(rr))
    model$penalty.factor <- reg
    model$featuremins <- apply(mm, 2, min)
    model$featuremaxs <- apply(mm, 2, max)
    vv <- (sapply(data, class) != "factor")
    model$varmin <- apply(data[, vv, drop = FALSE], 2, min)
    model$varmax <- apply(data[, vv, drop = FALSE], 2, max)
    means <- apply(data[p == 1, vv, drop = FALSE], 2, mean)
    majorities <- sapply(names(data)[!vv], function(n) which.max(table(data[p == 
        1, n, drop = FALSE])))
    names(majorities) <- names(data)[!vv]
    model$samplemeans <- unlist(c(means, majorities))
    model$levels <- lapply(data, levels)
    model
}
