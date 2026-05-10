#' Default Regularization for Maxent-like Models
#'
#' @param p Presence/background vector
#' @param m Model matrix
#' @export
default_regularization <- function(p, m) {
    isproduct <- function(x) grepl(":", x) & !grepl("\\(", x)
    isquadratic <- function(x) grepl("^I\\(.*\\^2\\)", x)
    ishinge <- function(x) grepl("^hinge\\(", x)
    isthreshold <- function(x) grepl("^thresholds\\(", x)
    iscategorical <- function(x) grepl("^categorical\\(", x)
    regtable <- function(name, default) {
        if (ishinge(name))
            return(list(c(0, 1), c(0.5, 0.5)))
        if (iscategorical(name))
            return(list(c(0, 10, 17), c(0.65, 0.5, 0.25)))
        if (isthreshold(name))
            return(list(c(0, 100), c(2, 1)))
        default
    }
    lregtable <- list(c(0, 10, 30, 100), c(1, 1, 0.2, 0.05))
    qregtable <- list(c(0, 10, 17, 30, 100), c(1.3, 0.8, 0.5,
        0.25, 0.05))
    pregtable <- list(c(0, 10, 17, 30, 100), c(2.6, 1.6, 0.9,
        0.55, 0.05))
    mm <- m[p == 1, , drop = FALSE]
    np <- nrow(mm)
    lqpreg <- lregtable
    if (sum(isquadratic(colnames(mm))))
        lqpreg <- qregtable
    if (sum(isproduct(colnames(mm))))
        lqpreg <- pregtable
    classregularization <- sapply(colnames(mm), function(n) {
        t <- regtable(n, lqpreg)
        stats::approx(t[[1]], t[[2]], np, rule = 2)$y
    })/sqrt(np)
    ishinge_vec <- grepl("^hinge\\(", colnames(mm))
    hmindev <- sapply(1:ncol(mm), function(i) {
        if (!ishinge_vec[i])
            return(0)
        avg <- mean(mm[, i])
        std <- max(stats::sd(mm[, i]), 1/sqrt(np))
        std * 0.5/sqrt(np)
    })
    tmindev <- sapply(1:ncol(mm), function(i) {
        ifelse(isthreshold(colnames(mm)[i]) && (sum(mm[, i]) ==
            0 || sum(mm[, i]) == nrow(mm)), 1, 0)
    })
    pmax(0.001 * (apply(m, 2, max) - apply(m, 2, min)), hmindev,
        tmindev, apply(as.matrix(mm), 2, stats::sd) * classregularization)
}

#' Predict method for glmnet_mx objects
#'
#' @param object A glmnet_mx model object
#' @param newdata Data frame of new environmental data
#' @param type Type of prediction: "link", "response", "exponent", "cloglog"
#' @param clamp Logical, whether to clamp predictors to the range used during training
#' @param ... Additional arguments
#' @export
predict.glmnet_mx <- function(object, newdata, type = c("link", "response", "exponent", "cloglog"), clamp = TRUE, ...) {
  # If type is not one of our custom types, delegate to next method (glmnet)
  if (!missing(type) && !type %in% c("link", "response", "exponent", "cloglog")) {
    return(NextMethod("predict"))
  }

  type <- match.arg(type)

  if (missing(newdata) || is.null(newdata)) {
    return(NextMethod("predict"))
  }

  if (clamp) {
    for (v in names(object$varmin)) {
      if (v %in% names(newdata) && !is.na(object$varmin[v])) {
        newdata[[v]][newdata[[v]] < object$varmin[v]] <- object$varmin[v]
        newdata[[v]][newdata[[v]] > object$varmax[v]] <- object$varmax[v]
      }
    }
  }

  f <- object$formula
  mm <- model.matrix(f, newdata)

  res <- stats::predict(object, newx = mm, s = object$lambda[200], type = "link")[, 1]
  res <- res + object$alpha

  if (type == "exponent") return(exp(res))
  if (type == "response") return(exp(res) / (1 + exp(res)))
  if (type == "cloglog") return(1 - exp(-exp(res)))
  return(res)
}

#' Maxent-like glmnet models
#'
#' @param p A vector of binary presence-background labels (1 presence, 0 background)
#' @param data A data.frame containing predictor variables
#' @param f A formula for the model
#' @param regmult Regularization multiplier
#' @param regfun Function to calculate regularization penalties
#' @param addsamplestobackground Logical, whether to add presence points to background
#' @param weights Numeric vector of weights
#' @param ... Additional arguments to glmnet
#' @importFrom glmnet glmnet glmnet.control
#' @export
glmnet_mx <- function(p,
                      data,
                      f,
                      regmult = 1.0,
                      regfun = default_regularization,
                      addsamplestobackground = TRUE,
                      weights = NULL,
                      ...) {

  if (anyNA(data)) stop("NA values in data table.")

  iniweight <- is.null(weights)
  if (iniweight) weights <- ifelse(p == 1, 1, 100)

  if (addsamplestobackground) {
    pdata <- data[p == 1, , drop = FALSE]
    ndata <- data[p == 0, , drop = FALSE]

    wadd <- !do.call(paste, pdata) %in% do.call(paste, ndata)
    if (sum(wadd) > 0) {
      p <- c(p, rep(0, sum(wadd)))
      data <- rbind(data, pdata[wadd, , drop = FALSE])
      if (iniweight) {
        weights <- c(weights, rep(100, sum(wadd)))
      } else {
        pweight <- weights[p == 1]
        weights <- c(weights, pweight[wadd])
      }
    }
  }

  mm <- model.matrix(f, data)

  reg <- regfun(p, mm) * regmult
  lambdas <- 10^(seq(4, 0, length.out = 200)) * sum(reg) / length(reg) * sum(p) / sum(weights)

  glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  model <- glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",
                          standardize = FALSE, penalty.factor = reg,
                          lambda = lambdas, weights = weights, ...)

  class(model) <- c("glmnet_mx", class(model))

  bb <- stats::coef(model, s = model$lambda[200])

  model$betas <- bb[-1, 1]
  model$betas <- model$betas[model$betas != 0]
  model$alpha <- 0
  model$formula <- f

  model$varmin <- apply(data, 2, function(x) if(is.numeric(x)) min(x) else NA)
  model$varmax <- apply(data, 2, function(x) if(is.numeric(x)) max(x) else NA)

  numeric_vars <- sapply(data, is.numeric)
  model$samplemeans <- colMeans(data[p == 1, numeric_vars, drop = FALSE])

  rr <- predict.glmnet_mx(model, data[p == 0, , drop = FALSE], type = "exponent", clamp = FALSE)
  rr <- rr + .Machine$double.eps
  model$entropy <- -sum((rr/sum(rr)) * log(rr/sum(rr)))
  model$alpha <- -log(sum(rr))

  return(model)
}
