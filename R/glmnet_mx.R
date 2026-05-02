#' Predict method for glmnet_mx objects
#'
#' @param object A glmnet_mx model object
#' @param newdata Data frame of new environmental data
#' @param type Type of prediction: "link", "response", "exponent", "cloglog"
#' @param clamp Logical, whether to clamp predictors to the range used during training
#' @param ... Additional arguments
#' @export
predict.glmnet_mx <- function(object, newdata, type = c("link", "response", "exponent", "cloglog"), clamp = TRUE, ...) {
  type <- match.arg(type)

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

  # Use NextMethod or direct call to avoid infinite recursion
  class(object) <- setdiff(class(object), "glmnet_mx")
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
                      regfun = maxnet::maxnet.default.regularization,
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

  # Store sample means for response curves
  numeric_vars <- sapply(data, is.numeric)
  model$samplemeans <- colMeans(data[p == 1, numeric_vars, drop = FALSE])

  rr <- predict.glmnet_mx(model, data[p == 0, , drop = FALSE], type = "exponent", clamp = FALSE)
  rr <- rr + .Machine$double.eps
  model$entropy <- -sum((rr/sum(rr)) * log(rr/sum(rr)))
  model$alpha <- -log(sum(rr))

  return(model)
}
