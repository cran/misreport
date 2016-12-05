#' Add two numbers
#' 
#' The function \code{logAdd} adds together two numbers using their log
#' to prevent under-/over-flow
#' 
#' @param x log of the first number.
#' @param y log of the second number.
#' @return Log of the sum of \code{exp(x)} and \code{exp(y)}.
#'
#' @keywords internal
logAdd <- function(lx, ly) {
  max(lx, ly) + log1p(exp(-abs(lx - ly)))
}


#' Misreport sub-model m-step
#'
#' The maximization step in the EM algorithm called by \code{\link{listExperiment}}
#' for the misreport sub-model.
#' 
#' @keywords internal
#'
#' @importFrom stats .getXlevels as.formula binomial coef
#'             cov dbinom model.matrix model.frame model.response
#'             na.pass plogis pt rnorm runif glm 
mstepMisreport <- function(y, x.misreport, w, treat,
                              misreport.treatment, weight) {

  lrep <- rep(c(1, 0), each = length(y)) # Misreport is the first column of w

  if(misreport.treatment == TRUE) {
    xrep <- as.matrix(rbind(cbind(x.misreport, treat), cbind(x.misreport, treat)))
  } else if(misreport.treatment == FALSE) {
    xrep <- as.matrix(rbind(x.misreport, x.misreport))
  }

  wrep <- c(w[, 1] + log(weight), w[, 2] + log(weight)) # Misreport is the first column of w

  lrep <- lrep[wrep > -Inf]
  xrep <- xrep[wrep > -Inf, , drop = FALSE]
  wrep <- wrep[wrep > -Inf]

  X <- xrep

  if(ncol(X) == 1) {
    fit.misreport <- glm(cbind(lrep, 1 - lrep) ~ 1,     weights = exp(wrep), family = binomial)
  } else if(ncol(X) > 1) {
    fit.misreport <- glm(cbind(lrep, 1 - lrep) ~ -1 + X, weights = exp(wrep), family = binomial)
  }

  coefs <- coef(fit.misreport)
  names(coefs) <- gsub("^X1|^X2|^X3|^X", "", names(coefs))

  return(coefs)

}


#' Sensitive-item sub-model m-step
#'
#' The maximization step in the EM algorithm called by \code{\link{listExperiment}}
#' for the sensitive-item sub-model.
#' 
#' @keywords internal
#'
mstepSensitive <- function(y, treat, x.sensitive, w, d, sensitive.response,
                           weight, model.misreport) {

  if(model.misreport == TRUE) {
    zrep <- rep(c(sensitive.response, abs(1 - sensitive.response)), each = length(y))
    xrep <- as.matrix(rbind(x.sensitive, x.sensitive))
    wrep <- c(apply(w[, 1:2], 1, function(x) logAdd(x[1], x[2])) + log(weight), w[, 3] + log(weight))

    # wrep <- c(apply(w[, 1:2], 1, sum) * weight, w[, 3] * weight)

    zrep <- zrep[wrep > -Inf]
    xrep <- xrep[wrep > -Inf, , drop = FALSE]
    wrep <- wrep[wrep > -Inf]

    X <- xrep

    if(ncol(X) == 1) fit <- glm(cbind(zrep, 1 - zrep) ~ 1,     weights = exp(wrep), family = binomial)
    if(ncol(X) > 1)  fit <- glm(cbind(zrep, 1 - zrep) ~ -1 + X, weights = exp(wrep), family = binomial)

    coefs <- coef(fit)
  } else {
    zrep <- rep(c(1, 0), each = length(y))
    xrep <- as.matrix(rbind(x.sensitive, x.sensitive))
    wrep <- c(w[, 1] + log(weight), w[, 2] + log(weight))

    zrep <- zrep[wrep > -Inf]
    xrep <- xrep[wrep > -Inf, , drop = FALSE]
    wrep <- wrep[wrep > -Inf]

    X <- xrep

    if(ncol(X) == 1) fit <- glm(cbind(zrep, 1 - zrep) ~ 1,     weights = exp(wrep), family = binomial)
    if(ncol(X) > 1)  fit <- glm(cbind(zrep, 1 - zrep) ~ -1 + X, weights = exp(wrep), family = binomial)

    coefs <- coef(fit)
  }

  names(coefs) <- gsub("^X1|^X2|^X3|^X", "", names(coefs))

  return(coefs)

}


#' Control-items sub-model m-step
#'
#' The maximization step in the EM algorithm called by \code{\link{listExperiment}}
#' for the control-items sub-model.
#' 
#' @keywords internal
#'
mstepControl <- function(y, treat, J, x.control, w, d, sensitive.response,
                         weight, model.misreport, control.constraint) {

  if(model.misreport == TRUE) {

    if(control.constraint == "none") {
      yrep <- c((y - treat * as.numeric(sensitive.response == 1)),
                (y - treat * as.numeric(sensitive.response == 1)),
                (y - treat * as.numeric(sensitive.response == 0)))
      xrep <- as.matrix(rbind(x.control, x.control, x.control))
      zrep1 <- rep(c(1, 0, 0), each = length(y)) # Misreport sensitive
      zrep2 <- rep(c(sensitive.response,
                     sensitive.response,
                     1 - sensitive.response), each = length(y)) # Truthful sensitive
      wrep <- c(w[, 1] + log(weight), w[, 2] + log(weight), w[, 3] + log(weight))

      yrep <- yrep[wrep > -Inf]
      xrep <- xrep[wrep > -Inf, , drop = FALSE]
      zrep1 <- zrep1[wrep > -Inf]
      zrep2 <- zrep2[wrep > -Inf]
      wrep <- wrep[wrep > -Inf]

      X <- cbind(xrep, U = zrep1, Z = zrep2)
      control.fit <- glm(cbind(yrep, J - yrep) ~ -1 + X, weights = exp(wrep), family = binomial)
      coefs <- coef(control.fit)

    } else if(control.constraint == "partial") { # U* = 0
      yrep <- c((y - treat * as.numeric(sensitive.response == 1)),
                (y - treat * as.numeric(sensitive.response == 0)))
      xrep <- as.matrix(rbind(x.control, x.control))
      zrep1 <- rep(c(sensitive.response, 1 - sensitive.response), each = length(y)) # Sensitive
      wrep <- c(apply(w[, 1:2], 1, function(x) logAdd(x[1], x[2])) + log(weight), w[, 3] + log(weight))

      yrep <- yrep[wrep > -Inf]
      xrep <- xrep[wrep > -Inf, , drop = FALSE]
      zrep1 <- zrep1[wrep > -Inf]
      wrep <- wrep[wrep > -Inf]

      X <- cbind(xrep, Z = zrep1)
      control.fit <- glm(cbind(yrep, J - yrep) ~ -1 + X, weights = exp(wrep), family = binomial)
      coefs <- coef(control.fit)

    } else if(control.constraint == "full") { # U* = Z* = 0
      yrep <- c((y - treat * as.numeric(sensitive.response == 1)),
                (y - treat * as.numeric(sensitive.response == 0)))
      xrep <- as.matrix(rbind(x.control, x.control))
      wrep <- c(apply(w[, 1:2], 1, function(x) logAdd(x[1], x[2])) + log(weight), w[, 3] + log(weight))

      yrep <- yrep[wrep > -Inf]
      xrep <- xrep[wrep > -Inf, , drop = FALSE]
      wrep <- wrep[wrep > -Inf]

      X <- xrep
      
      if(ncol(X) == 1) control.fit <- glm(cbind(yrep, J - yrep) ~ 1    , weights = exp(wrep), family = binomial)
      if(ncol(X) > 1)  control.fit <- glm(cbind(yrep, J - yrep) ~ -1 + X, weights = exp(wrep), family = binomial)
      coefs <- coef(control.fit)
    }
  } else {

      yrep <- c(y - treat, y)
      xrep <- as.matrix(rbind(x.control, x.control))
      zrep1 <- rep(c(1, 0), each = length(y))
      zrep2 <- rep(c(0, 1), each = length(y))
      wrep <- c(w[, 1] + log(weight), w[, 2] + log(weight))

      yrep <- yrep[wrep > -Inf]
      xrep <- xrep[wrep > -Inf, , drop = FALSE]
      zrep1 <- zrep1[wrep > -Inf]
      zrep2 <- zrep2[wrep > -Inf]
      wrep <- wrep[wrep > -Inf]

    if(control.constraint == "none") {
      X <- cbind(xrep, Z = zrep1)
      fit.partial <- glm(cbind(yrep, J - yrep) ~ -1 + X, weights = exp(wrep), family = binomial)
      coefs <- c(coef(fit.partial))
    }

    if(control.constraint == "full") {
      X <- xrep
      if(ncol(X) == 1) fit.full <- glm(cbind(yrep, J - yrep) ~ 1    , weights = exp(wrep), family = binomial)
      if(ncol(X) > 1)  fit.full <- glm(cbind(yrep, J - yrep) ~ -1 + X, weights = exp(wrep), family = binomial)
      coefs <- c(coef(fit.full))
    }
  }

  names(coefs) <- gsub("^X1|^X2|^X3|^X", "", names(coefs))
  names(coefs)[names(coefs) == "(Intercept):1"] <- "(Intercept)"

  return(coefs)

}


#' Outcome sub-model m-step
#'
#' The maximization step in the EM algorithm called by \code{\link{listExperiment}}
#' for the outcome sub-model.
#' 
#' @keywords internal
#'
mstepOutcome <- function(y, treat, x.outcome, w, d, sensitive.response, o,
                         trials, weight, outcome.model, model.misreport,
                         outcome.constrained, control.constraint) {

  coefs.aux <- NULL

  if(outcome.constrained == TRUE) {
    if(model.misreport == TRUE) {
      xrep <- as.matrix(rbind(x.outcome, x.outcome))
      zrep <- rep(c(1, 0), each = length(y))
      orep <- as.matrix(c(o, o))
      trialsrep <- as.matrix(c(trials, trials))
      wrep <- c(apply(w[, 1:2], 1, function(x) logAdd(x[1], x[2])) + log(weight), w[, 3] + log(weight))
    } else {
      xrep <- as.matrix(rbind(x.outcome, x.outcome))
      zrep <- rep(c(1, 0), each = length(y))
      orep <- as.matrix(c(o, o))
      trialsrep <- as.matrix(c(trials, trials))
      wrep <- c(w[, 1] + log(weight), w[, 2] + log(weight))
    }

    xrep <- xrep[wrep > -Inf, , drop = FALSE]
    zrep <- zrep[wrep > -Inf]
    orep <- orep[wrep > -Inf]
    trialsrep <- trialsrep[wrep > -Inf]
    wrep <- wrep[wrep > -Inf]

    X <- cbind(xrep, Z = zrep)[, -1, drop = FALSE]

    if(outcome.model == "logistic") {
      fit.constrained.logistic <- glm(cbind(orep, 1 - orep) ~ 1 + X, weights = exp(wrep), family = binomial)
      coefs <- coef(fit.constrained.logistic)
    } else if(outcome.model == "binomial") {
      fit.constrained.binomial <- glm(cbind(orep, trialsrep - orep) ~ 1 + X, family = binomial, weights = exp(wrep))
      coefs <- coef(fit.constrained.binomial)
    } else if(outcome.model == "betabinomial") {
      fit.constrained.betabinomial <- VGAM::vglm(cbind(orep, trialsrep - orep) ~ 1 + X, VGAM::betabinomial, weights = exp(wrep))
      coefs <- coef(fit.constrained.betabinomial)[-2]
      coefs.aux <- c(rho = mean(fit.constrained.betabinomial@misc$rho))
    }
  } else if(outcome.constrained == FALSE) {
    if(model.misreport == TRUE) {
      xrep <- as.matrix(rbind(x.outcome, x.outcome, x.outcome))
      zrep1 <- rep(c(1, 0, 0), each = length(y))
      zrep2 <- rep(c(1, 1, 0), each = length(y))
      orep <- as.matrix(c(o, o, o))
      trialsrep <- as.matrix(c(trials, trials, trials))
      wrep <- c(w[, 1] + log(weight), w[, 2] + log(weight), w[, 3] + log(weight))
    } else {
      stop("\noutcome.constrained = TRUE is only possible when a direct question is included.")
    }

    xrep <- xrep[wrep > -Inf, , drop = FALSE]
    zrep1 <- zrep1[wrep > -Inf]
    zrep2 <- zrep2[wrep > -Inf]
    orep <- orep[wrep > -Inf]
    trialsrep <- trials[wrep > -Inf]
    wrep <- wrep[wrep > -Inf]

    X <- cbind(xrep, U = zrep1, Z = zrep2)[, -1, drop = FALSE]

    if(outcome.model == "logistic") {
      fit.unconstrained.logitistic <- glm(cbind(orep, 1 - orep) ~ 1 + X, weights = log(wrep), family = binomial)
      coefs <- coef(fit.unconstrained.logitistic)
    } else if(outcome.model == "binomial") {
      fit.unconstrained.binomial <- glm(cbind(orep, trialsrep - orep) ~ 1 + X, family = binomial, weights = log(wrep))
      coefs <- coef(fit.unconstrained.binomial)
    } else if(outcome.model == "betabinomial") {
      fit.constrained.betabinomial <- VGAM::vglm(cbind(orep, trialsrep - orep) ~ 1 + X, VGAM::betabinomial, weights = log(wrep))
      coefs <- coef(fit.constrained.betabinomial)[-2]
      coefs.aux <- c(rho = mean(fit.constrained.betabinomial@misc$rho))
    }
  }

  names(coefs) <- gsub("^X", "", names(coefs))
  names(coefs)[names(coefs) == ""] <- "Z"
  names(coefs)[names(coefs) == "(Intercept):1"] <- "(Intercept)"

  return(list(coefs = coefs, coefs.aux = coefs.aux))

}


#' E-step
#'
#' The expectation step in the EM algorithm called by \code{\link{listExperiment}}.
#' 
#' @keywords internal
#'
estep <- function(y, w, x.control, x.sensitive, x.outcome, x.misreport, treat, J,
                  par.sensitive, par.control, par.outcome,
                  par.outcome.aux, par.misreport,
                  d, sensitive.response, model.misreport,
                  o, trials, outcome.model, weight,
                  outcome.constrained, control.constraint, respondentType,
                  misreport.treatment) {

  log.lik <- rep(as.numeric(NA), length(y))

  if(model.misreport == TRUE) {

    # CONTROL ITEMS
    if(control.constraint == "none") {
      hX.misreport.sensitive <-   plogis(cbind(x.control, 1,     sensitive.response) %*% par.control, log.p = TRUE)
      hX.truthful.sensitive <-    plogis(cbind(x.control, 0,     sensitive.response) %*% par.control, log.p = TRUE)
      hX.truthful.nonsensitive <- plogis(cbind(x.control, 0, 1 - sensitive.response) %*% par.control, log.p = TRUE)
    }

    if(control.constraint == "partial") {
      hX.misreport.sensitive <- plogis(cbind(x.control, sensitive.response) %*% par.control, log.p = TRUE)
      hX.truthful.sensitive <- plogis(cbind(x.control, sensitive.response) %*% par.control, log.p = TRUE)
      hX.truthful.nonsensitive <- plogis(cbind(x.control, 1 - sensitive.response) %*% par.control, log.p = TRUE)
    }

    if(control.constraint == "full") {
      hX.misreport.sensitive <- plogis(x.control %*% par.control, log.p = TRUE)
      hX.truthful.sensitive <- plogis(x.control %*% par.control, log.p = TRUE)
      hX.truthful.nonsensitive <- plogis(x.control %*% par.control, log.p = TRUE)
    }

    hX.misreport.sensitive <- dbinom((y - treat * as.numeric(sensitive.response == 1)), size = J, prob = exp(hX.misreport.sensitive), log = TRUE)
    hX.truthful.sensitive <- dbinom((y - treat * as.numeric(sensitive.response == 1)), size = J, prob = exp(hX.truthful.sensitive), log = TRUE)
    hX.truthful.nonsensitive <- dbinom((y - treat * as.numeric(sensitive.response == 0)), size = J, prob = exp(hX.truthful.nonsensitive), log = TRUE)

    # OUTCOME
    if(outcome.model != "none") {
      if(outcome.constrained == TRUE) {
        if(outcome.model %in% c("logistic", "binomial", "betabinomial")) {
          fX.misreport.sensitive <- plogis(cbind(x.outcome, 1) %*% par.outcome, log.p = TRUE)
          fX.truthful.sensitive <- plogis(cbind(x.outcome, 1) %*% par.outcome, log.p = TRUE)
          fX.truthful.nonsensitive <- plogis(cbind(x.outcome, 0) %*% par.outcome, log.p = TRUE)
        }
      } else {
        if(outcome.model %in% c("logistic", "binomial", "betabinomial")) {
          fX.misreport.sensitive <- plogis(cbind(x.outcome, 1, 1) %*% par.outcome, log.p = TRUE)
          fX.truthful.sensitive <- plogis(cbind(x.outcome, 0, 1) %*% par.outcome, log.p = TRUE)
          fX.truthful.nonsensitive <- plogis(cbind(x.outcome, 0, 0) %*% par.outcome, log.p = TRUE)
        }
      }
    } else {
        fX.misreport.sensitive <- rep(0, length(y))
        fX.truthful.sensitive <- rep(0, length(y))
        fX.truthful.nonsensitive <- rep(0, length(y))
    }

    if(outcome.model == "logistic") {
      fX.misreport.sensitive <- dbinom(o, size = 1, prob = exp(fX.misreport.sensitive), log = TRUE)
      fX.truthful.sensitive <- dbinom(o, size = 1, prob = exp(fX.truthful.sensitive), log = TRUE)
      fX.truthful.nonsensitive <- dbinom(o, size = 1, prob = exp(fX.truthful.nonsensitive), log = TRUE)
    } else if(outcome.model == "binomial") {
      fX.misreport.sensitive <- dbinom(o, size = trials, prob = exp(fX.misreport.sensitive), log = TRUE)
      fX.truthful.sensitive <- dbinom(o, size = trials, prob = exp(fX.truthful.sensitive), log = TRUE)
      fX.truthful.nonsensitive <- dbinom(o, size = trials, prob = exp(fX.truthful.nonsensitive), log = TRUE)
    } else if(outcome.model == "betabinomial") {
      fX.misreport.sensitive <- VGAM::dbetabinom(o, size = trials, prob = exp(fX.misreport.sensitive), rho = par.outcome.aux["rho"], log = TRUE)
      fX.truthful.sensitive <- VGAM::dbetabinom(o, size = trials, prob = exp(fX.truthful.sensitive), rho = par.outcome.aux["rho"], log = TRUE)
      fX.truthful.nonsensitive <- VGAM::dbetabinom(o, size = trials, prob = exp(fX.truthful.nonsensitive), rho = par.outcome.aux["rho"], log = TRUE)
    }

    # SENSITIVE ITEM
    if(sensitive.response == 1) {
      gX.misreport.sensitive <- plogis(x.sensitive %*% par.sensitive, log.p = TRUE)
      gX.truthful.sensitive <- plogis(x.sensitive %*% par.sensitive, log.p = TRUE)
      gX.truthful.nonsensitive <- log1p(-exp(plogis(x.sensitive %*% par.sensitive, log.p = TRUE))) # log(1 - plogis(x.sensitive %*% par.sensitive))
    } else {
      gX.misreport.sensitive <- log1p(-exp(plogis(x.sensitive %*% par.sensitive, log.p = TRUE))) # log(1 - plogis(x.sensitive %*% par.sensitive))
      gX.truthful.sensitive <- log1p(-exp(plogis(x.sensitive %*% par.sensitive, log.p = TRUE))) # log(1 - plogis(x.sensitive %*% par.sensitive))
      gX.truthful.nonsensitive <- plogis(x.sensitive %*% par.sensitive, log.p = TRUE)
    }

    # MISREPORTING
    if(misreport.treatment == TRUE) {
      lX.misreport.sensitive <- plogis(cbind(x.misreport, treat) %*% par.misreport, log.p = TRUE)
      lX.truthful.sensitive <- log1p(-exp(plogis(cbind(x.misreport, treat) %*% par.misreport, log.p = TRUE))) # log(1 - exp(plogis(x.misreport %*% par.misreport, log.p = TRUE)))
      lX.truthful.nonsensitive <- log(rep(1, length(y))) # Non-sensitive don't misreport it (monotonicity)
    } else {
      lX.misreport.sensitive <- plogis(x.misreport %*% par.misreport, log.p = TRUE)
      lX.truthful.sensitive <- log1p(-exp(plogis(x.misreport %*% par.misreport, log.p = TRUE))) # log(1 - exp(plogis(x.misreport %*% par.misreport, log.p = TRUE)))
      lX.truthful.nonsensitive <- log(rep(1, length(y))) # Non-sensitive don't misreport it  (monotonicity)
    }

    w[, 1] <- lX.misreport.sensitive +   gX.misreport.sensitive +   hX.misreport.sensitive +   fX.misreport.sensitive
    w[, 2] <- lX.truthful.sensitive +    gX.truthful.sensitive +    hX.truthful.sensitive +    fX.truthful.sensitive
    w[, 3] <- lX.truthful.nonsensitive + gX.truthful.nonsensitive + hX.truthful.nonsensitive + fX.truthful.nonsensitive

    w[respondentType == "Misreport sensitive", 1] <- log(1)
    w[respondentType == "Misreport sensitive", 2] <- log(0)
    w[respondentType == "Misreport sensitive", 3] <- log(0)

    w[respondentType == "Truthful sensitive", 1] <- log(0)
    w[respondentType == "Truthful sensitive", 2] <- log(1)
    w[respondentType == "Truthful sensitive", 3] <- log(0)

    w[respondentType == "Non-sensitive", 1] <- log(0)
    w[respondentType == "Non-sensitive", 2] <- log(0)
    w[respondentType == "Non-sensitive", 3] <- log(1)

    w[respondentType == "Non-sensitive or misreport sensitive", 2] <- log(0)

    denominator <- apply(w, 1, function(x) logAdd(logAdd(x[1], x[2]), x[3]))

    w[, 1] <- w[, 1] - denominator
    w[, 2] <- w[, 2] - denominator
    w[, 3] <- w[, 3] - denominator

    w[respondentType == "Misreport sensitive", 1] <- log(1)
    w[respondentType == "Misreport sensitive", 2] <- log(0)
    w[respondentType == "Misreport sensitive", 3] <- log(0)

    w[respondentType == "Truthful sensitive", 1] <- log(0)
    w[respondentType == "Truthful sensitive", 2] <- log(1)
    w[respondentType == "Truthful sensitive", 3] <- log(0)

    w[respondentType == "Non-sensitive", 1] <- log(0)
    w[respondentType == "Non-sensitive", 2] <- log(0)
    w[respondentType == "Non-sensitive", 3] <- log(1)

    w[respondentType == "Non-sensitive or misreport sensitive", 2] <- log(0)

    # w <- exp(w) / apply(exp(w), 1, sum)

    # Non-sensitive or misreport sensitive
    log.lik[respondentType == "Non-sensitive or misreport sensitive"] <- apply(data.frame(lX.truthful.nonsensitive[respondentType == "Non-sensitive or misreport sensitive"] +
                                                                                          gX.truthful.nonsensitive[respondentType == "Non-sensitive or misreport sensitive"] +
                                                                                          hX.truthful.nonsensitive[respondentType == "Non-sensitive or misreport sensitive"] +
                                                                                          fX.truthful.nonsensitive[respondentType == "Non-sensitive or misreport sensitive"],
                                                                                          lX.misreport.sensitive[respondentType == "Non-sensitive or misreport sensitive"] +
                                                                                          gX.misreport.sensitive[respondentType == "Non-sensitive or misreport sensitive"] +
                                                                                          hX.misreport.sensitive[respondentType == "Non-sensitive or misreport sensitive"] +
                                                                                          fX.misreport.sensitive[respondentType == "Non-sensitive or misreport sensitive"]),
                                                                                  1, function(x) logAdd(x[1], x[2]))

    # Truthful sensitive
    log.lik[respondentType == "Truthful sensitive"] <- lX.truthful.sensitive[respondentType == "Truthful sensitive"] +
                                                       gX.truthful.sensitive[respondentType == "Truthful sensitive"] +
                                                       hX.truthful.sensitive[respondentType == "Truthful sensitive"] +
                                                       fX.truthful.sensitive[respondentType == "Truthful sensitive"]

    # Non-sensitive
    log.lik[respondentType == "Non-sensitive"] <- lX.truthful.nonsensitive[respondentType == "Non-sensitive"] +
                                                  gX.truthful.nonsensitive[respondentType == "Non-sensitive"] +
                                                  hX.truthful.nonsensitive[respondentType == "Non-sensitive"] +
                                                  fX.truthful.nonsensitive[respondentType == "Non-sensitive"]

    # Misreport sensitive
    log.lik[respondentType == "Misreport sensitive"] <- lX.misreport.sensitive[respondentType == "Misreport sensitive"] +
                                                        gX.misreport.sensitive[respondentType == "Misreport sensitive"] +
                                                        hX.misreport.sensitive[respondentType == "Misreport sensitive"] +
                                                        fX.misreport.sensitive[respondentType == "Misreport sensitive"]

  }

  if(model.misreport == FALSE) {

    # CONTROL ITEMS
    if(control.constraint == "none") {
      hX.1 <- plogis(cbind(x.control, 1) %*% par.control)
      hX.0 <- plogis(cbind(x.control, 0) %*% par.control)
    }

    if(control.constraint == "full") {
      hX.1 <- plogis(x.control %*% par.control)
      hX.0 <- plogis(x.control %*% par.control)
    }

    hX.1 <- dbinom((y - treat), size = J, prob = hX.1, log = TRUE)
    hX.0 <- dbinom(y, size = J, prob = hX.0, log = TRUE)

    # OUTCOME
    if(outcome.model %in% c("logistic", "binomial", "betabinomial")) {
      fX.1 <- plogis(cbind(x.outcome, 1) %*% par.outcome)
      fX.0 <- plogis(cbind(x.outcome, 0) %*% par.outcome)
    } else {
      fX.1 <- rep(0, length(y))
      fX.0 <- rep(0, length(y))
    }

    if(outcome.model == "logistic") {
      fX.1 <- dbinom(o, size = 1, prob = fX.1, log = TRUE)
      fX.0 <- dbinom(o, size = 1, prob = fX.0, log = TRUE)
    } else if(outcome.model == "binomial") {
      fX.1 <- dbinom(o, size = trials, prob = fX.1, log = TRUE)
      fX.0 <- dbinom(o, size = trials, prob = fX.0, log = TRUE)
    } else if(outcome.model == "betabinomial") {
      fX.1 <- VGAM::dbetabinom(o, size = trials, prob = fX.1, rho = par.outcome.aux["rho"], log = TRUE)
      fX.0 <- VGAM::dbetabinom(o, size = trials, prob = fX.0, rho = par.outcome.aux["rho"], log = TRUE)
    }

    # SENSITIVE ITEM
    gX.1 <- plogis(x.sensitive %*% par.sensitive, log.p = TRUE)
    gX.0 <- log(1 - exp(gX.1))

    w[, 1] <- gX.1 + hX.1 + fX.1
    w[, 2] <- gX.0 + hX.0 + fX.0

    w[respondentType == "1", 1] <- log(1)
    w[respondentType == "1", 2] <- log(0)

    w[respondentType == "0", 1] <- log(0)
    w[respondentType == "0", 2] <- log(1)

    denominator <- apply(w, 1, function(x) logAdd(x[1], x[2]))

    w[, 1] <- w[, 1] - denominator
    w[, 2] <- w[, 2] - denominator

    w[respondentType == "1", 1] <- log(1)
    w[respondentType == "1", 2] <- log(0)

    w[respondentType == "0", 1] <- log(0)
    w[respondentType == "0", 2] <- log(1)

    # Log likelihood
    log.lik[respondentType == "0"] <- gX.0[respondentType == "0"] +
                                      hX.0[respondentType == "0"] +
                                      fX.0[respondentType == "0"]

    log.lik[respondentType == "1"] <- gX.1[respondentType == "1"] +
                                      hX.1[respondentType == "1"] +
                                      fX.1[respondentType == "1"]

    log.lik[respondentType == "0 or 1"] <- apply(data.frame(gX.1[respondentType == "0 or 1"] +
                                                            hX.1[respondentType == "0 or 1"] +
                                                            fX.1[respondentType == "0 or 1"],
                                                            gX.0[respondentType == "0 or 1"] +
                                                            hX.0[respondentType == "0 or 1"] +
                                                            fX.0[respondentType == "0 or 1"]),
                                                 1, function(x) logAdd(x[1], x[2]))

    log.lik[respondentType == "0 or 1"] <- log(exp(gX.1[respondentType == "0 or 1"] +
                                                   hX.1[respondentType == "0 or 1"] +
                                                   fX.1[respondentType == "0 or 1"]) +
                                               exp(gX.0[respondentType == "0 or 1"] +
                                                   hX.0[respondentType == "0 or 1"] +
                                                   fX.0[respondentType == "0 or 1"]))
  }

  return(list(w = w, ll = sum(weight * log.lik)))
}


#' List experiment regression
#'
#' Regression analysis for sensitive survey questions using a list experiment and direct question.
#' 
#' @param formula An object of class "\code{\link{formula}}": a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables to be used in the model.
#' @param treatment A string indicating the name of the treatment indicator in the data. This variable must be coded as a binary, where 1 indicates assignment to treatment and 0 indicates assignment to control.
#' @param J An integer indicating the number of control items in the list experiment.
#' @param direct A string indicating the name of the direct question response in the data. The direct question must be coded as a binary variable. If NULL (default), a misreport sub-model is not fit.
#' @param sensitive.response A value 0 or 1 indicating whether the response that is considered sensitive in the list experiment/direct question is 0 or 1.
#' @param outcome A string indicating the variable name in the data to use as the outcome in an outcome sub-model. If NULL (default), no outcome sub-model is fit. [\emph{experimental}]
#' @param outcome.trials An integer indicating the number of trials in a binomial/betabinomial model if both an outcome sub-model is used and if the argument \code{outcome.model} is set to "binomial" or "betabinomial". [\emph{experimental}]
#' @param outcome.model A string indicating the model type to fit for the outcome sub-model ("logistic", "binomial", "betabinomial"). [\emph{experimental}]
#' @param outcome.constrained A logical value indicating whether to constrain U* = 0 in the outcome sub-model. Defaults to TRUE. [\emph{experimental}]
#' @param control.constraint A string indicating the constraint to place on Z* and U* in the control-items sub-model:
#'    \describe{
#'         \item{"none" (default)}{Estimate separate parameters for Z* and U*.}
#'         \item{"partial"}{Constrain U* = 0.}
#'         \item{"full"}{Constrain U* = Z* = 0.}
#'     }
#' @param misreport.treatment A logical value indicating whether to include a parameter for the treatment indicator in the misreport sub-model. Defaults to TRUE.
#' @param weights A string indicating the variable name of survey weights in the data (note: standard errors are not currently output when survey weights are used).
#' @param se A logical value indicating whether to calculate standard errors. Defaults to TRUE.
#' @param tolerance The desired accuracy for EM convergence. The EM loop breaks after the change in the log-likelihood is less than the value of \code{tolerance}. Defaults to 1e-08.
#' @param max.iter The maximum number of iterations for the EM algorithm. Defaults to 10000.
#' @param n.runs The total number of times that the EM algorithm is run (can potentially help avoid local maxima). Defaults to 1.
#' @param verbose A logical value indicating whether to print information during model fitting. Defaults to TRUE.
#' @param get.data For internal use. Used by wrapper function \code{\link{bootListExperiment}}.
#' @param par.control A vector of starting parameters for the control-items sub-model. Must be in the order of the parameters in the resulting regression output. If NULL (default), randomly generated starting points are used, drawn from uniform(-2, 2).
#' @param par.sensitive A vector of starting parameters for the sensitive-item sub-model. Must be in the order of the parameters in the resulting regression output. If NULL (default), randomly generated starting points are used, drawn from uniform(-2, 2).
#' @param par.misreport A vector of starting parameters for the misreport sub-model. Must be in the order of the parameters in the resulting regression output. If NULL (default), randomly generated starting points are used, drawn from uniform(-2, 2).
#' @param par.outcome A vector of starting parameters for the outcome sub-model.  Must be in the order of the parameters in the resulting regression output. If NULL (default), randomly generated starting points are used, drawn from uniform(-2, 2). [experimental]
#' @param par.outcome.aux A vector of starting parameters for the outcome sub-model in which \code{outcome.model} is "betabinomial". i.e. c(alpha, beta). If NULL (default), randomly generated starting points are used, drawn from uniform(0, 1). [experimental]
#' @param formula.control An object of class "\code{\link{formula}}" used to specify a control-items sub-model that is different from that given in \code{formula}. (e.g. ~ x1 + x2)
#' @param formula.sensitive An object of class "\code{\link{formula}}" used to specify a sensitive-item sub-model that is different from that given in \code{formula}. (e.g. ~ x1 + x2)
#' @param formula.misreport An object of class "\code{\link{formula}}" used to specify a misreport sub-model that is different from that given in \code{formula}. (e.g. ~ x1 + x2)
#' @param formula.outcome An object of class "\code{\link{formula}}" used to specify an outcome sub-model that is different from that given in \code{formula}. (e.g. ~ x1 + x2) [\emph{experimental}]
#' @param get.boot For internal use. An integer, which if greater than 0 requests that \code{listExperiment()} generate a non-parametric bootstrap sample and fit a model to that sample. Used by the function \code{\link{bootListExperiment}}.
#' @param ... Additional options.
#' 
#' @details The \code{listExperiment} function allows researchers to fit a model
#'          for a list experiment and direct question simultaneously, as
#'          described in Eady (2016). The primary aim of the function is
#'          to allow researchers to model the probability that respondents
#'          provides one response to the sensitive item in a list experiment
#'          but respond otherwise when asked about the same sensitive item on a
#'          direct question. When a direct question response is excluded from
#'          the function, the model is functionally equivalent to that proposed
#'          by Imai (2011), as implemented as the \code{\link[list]{ictreg}} function
#'          in the \code{list} package (\url{https://CRAN.R-project.org/package=list}).
#' 
#' @return \code{listExperiment} returns an object of class "listExperiment".
#'     A summary of this object is given using the \code{\link{summary.listExperiment}}
#'     function. All components in the "listExperiment" class are listed below.
#' @slot par.control A named vector of coefficients from the control-items sub-model.
#' @slot par.sensitive A named vector of coefficients from the sensitive-item sub-model.
#' @slot par.misreport A named vector of coefficients from the misreport sub-model.
#' @slot par.outcome A named vector of coefficients from the outcome sub-model.
#' @slot par.outcome.aux A named vector of (auxiliary) coefficients from the outcome sub-model (if \code{outcome.model} = "betabinomial").
#' @slot df Degrees of freedom.
#' @slot se.sensitive Standard errors for parameters in the sensitive-item sub-model.
#' @slot se.control Standard errors for parameters in the control-items sub-model.
#' @slot se.misreport Standard errors for parameters in the misreport sub-model.
#' @slot se.outcome Standard errors for parameters in the outcome sub-model.
#' @slot se.outcome.aux Standard errors for the auxiliary parameters in the outcome sub-model (if \code{outcome.model} = "betabinomial").
#' @slot vcov.mle Variance-covariance matrix.
#' @slot w The matrix of posterior predicted probabilities for each observation in the data used for model fitting.
#' @slot data The data frame used for model fitting.
#' @slot direct The string indicating the variable name of the direct question.
#' @slot treatment The string indicating the variable name of the treatment indicator.
#' @slot model.misreport A logical value indicating whether a misreport sub-model was fit.
#' @slot outcome.model The type of model used as the outcome sub-model.
#' @slot outcome.constrained A logical value indicating whether the parameter U* was constrained to 0 in the outcome sub-model.
#' @slot control.constraint A string indicating the constraints placed on the parameters Z* and U* in the control-items sub-model.
#' @slot misreport.treatment A logical value indicating whether a treatment indicator was included in the misreport sub-model.
#' @slot weights A string indicating the variable name of the survey weights.
#' @slot formula The model formula.
#' @slot formula.control The model specification of the control-items sub-model.
#' @slot formula.sensitive The model specification of the sensitive-item sub-model.
#' @slot formula.misreport The model specification of the misreport sub-model.
#' @slot formula.outcome The model specification of the outcome sub-model.
#' @slot sensitive.response The value 0 or 1 indicating the response to the list experiment/direct question that is considered sensitive.
#' @slot xlevels The factor levels of the variables used in the model.
#' @slot llik The model log-likelihood.
#' @slot n The sample size of the data used for model fitting (this value excludes rows removed through listwise deletion).
#' @slot J The number of control items in the list experiment.
#' @slot se A logical value indicating whether standard errors were calculated.
#' @slot runs The parameter estimates from each run of the EM algorithm (note: the parameters that result in the highest log-likelihood are used as the model solution).
#' @slot call The method call.
#' @slot boot A logical value indicating whether non-parametric bootstrapping was used to calculate model parameters and standard errors.
#' 
#' @references Eady, Gregory. 2016 "The Statistical Analysis of Misreporting on Sensitive Survey Questions."
#' @references Imai, Kosuke. 2011. "Multivariate Regression Analysis for the Item Count Technique." \emph{Journal of the American Statistical Association} 106 (494): 407-416.
#' 
#' @examples
#'
#' ## EXAMPLE 1: Simulated list experiment and direct question
#' n <- 10000
#' J <- 4
#'
#' # Covariates
#' x <- cbind(intercept = rep(1, n), continuous1 = rnorm(n),
#'            continuous2 = rnorm(n), binary1 = rbinom(n, 1, 0.5))
#'
#' treatment <- rbinom(n, 1, 0.5)
#'
#' # Simulate Z*
#' param_sensitive <- c(0.25, -0.25, 0.5, 0.25)
#' prob_sensitive <- plogis(x %*% param_sensitive)
#' true_belief <- rbinom(n, 1, prob = prob_sensitive)
#'
#' # Simulate whether respondent misreports (U*)
#' param_misreport <- c(-0.25, 0.25, -0.5, 0.5)
#' prob_misreport <- plogis(x %*% param_misreport) * true_belief
#' misreport <- rbinom(n, 1, prob = prob_misreport)
#'
#' # Simulate control items Y*
#' param_control <- c(0.25, 0.25, -0.25, 0.25, U = -0.5, Z = 0.25)
#' prob.control <- plogis(cbind(x, misreport, true_belief) %*% param_control)
#' control_items <- rbinom(n, J, prob.control)
#'
#' # List experiment and direct question responses
#' direct <- true_belief
#' direct[misreport == 1] <- 0
#' y <- control_items + true_belief * treatment
#'
#' A <- data.frame(y, direct, treatment,
#'                 continuous1 = x[, "continuous1"],
#'                 continuous2 = x[, "continuous2"],
#'                 binary1 = x[, "binary1"])
#'
#' \dontrun{
#' model.sim <- listExperiment(y ~ continuous1 + continuous2 + binary1,
#'                             data = A, treatment = "treatment", direct = "direct",
#'                             J = 4, control.constraint = "none",
#'                             sensitive.response = 1)
#' summary(model.sim, digits = 3)
#' }
#'
#'
#' ## EXAMPLE 2: Data from Eady (2016)
#' data(gender)
#'
#' \dontrun{
#' # Note: substantial computation time
#' model.gender <- listExperiment(y ~ gender + ageGroup + education +
#'                                    motherTongue + region + selfPlacement,
#'                                data = gender, J = 4,
#'                                treatment = "treatment", direct = "direct",
#'                                control.constraint = "none",
#'                                sensitive.response = 0,
#'                                misreport.treatment = TRUE)
#' summary(model.gender)
#' }
#'
#' @export
#' 
listExperiment <- function(formula, data, treatment, J, 
                           direct = NULL, sensitive.response = NULL,
                           outcome = NULL, outcome.trials = NULL,
                           outcome.model = "logistic",
                           outcome.constrained = TRUE,
                           control.constraint = "none",
                           misreport.treatment = TRUE,
                           weights = NULL, se = TRUE, tolerance = 1E-8,
                           max.iter = 10000, n.runs = 1, verbose = TRUE,
                           get.data = FALSE,
                           par.control = NULL, par.sensitive = NULL,
                           par.misreport = NULL, par.outcome = NULL,
                           par.outcome.aux = NULL,
                           formula.control = NULL, formula.sensitive = NULL,
                           formula.misreport = NULL, formula.outcome = NULL,
                           get.boot = 0, ...) {

  function.call <- match.call(expand.dots = FALSE)

  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)

  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.pass"
  mf[[1]] <- quote(model.frame)
  mt <- attr(eval(mf, parent.frame()), "terms")
  xlevels.formula <- .getXlevels(attr(eval(mf, parent.frame()), "terms"), eval(mf, parent.frame()))

  if(!is.null(formula.control)) {
    mf.control <- mf
    mf.control$formula <- formula.control
    xlevels.formula.control <- .getXlevels(attr(eval(mf.control, parent.frame()), "terms"), eval(mf.control, parent.frame()))
    mf.control <- eval(mf.control, parent.frame())
    x.control <- model.matrix(attr(mf.control, "terms"), data = mf.control)
  } else {
    formula.control <- as.formula(mf$formula)
    xlevels.formula.control <- xlevels.formula
    mf.control <- eval(mf, parent.frame())
    x.control <- model.matrix(attr(mf.control, "terms"), data = mf.control)
  }

  if(!is.null(formula.sensitive)) {
    mf.sensitive <- mf
    mf.sensitive$formula <- formula.sensitive
    xlevels.formula.sensitive <- .getXlevels(attr(eval(mf.sensitive, parent.frame()), "terms"), eval(mf.sensitive, parent.frame()))
    mf.sensitive <- eval(mf.sensitive, parent.frame())
    x.sensitive <- model.matrix(attr(mf.sensitive, "terms"), data = mf.sensitive)
  } else {
    formula.sensitive <- as.formula(mf$formula)
    xlevels.formula.sensitive <- xlevels.formula
    mf.sensitive <- eval(mf, parent.frame())
    x.sensitive <- model.matrix(attr(mf.sensitive, "terms"), data = mf.sensitive)
  }

  if(!is.null(formula.misreport)) {
    mf.misreport <- mf
    mf.misreport$formula <- formula.misreport
    xlevels.formula.misreport <- .getXlevels(attr(eval(mf.misreport, parent.frame()), "terms"), eval(mf.misreport, parent.frame()))
    mf.misreport <- eval(mf.misreport, parent.frame())
    x.misreport <- model.matrix(attr(mf.misreport, "terms"), data = mf.misreport)
  } else {
    formula.misreport <- as.formula(mf$formula)
    xlevels.formula.misreport <- xlevels.formula
    mf.misreport <- eval(mf, parent.frame())
    x.misreport <- model.matrix(attr(mf.misreport, "terms"), data = mf.misreport)
  }

  if(!is.null(formula.outcome)) {
    mf.outcome <- mf
    mf.outcome$formula <- formula.outcome
    xlevels.formula.outcome <- .getXlevels(attr(eval(mf.outcome, parent.frame()), "terms"), eval(mf.outcome, parent.frame()))
    mf.outcome <- eval(mf.outcome, parent.frame())
    x.outcome <- model.matrix(attr(mf.outcome, "terms"), data = mf.outcome)
   } else {
    formula.outcome <- as.formula(mf$formula)
    xlevels.formula.outcome <- xlevels.formula
    mf.outcome <- eval(mf, parent.frame())
    x.outcome <- model.matrix(attr(mf.outcome, "terms"), data = mf.outcome)
  }

  mf <- eval(mf, parent.frame())
  y <- model.response(mf, type = "any")
  treat <- data[, paste(treatment)]

  xlevels <- c(xlevels.formula,
               xlevels.formula.control,
               xlevels.formula.sensitive,
               xlevels.formula.misreport,
               xlevels.formula.outcome)
  xlevels <- xlevels[-which(duplicated(xlevels))]


  # x.na <- apply(x, 1, function(X) all(!is.na(X)))
  x.control.na <- apply(x.control, 1, function(X) all(!is.na(X)))
  x.sensitive.na <- apply(x.sensitive, 1, function(X) all(!is.na(X)))
  x.misreport.na <- apply(x.misreport, 1, function(X) all(!is.na(X)))
  x.outcome.na <- apply(x.outcome, 1, function(X) all(!is.na(X)))

  y.na <- !is.na(y)
  treat.na <- !is.na(treat)

  if(!is.null(direct)) {
    d <- data[, paste(direct)]
    d.na <- !is.na(d)
    model.misreport <- TRUE
  } else {
    model.misreport <- FALSE
    d <- rep(NA, length(y))
    d.na <- rep(TRUE, length(y))
  }

  if(!is.null(outcome) & outcome.model %in% c("logistic")) {
    o <- data[, paste(outcome)]
    trials <- rep(NA, length(y))
    o.na <- !is.na(o)
  } else if(!is.null(outcome) & outcome.model %in% c("binomial", "betabinomial")) {
    o <- data[, paste(outcome)]
    trials <- data[, paste(outcome.trials)]
    o.na <- !is.na(o) & !is.na(trials)
  } else {
    o <- rep(NA, length(y))
    trials <- rep(NA, length(y))
    o.na <- rep(TRUE, length(y))
    outcome.model <- "none"
  }

  if(!is.null(weights)) {
    weight <- data[, paste(weights)]
    weight.na <- !is.na(weight)
  } else {
    weight <- rep(1, length(y))
    weight.na <- !is.na(weight)
  }

  y <- y[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na]
  x.control <- as.matrix(x.control[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na, , drop = FALSE])
  x.sensitive <- as.matrix(x.sensitive[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na, , drop = FALSE])
  x.outcome <- as.matrix(x.outcome[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na, , drop = FALSE])
  x.misreport <- as.matrix(x.misreport[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na, , drop = FALSE])
  treat <- treat[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na]
  d <- d[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na]
  o <- o[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na]
  trials <- trials[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na]
  weight <- weight[y.na & x.control.na & x.sensitive.na & x.outcome.na & x.misreport.na & treat.na & d.na & weight.na]

  n <- nrow(x.control)

  # For testing whether arguments are correctly interpreted:

  # return(list(y = y, x.control = x.control, x.sensitive = x.sensitive,
  #             x.outcome = x.outcome, x.misreport = x.misreport, treat = treat,
  #             d = d, o = o, trials = trials, weight = weight,
  #             control.constraint = control.constraint, misreport.treatment = misreport.treatment,
  #             model.misreport = model.misreport, outcome.model = outcome.model, se = se,
  #             sensitive.response = sensitive.response, J = J,
  #             par.control = par.control, par.sensitive = par.sensitive,
  #             par.outcome = par.outcome, par.outcome.aux = par.outcome.aux, par.misreport = par.misreport))}

  # y <- model$y
  # x.control <- model$x.control
  # x.sensitive <- model$x.sensitive
  # x.outcome <- model$x.outcome
  # x.misreport <- model$x.misreport
  # treat <- model$treat
  # d <- model$d
  # o <- model$o
  # trials <- model$trials
  # weight <- model$weight
  # control.constraint <- model$control.constraint
  # misreport.treatment <- model$misreport.treatment
  # model.misreport <- model$model.misreport
  # outcome.model <- model$outcome.model
  # se <- model$se
  # sensitive.response <- model$sensitive.response
  # J <- model$J
  # max.iter <- 2000
  # par.control <- NULL
  # par.sensitive <- NULL
  # par.outcome <- NULL
  # par.outcome.aux <- NULL
  # par.misreport <- NULL
  ###########
  # max.iter <- 5000
  # verbose <- TRUE
  # tolerance <- 1E-08
  # j <- 1
  # i <- 1
  # n.runs <- 1
  # get.boot <- 0
  # get.data <- FALSE

  # Draw a non-parametric boot-strap sample if
  # requested by the bootListExperiment wrapper
  if(get.boot > 0) {
    boot.sample <- sample(1:length(weight), prob = weight, replace = TRUE)
    y <- as.matrix(y)[boot.sample, , drop = FALSE]
    x.control <- as.matrix(x.control)[boot.sample, , drop = FALSE]
    x.sensitive <- as.matrix(x.sensitive)[boot.sample, , drop = FALSE]
    x.outcome <- as.matrix(x.outcome)[boot.sample, , drop = FALSE]
    x.misreport <- as.matrix(x.misreport)[boot.sample, , drop = FALSE]
    treat <- as.matrix(treat)[boot.sample]
    d <- as.matrix(d)[boot.sample, , drop = FALSE]
    o <- as.matrix(o)[boot.sample, , drop = FALSE]
    trials <- as.matrix(trials)[boot.sample, , drop = FALSE]
    weight <- rep(1, length(y))
    se <- FALSE
  }

  respondentType <- rep(as.character(NA), length(y))

  if(model.misreport == TRUE) {

    # Treat == 0, Y == control only, direct == Non-sensitive
    respondentType[treat == 0 & d != sensitive.response] <- "Non-sensitive or misreport sensitive"

    # Treat == 0, Y == control only, direct == Sensitive
    respondentType[treat == 0 & d == sensitive.response] <- "Truthful sensitive"

    # Treat == 1, Y == (J+1) or 0, direct == Non-sensitive
    if(sensitive.response == 1) respondentType[treat == 1 & y == 0 & d != sensitive.response] <- "Non-sensitive"
    if(sensitive.response == 0) respondentType[treat == 1 & y == (J + 1) & d != sensitive.response] <- "Non-sensitive"

    # Treat == 1, 0 < Y < (J + 1), direct == Non-sensitive
    respondentType[treat == 1 & y > 0 & y < (J + 1) & d != sensitive.response] <- "Non-sensitive or misreport sensitive"

    # Treat == 1, Y == (J + 1) or 0, direct == Non-sensitive
    if(sensitive.response == 1) respondentType[treat == 1 & y == (J + 1) & d != sensitive.response] <- "Misreport sensitive"
    if(sensitive.response == 0) respondentType[treat == 1 & y == 0 & d != sensitive.response] <- "Misreport sensitive"

    # Treat == 1, Y == (J + 1) or 0, direct == Sensitive
    if(sensitive.response == 1) respondentType[treat == 1 & y == (J + 1) & d == sensitive.response] <- "Truthful sensitive"
    if(sensitive.response == 0) respondentType[treat == 1 & y == 0 & d == sensitive.response] <- "Truthful sensitive"

    # Treat == 1, Y == (J + 1) or 0, direct == Sensitive (not possible by assumption; error check for this)
    if(sensitive.response == 1) respondentType[treat == 1 & y == 0 & d == sensitive.response] <- "Violates assumption"
    if(sensitive.response == 0) respondentType[treat == 1 & y == (J + 1) & d == sensitive.response] <- "Violates assumption"

    # Treat == 1, 0 < Y < (J + 1), direct == Sensitive
    respondentType[treat == 1 & y > 0 & y < (J + 1) & d == sensitive.response] <- "Truthful sensitive"
  } else {
    # Treat == 1 0 < Y < (J + 1) is a "0 or a 1"
    respondentType[treat == 1 & y > 0 & y < (J + 1)] <- "0 or 1"

    # Treat == 0 Y == 0 is a 0 or a "1"
    respondentType[treat == 0] <- "0 or 1"

    # Treat == 1 Y == 0 is a "0"
    respondentType[treat == 1 & y == 0] <- "0"

    # Treat == 1 Y == (J + 1) is a "1"
    respondentType[(treat == 1 & y == (J + 1))] <- "1"
  }

  if("Violates assumption" %in% respondentType) {
    stop("\nSome observations violate the monotonicity assumption.")
  }

  # SET UP THE POSTERIOR PROBABILITIES
  if(model.misreport == TRUE) {
    w <- as.matrix(data.frame(as.numeric(respondentType %in% c("Non-sensitive or misreport sensitive", "Misreport sensitive")),
                              as.numeric(respondentType == "Truthful sensitive"),
                              as.numeric(respondentType %in% c("Non-sensitive or misreport sensitive", "Non-sensitive"))))
    w <- w / apply(w, 1, sum)
    colnames(w) <- c("Misreport sensitive", "Truthful sensitive", "Non-sensitive")
  } else {
    w <- as.matrix(data.frame(as.numeric(respondentType %in% c("1", "0 or 1")),
                              as.numeric(respondentType %in% c("0", "0 or 1"))))
    w <- w / apply(w, 1, sum)
    colnames(w) <- c("1", "0")
  }

  w <- log(w)

  if(get.data == TRUE) {

    estep.out <- estep(y = y, w = w, x.control = x.control, x.sensitive = x.sensitive, x.outcome = x.outcome, x.misreport = x.misreport, treat = treat, J = J,
                       par.sensitive = par.sensitive, par.control = par.control, par.outcome = par.outcome, par.outcome.aux = par.outcome.aux, par.misreport = par.misreport,
                       d = d, sensitive.response = sensitive.response, model.misreport = model.misreport,
                       o = o, trials = trials, outcome.model = outcome.model, 
                       weight = weight, respondentType = respondentType,
                       outcome.constrained = outcome.constrained,
                       control.constraint = control.constraint,
                       misreport.treatment = misreport.treatment)

    return(list(w = estep.out$w,
                ll = estep.out$ll,
                x.control = x.control,
                x.sensitive = x.sensitive,
                x.misreport = x.misreport,
                x.outcome = x.outcome))

  }

  if(model.misreport == TRUE) {
    if(misreport.treatment == TRUE) par.misreport <- rep(0, ncol(x.misreport) + 1) # +1 for treatment (consistency bias)
    if(misreport.treatment == FALSE) par.misreport <- rep(0, ncol(x.misreport))
  } else {
    par.misreport <- NULL
  }

  if(is.null(par.sensitive)) par.sensitive <- rep(0, ncol(x.sensitive))

  if(is.null(par.control)) {
    if(control.constraint == "none" & model.misreport == FALSE) {
      par.control <- rep(0, ncol(x.control) + 1)
    } else if(control.constraint == "none" & model.misreport == TRUE) {
      par.control <- rep(0, ncol(x.control) + 2)
    } else if(control.constraint == "partial" & model.misreport == FALSE) {
      stop("If not modeling misreporting, set argument control.constraint to 'none' or 'full'")
    } else if(control.constraint == "partial" & model.misreport == TRUE) {
      par.control <- rep(0, ncol(x.control) + 1)
    } else if(control.constraint == "full") {
      par.control <- rep(0, ncol(x.control))
    }
  }

  if(is.null(par.outcome)) {
    if(outcome.model != "none") {
      if(outcome.constrained == TRUE) par.outcome <- rep(0, ncol(x.outcome) + 1)
      if(outcome.constrained == FALSE) par.outcome <- rep(0, ncol(x.outcome) + 2)
    } else {
      par.outcome <- NULL
    }
  }

  if(is.null(par.outcome.aux)) {
    if(outcome.model %in% c("none", "logistic")) {
      par.outcome.aux <- NULL
    } else if(outcome.model == "betabinomial") {
      par.outcome.aux <- list(rho = 0)
    }
  }
  
  runs <- list()

  # EXPECTATION MAXIMIZATION LOOP
  for(j in 1:n.runs) {
    if(j > 1 & verbose == TRUE) cat("\n")

    logLikelihood <- rep(as.numeric(NA), max.iter)

    # Get starting points on uniform(-2, 2)
    while(TRUE) {
      par.control <- runif(length(par.control), -2, 2)
      par.sensitive <- runif(length(par.sensitive), -2, 2)
      if(model.misreport == TRUE) par.misreport <- runif(length(par.misreport), -2, 2)
      if(outcome.model != "none") par.outcome <- runif(length(par.outcome), -2, 2)
      if(outcome.model != "none" & length(par.outcome.aux) > 0) par.outcome.aux <- runif(length(par.outcome.aux), 0, 1)

      templl <- estep(y = y, w = w, x.control = x.control, x.sensitive = x.sensitive, x.outcome = x.outcome, x.misreport = x.misreport, treat = treat, J = J,
                      par.sensitive = par.sensitive, par.control = par.control, par.outcome = par.outcome, par.outcome.aux = par.outcome.aux, par.misreport = par.misreport,
                      d = d, sensitive.response = sensitive.response, model.misreport = model.misreport,
                      o = o, trials = trials, outcome.model = outcome.model, 
                      weight = weight, respondentType = respondentType,
                      outcome.constrained = outcome.constrained,
                      control.constraint = control.constraint,
                      misreport.treatment = misreport.treatment)$ll
      templl

      if(!is.nan(templl) & templl > -Inf) break()
    }

    for(i in 1:max.iter) { # E-M loop

      estep.out <- estep(y = y, w = w, x.control = x.control, x.sensitive = x.sensitive, x.outcome = x.outcome, x.misreport = x.misreport, treat = treat, J = J,
                         par.sensitive = par.sensitive, par.control = par.control, par.outcome = par.outcome, par.outcome.aux = par.outcome.aux, par.misreport = par.misreport,
                         d = d, sensitive.response = sensitive.response, model.misreport = model.misreport,
                         o = o, trials = trials, outcome.model = outcome.model, 
                         weight = weight, respondentType = respondentType,
                         outcome.constrained = outcome.constrained,
                         control.constraint = control.constraint,
                         misreport.treatment = misreport.treatment)

      w <- estep.out$w
      logLikelihood[i] <- estep.out$ll

      if(i > 1 & verbose == TRUE & get.boot == 0) {
        cat("\r\rRun:", paste0(j, "/", n.runs), "Iter:", i,
            "llik:", sprintf("%.2f", logLikelihood[i]),
            "llik change:", sprintf("%.8f", (logLikelihood[i] - logLikelihood[i-1])),
            "(tol =", paste0(as.character(tolerance), ")       "))
      }

      if(i > 1 & verbose == TRUE & get.boot > 0) {
        cat("\r\rBoot:", get.boot, "Run:", paste0(j, "/", n.runs), "Iter:", i,
            "llik:", sprintf("%.2f", logLikelihood[i]),
            "llik change:", sprintf("%.8f", (logLikelihood[i] - logLikelihood[i-1])),
            "(tol =", paste0(as.character(tolerance), ")       "))
      }

      if(i > 1 && (logLikelihood[i] - logLikelihood[i - 1]) < 0) {
        stop("Log-likelihood increasing.")
      }
      if(i > 1 && (logLikelihood[i] - logLikelihood[i - 1]) < tolerance) {
        break()
      }

      par.sensitive <- mstepSensitive(y = y, treat = treat, x.sensitive = x.sensitive, w = w,
                                      d = d, sensitive.response = sensitive.response,
                                      weight = weight, model.misreport = model.misreport)

      par.control <- mstepControl(y = y, J = J, treat = treat, x.control = x.control, w = w,
                                  d = d, sensitive.response = sensitive.response,
                                  weight = weight, model.misreport = model.misreport,
                                  control.constraint = control.constraint)

      if(outcome.model != "none") {
        outcome <- mstepOutcome(y = y, treat = treat, x.outcome = x.outcome, w = w,
                                d = d, sensitive.response = sensitive.response,
                                o = o, trials = trials, weight = weight,
                                model.misreport = model.misreport,
                                outcome.model = outcome.model,
                                outcome.constrained = outcome.constrained,
                                control.constraint = control.constraint)
        par.outcome <- outcome$coefs
        par.outcome.aux <- outcome$coefs.aux
      }

      if(model.misreport == TRUE) {
        par.misreport <- mstepMisreport(y = y, x.misreport = x.misreport,
                                              w = w, treat = treat,
                                              misreport.treatment = misreport.treatment,
                                              weight = weight)
      }
    }
  
    runs[[j]] <- list(logLikelihood = logLikelihood[i],
                      par.control = par.control,
                      par.sensitive = par.sensitive,
                      par.misreport = par.misreport,
                      par.outcome = par.outcome,
                      par.outcome.aux = par.outcome.aux)

  }

  if(verbose == TRUE) cat("\n")

  max.ll <- which(sapply(runs, function(X) X$logLikelihood) == max(sapply(runs, function(X) X$logLikelihood)))
  llik <- runs[[max.ll]]$logLikelihood
  par.control <- runs[[max.ll]]$par.control
  par.sensitive <- runs[[max.ll]]$par.sensitive
  par.misreport <- runs[[max.ll]]$par.misreport
  par.outcome <- runs[[max.ll]]$par.outcome
  par.outcome.aux <- runs[[max.ll]]$par.outcome.aux

  par <- c(par.control, par.sensitive, par.misreport, par.outcome, par.outcome.aux)
  num <- c(length(par.control), length(par.sensitive), length(par.misreport), length(par.outcome), length(par.outcome.aux))

  llik.wrapper <- function(par, num, y, w,
                           x.control, x.sensitive, x.outcome, x.misreport, treat, J,
                           d, sensitive.response, model.misreport,
                           o, trials, outcome.model,
                           weight, respondentType,
                           outcome.constrained,
                           control.constraint,
                           misreport.treatment) {
    par.control <- par[1:num[1]]
    par.sensitive <- par[(num[1]+1):sum(num[1:2])]
    if(model.misreport == TRUE) {
      par.misreport <- par[(sum(num[1:2])+1):sum(num[1:3])]
    } else{
      par.misreport <- NULL
    }
    if(outcome.model != "none") {
      par.outcome <- par[(sum(num[1:3])+1):sum(num[1:4])]
      if(outcome.model %in% c("betabinomial", "linear")) {
        par.outcome.aux <- par[(sum(num[1:4])+1):sum(num[1:5])]
      } else {
        par.outcome.aux <- NULL
      }
    } else {
      par.outcome <- NULL
    }

    llik <- estep(y = y, w = w, x.control = x.control, x.sensitive = x.sensitive, x.outcome = x.outcome, x.misreport = x.misreport, treat = treat, J = J,
                  par.sensitive = par.sensitive, par.control = par.control, par.outcome = par.outcome, par.outcome.aux = par.outcome.aux, par.misreport = par.misreport,
                  d = d, sensitive.response = sensitive.response, model.misreport = model.misreport,
                  o = o, trials = trials, outcome.model = outcome.model,
                  weight = weight, respondentType = respondentType,
                  outcome.constrained = outcome.constrained,
                  control.constraint = control.constraint,
                  misreport.treatment)$ll
    return(llik)
  }

# For testing:
# par = c(par.control, par.sensitive, par.misreport, par.outcome, par.outcome.aux)
# num = c(length(par.control), length(par.sensitive), length(par.misreport), length(par.outcome))
# J = J
# y = y
# treat = treat
# x = x
# x.misreport = x.misreport
# d = d
# sensitive.response = sensitive.response
# model.misreport = model.misreport
# o = o
# outcome.model = outcome.model
# weight = weight
# respondentType = respondentType
# control.constraint = control.constraint

# llik.wrapper(par = par, num = num, y = y, 
#                  x.control = x.control, x.sensitive = x.sensitive,
#                  x.outcome = x.outcome, x.misreport = x.misreport, treat = treat,
#                  J = J, d = d, sensitive.response = sensitive.response,
#                  model.misreport = model.misreport, o = o, outcome.model = outcome.model,
#                  outcome.constrained = outcome.constrained, weight = weight, respondentType = respondentType,
#                  control.constraint = control.constraint)

  if(se == TRUE & all(weight == 1)) {
    # hess <- optimHess(c(par.control, par.sensitive, par.misreport, par.outcome), obs.llik.wrapper,
    #                      num = c(length(par.control), length(par.sensitive), length(par.misreport), length(par.outcome)),
    #                      J = J, y = y, treat = treat,
    #                      x.control = x.control, x.sensitive = x.sensitive, x.outcome = x.outcome, x.misreport = x.misreport,
    #                      d = d, sensitive.response = sensitive.response, model.misreport = model.misreport,
    #                      o = o, outcome.model = outcome.model,
    #                      outcome.constrained = outcome.constrained,
    #                      weight = weight,
    #                      respondentType = respondentType,
    #                      control.constraint = control.constraint,
    #                      control = list(reltol = 1E-16))

    num <- c(length(par.control),
             length(par.sensitive),
             length(par.misreport),
             length(par.outcome),
             length(par.outcome.aux))

    hess <- numDeriv::hessian(llik.wrapper, c(par.control, par.sensitive, par.misreport, par.outcome, par.outcome.aux),
                              num = num, J = J, y = y, w = w, treat = treat,
                              x.control = x.control, x.sensitive = x.sensitive, x.outcome = x.outcome, x.misreport = x.misreport,
                              d = d, sensitive.response = sensitive.response, model.misreport = model.misreport,
                              o = o, outcome.model = outcome.model,
                              outcome.constrained = outcome.constrained,
                              weight = weight,
                              respondentType = respondentType,
                              control.constraint = control.constraint,
                              misreport.treatment = misreport.treatment,
                              method.args = list(zero.tol = 1e-10))

    vcov.mle <- solve(-hess)
    se.mle <- sqrt(diag(vcov.mle))
    
    se.control <- se.mle[1:num[1]]
    names(se.control) <- names(par.control)
    
    se.sensitive <- se.mle[(num[1]+1):sum(num[1:2])]
    names(se.sensitive) <- names(par.sensitive)
    
    if(model.misreport == TRUE) {
      se.misreport <- se.mle[(sum(num[1:2])+1):sum(num[1:3])]
      names(se.misreport) <- names(par.misreport)
    } else {
      se.misreport <- NULL
    }
    if(outcome.model != "none") {
      se.outcome <- se.mle[(sum(num[1:3])+1):sum(num[1:4])]
      names(se.outcome) <- names(par.outcome)
      if(outcome.model %in% c("linear", "betabinomial")) {
        se.outcome.aux <- se.mle[(sum(num[1:4])+1):sum(num[1:5])]
        names(se.outcome.aux) <- names(par.outcome.aux)
      } else {
        se.outcome.aux <- NULL
      }
    } else {
     se.outcome <- NULL
     se.outcome.aux <- NULL
    }

  } else {
    se.control <- se.sensitive <- se.misreport <- se.outcome <- se.outcome.aux <- vcov.mle <- NULL
    if(se == TRUE) {
      warning("Standard errors are not implemented for models with survey weights.")
      se <- FALSE
    }
  }

  return.object <- list("par.control" = par.control,
                        "par.sensitive" = par.sensitive,
                        "par.misreport" = par.misreport,
                        "par.outcome" = par.outcome,
                        "par.outcome.aux" = par.outcome.aux,
                        "df" = n - length(c(par.control, par.sensitive, par.misreport, par.outcome, par.outcome.aux)),
                        "se.sensitive" = se.sensitive,
                        "se.control" = se.control,
                        "se.misreport" = se.misreport,
                        "se.outcome" = se.outcome,
                        "se.outcome.aux" = se.outcome.aux,
                        "vcov.mle" = vcov.mle,
                        "w" = exp(w), # Convert log posterior predicted probabilities
                        "data" = data,
                        "direct" = direct,
                        "treatment" = treatment,
                        "model.misreport" = model.misreport,
                        "outcome.model" = outcome.model,
                        "outcome.constrained" = outcome.constrained,
                        "control.constraint" = control.constraint,
                        "misreport.treatment" = misreport.treatment,
                        "weights" = weights,
                        "formula" = formula,
                        "formula.control" = formula.control,
                        "formula.sensitive" = formula.sensitive,
                        "formula.misreport" = formula.misreport,
                        "formula.outcome" = formula.outcome,
                        "sensitive.response" = sensitive.response,
                        "xlevels" = xlevels,
                        "llik" = llik,
                        "n" = n,
                        "J" = J,
                        "se" = se,
                        "runs" = runs,
                        "call" = function.call,
                        "boot" = FALSE)

  class(return.object) <- "listExperiment"
  return(return.object)

}


#' List experiment regression with bootstrapped standard errors
#' 
#' A wrapper function that makes repeated calls to \code{\link{listExperiment}}
#' to calculate parameter estimates and standard errors through non-parametric boot-strapping.
#' 
#' @param formula An object of class "\code{\link{formula}}": a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables to be used in the model.
#' @param treatment A string indicating the name of the treatment indicator in the data. This variable must be coded as a binary, where 1 indicates assignment to treatment and 0 indicates assignment to control.
#' @param J An integer indicating the number of control items in the list experiment.
#' @param direct A string indicating the name of the direct question response in the data. The direct question must be coded as a binary variable. If NULL (default), a misreport sub-model is not fit.
#' @param sensitive.response A value 0 or 1 indicating whether the response that is considered sensitive in the list experiment/direct question is 0 or 1.
#' @param outcome A string indicating the variable name in the data to use as the outcome in an outcome sub-model. If NULL (default), no outcome sub-model is fit. [\emph{experimental}]
#' @param outcome.trials An integer indicating the number of trials in a binomial/betabinomial model if both an outcome sub-model is used and if the argument \code{outcome.model} is set to "binomial" or "betabinomial". [\emph{experimental}]
#' @param outcome.model A string indicating the model type to fit for the outcome sub-model ("logistic", "binomial", "betabinomial"). [\emph{experimental}]
#' @param outcome.constrained A logical value indicating whether to constrain U* = 0 in the outcome sub-model. Defaults to TRUE. [\emph{experimental}]
#' @param control.constraint A string indicating the constraint to place on Z* and U* in the control-items sub-model:
#'    \describe{
#'         \item{"none" (default)}{Estimate separate parameters for Z* and U*.}
#'         \item{"partial"}{Constrain U* = 0.}
#'         \item{"full"}{Constrain U* = Z* = 0.}
#'     }
#' @param misreport.treatment A logical value indicating whether to include a parameter for the treatment indicator in the misreport sub-model. Defaults to TRUE.
#' @param weights A string indicating the variable name of survey weights in the data (note: standard errors are not currently output when survey weights are used).
#' @param se A logical value indicating whether to calculate standard errors. Defaults to TRUE.
#' @param tolerance The desired accuracy for EM convergence. The EM loop breaks after the change in the log-likelihood is less than the value of \code{tolerance}. Defaults to 1e-08.
#' @param max.iter The maximum number of iterations for the EM algorithm. Defaults to 10000.
#' @param n.runs The total number of times that the EM algorithm is run (can potentially help avoid local maxima). Defaults to 1.
#' @param verbose A logical value indicating whether to print information during model fitting. Defaults to TRUE.
#' @param get.data For internal use. Used by wrapper function \code{\link{bootListExperiment}}.
#' @param par.control A vector of starting parameters for the control-items sub-model. Must be in the order of the parameters in the resulting regression output. If NULL (default), randomly generated starting points are used, drawn from uniform(-2, 2).
#' @param par.sensitive A vector of starting parameters for the sensitive-item sub-model. Must be in the order of the parameters in the resulting regression output. If NULL (default), randomly generated starting points are used, drawn from uniform(-2, 2).
#' @param par.misreport A vector of starting parameters for the misreport sub-model. Must be in the order of the parameters in the resulting regression output. If NULL (default), randomly generated starting points are used, drawn from uniform(-2, 2).
#' @param par.outcome A vector of starting parameters for the outcome sub-model.  Must be in the order of the parameters in the resulting regression output. If NULL (default), randomly generated starting points are used, drawn from uniform(-2, 2). [experimental]
#' @param par.outcome.aux A vector of starting parameters for the outcome sub-model in which \code{outcome.model} is "betabinomial". i.e. c(alpha, beta). If NULL (default), randomly generated starting points are used, drawn from uniform(0, 1). [experimental]
#' @param formula.control An object of class "\code{\link{formula}}" used to specify a control-items sub-model that is different from that given in \code{formula}. (e.g. ~ x1 + x2)
#' @param formula.sensitive An object of class "\code{\link{formula}}" used to specify a sensitive-item sub-model that is different from that given in \code{formula}. (e.g. ~ x1 + x2)
#' @param formula.misreport An object of class "\code{\link{formula}}" used to specify a misreport sub-model that is different from that given in \code{formula}. (e.g. ~ x1 + x2)
#' @param formula.outcome An object of class "\code{\link{formula}}" used to specify an outcome sub-model that is different from that given in \code{formula}. (e.g. ~ x1 + x2) [\emph{experimental}]
#' @param boot.iter The number of boot strap samples to generate.
#' @param parallel A logical value indicating whether to run bootstraping in parallel on a multi-core computer.
#' @param n.cores The number of cores/threads on which to generate bootstrap samples (when \code{parallel} = TRUE). Defaults to 2.
#' @param cluster An optional cluster object using makeCluster() from the \code{parallel} package (useful if running on an MPI server).
#' 
#' @details \code{bootListExperiment} is a wrapper for the function
#'          \code{listExperiment} that allows researchers to fit a bootstrapped
#'          model. The arguments for this function include those for the
#'          \code{\link{listExperiment}} function, in addition to a small number
#'          of arguments specific to the bootstrap.
#' 
#' @return \code{listExperiment} returns an object of class "listExperiment".
#'     A summary of this object is given using the \code{\link{summary.listExperiment}}
#'     function. All components in the "listExperiment" class are listed below.
#' @slot par.control A named vector of coefficients from the control-items sub-model.
#' @slot par.sensitive A named vector of coefficients from the sensitive-item sub-model.
#' @slot par.misreport A named vector of coefficients from the misreport sub-model.
#' @slot par.outcome A named vector of coefficients from the outcome sub-model.
#' @slot par.outcome.aux A named vector of (auxiliary) coefficients from the outcome sub-model (if \code{outcome.model} = "betabinomial").
#' @slot df Degrees of freedom.
#' @slot se.sensitive Standard errors for parameters in the sensitive-item sub-model.
#' @slot se.control Standard errors for parameters in the control-items sub-model.
#' @slot se.misreport Standard errors for parameters in the misreport sub-model.
#' @slot se.outcome Standard errors for parameters in the outcome sub-model.
#' @slot se.outcome.aux Standard errors for the auxiliary parameters in the outcome sub-model (if \code{outcome.model} = "betabinomial").
#' @slot vcov.mle Variance-covariance matrix.
#' @slot w The matrix of posterior predicted probabilities for each observation in the data used for model fitting.
#' @slot data The data frame used for model fitting.
#' @slot direct The string indicating the variable name of the direct question.
#' @slot treatment The string indicating the variable name of the treatment indicator.
#' @slot model.misreport A logical value indicating whether a misreport sub-model was fit.
#' @slot outcome.model The type of model used as the outcome sub-model.
#' @slot outcome.constrained A logical value indicating whether the parameter U* was constrained to 0 in the outcome sub-model.
#' @slot control.constraint A string indicating the constraints placed on the parameters Z* and U* in the control-items sub-model.
#' @slot misreport.treatment A logical value indicating whether a treatment indicator was included in the misreport sub-model.
#' @slot weights A string indicating the variable name of the survey weights.
#' @slot formula The model formula.
#' @slot formula.control The model specification of the control-items sub-model.
#' @slot formula.sensitive The model specification of the sensitive-item sub-model.
#' @slot formula.misreport The model specification of the misreport sub-model.
#' @slot formula.outcome The model specification of the outcome sub-model.
#' @slot sensitive.response The value 0 or 1 indicating the response to the list experiment/direct question that is considered sensitive.
#' @slot xlevels The factor levels of the variables used in the model.
#' @slot llik The model log-likelihood.
#' @slot n The sample size of the data used for model fitting (this value excludes rows removed through listwise deletion).
#' @slot J The number of control items in the list experiment.
#' @slot se A logical value indicating whether standard errors were calculated.
#' @slot runs The parameter estimates from each run of the EM algorithm (note: the parameters that result in the highest log-likelihood are used as the model solution).
#' @slot call The method call.
#' @slot boot A logical value indicating whether non-parametric bootstrapping was used to calculate model parameters and standard errors.
#' 
#' @references Eady, Gregory. 2016 "The Statistical Analysis of Misreporting on Sensitive Survey Questions."
#' @references Imai, Kosuke. 2011. "Multivariate Regression Analysis for the Item Count Technique." \emph{Journal of the American Statistical Association} 106 (494): 407-416.
#' 
#' @examples
#'
#' ## Simulated list experiment and direct question
#' n <- 10000
#' J <- 4
#'
#' # Covariates
#' x <- cbind(intercept = rep(1, n), continuous1 = rnorm(n),
#'            continuous2 = rnorm(n), binary1 = rbinom(n, 1, 0.5))
#'
#' treatment <- rbinom(n, 1, 0.5)
#'
#' # Simulate Z*
#' param_sensitive <- c(0.25, -0.25, 0.5, 0.25)
#' prob_sensitive <- plogis(x %*% param_sensitive)
#' true_belief <- rbinom(n, 1, prob = prob_sensitive)
#'
#' # Simulate whether respondent misreports (U*)
#' param_misreport <- c(-0.25, 0.25, -0.5, 0.5)
#' prob_misreport <- plogis(x %*% param_misreport) * true_belief
#' misreport <- rbinom(n, 1, prob = prob_misreport)
#'
#' # Simulate control items Y*
#' param_control <- c(0.25, 0.25, -0.25, 0.25, U = -0.5, Z = 0.25)
#' prob.control <- plogis(cbind(x, misreport, true_belief) %*% param_control)
#' control_items <- rbinom(n, J, prob.control)
#'
#' # List experiment and direct question responses
#' direct <- true_belief
#' direct[misreport == 1] <- 0
#' y <- control_items + true_belief * treatment
#'
#' A <- data.frame(y, direct, treatment,
#'                 continuous1 = x[, "continuous1"],
#'                 continuous2 = x[, "continuous2"],
#'                 binary1 = x[, "binary1"])
#'
#' \dontrun{
#' # Note: substantial computation time
#' model.sim <- bootListExperiment(y ~ continuous1 + continuous2 + binary1,
#'                                 data = A, treatment = "treatment",
#'                                 direct = "direct",
#'                                 J = 4, control.constraint = "none",
#'                                 sensitive.response = 1,
#'                                 boot.iter = 500, parallel = TRUE, n.cores = 2)
#' summary(model.sim, digits = 3)
#' }
#'
#' @export
#' 
bootListExperiment <- function(formula, data, treatment, J, 
                               direct = NULL, sensitive.response = NULL,
                               outcome = NULL, outcome.trials = NULL, outcome.model = "logistic",
                               outcome.constrained = TRUE, control.constraint = "partial",
                               misreport.treatment = TRUE,
                               weights = NULL, se = TRUE, tolerance = 1E-8, max.iter = 5000,
                               n.runs = 10, verbose = TRUE, get.data = FALSE,
                               par.control = NULL, par.sensitive = NULL, par.misreport = NULL,
                               par.outcome = NULL, par.outcome.aux = NULL,
                               formula.control = NULL, formula.sensitive = NULL,
                               formula.misreport = NULL, formula.outcome = NULL,
                               boot.iter = 1000, parallel = FALSE, n.cores = 2, cluster = NULL) {

  function.call <- match.call()
  args.call <- as.list(function.call)[-1]
  args.call$se <- FALSE
  args.call$get.boot <- 1

  args.call <- lapply(args.call, eval)

  data <- args.call$data
  args.call$data <- as.name("data")
  
  # For testing:
  # return(args.call)}

  if(parallel == FALSE) {
    boot.out <- list()
    for(i in 1:boot.iter) {   
      args.call$get.boot <- i
      boot.out[[i]] <- do.call(listExperiment, args.call)
    }
  }

  if(parallel == TRUE) {
      
    args.call$verbose <- FALSE

    cat("Running bootstrap in parallel on ", n.cores, " cores/threads (", parallel::detectCores(), " available)...\n", sep = ""); Sys.sleep(0.2)

    if(!is.null(cluster)) cl <- cluster
    if(is.null(cluster)) cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl,
                  list("args.call", "data", "listExperiment", "logAdd", "estep",
                       "mstepControl", "mstepSensitive", "mstepMisreport", "mstepOutcome"),
                  envir = environment())

    boot.out <- parallel::parLapply(cl, 1:boot.iter, function(x) do.call(listExperiment, args.call))
    parallel::stopCluster(cl)
  }

  getPars <- function(varName) {
    X <- do.call(rbind, sapply(boot.out, function(x) x[varName]))
    cov.var <- cov(X)
    par.var <- colMeans(X)
    se.var <- as.vector(as.matrix(sqrt(diag(cov.var))))
    names(se.var) <- row.names(cov.var)
    return(list(par = par.var, se = se.var))
  }

  # Control items
  par.control <- getPars("par.control")$par
  se.control <- getPars("par.control")$se
  
  # Sensitive items
  par.sensitive <- getPars("par.sensitive")$par
  se.sensitive <- getPars("par.sensitive")$se
  
  # Misreport
  if(!is.null(boot.out[[1]]$par.misreport)) {
      par.misreport <- getPars("par.misreport")$par
      se.misreport <- getPars("par.misreport")$se
  } else {
    par.misreport <- se.misreport <- NULL
  }

  # Outcome
  if(!is.null(boot.out[[1]]$par.outcome)) {
    par.outcome <- getPars("par.outcome")$par
    se.outcome <- getPars("par.outcome")$se
  } else {
    par.outcome <- se.outcome <- NULL
  }

  # Outcome (auxiliary parameters)
  if(!is.null(boot.out[[1]]$outcome.model.aux)) {
    par.outcome <- getPars("par.outcome.aux")$par
    se.outcome <- getPars("par.outcome.aux")$se
  } else {
    par.outcome.aux <- se.outcome.aux <- NULL
  }

  se <- TRUE

  # Get log-likelihood and posterior probabilities with bootstrap estimates
  args.call$get.boot <- 0
  args.call$get.data <- TRUE
  args.call$par.control <- par.control
  args.call$par.sensitive <- par.sensitive
  args.call$par.misreport <- par.misreport
  args.call$par.outcome <- par.outcome
  args.call$par.outcome.aux <- par.outcome.aux
  llik <- do.call(listExperiment, args.call)$ll
  w <- do.call(listExperiment, args.call)$w

  return.object <- boot.out[[1]] # Use the first iteration object as a container
  return.object$par.control <- par.control
  return.object$par.sensitive <- par.sensitive
  return.object$par.misreport <- par.misreport
  return.object$par.outcome <- par.outcome
  return.object$par.outcome.aux <- par.outcome.aux
  return.object$df <- return.object$n - length(c(par.control, par.sensitive, par.misreport, par.outcome, par.outcome.aux))
  return.object$se.control <- se.control
  return.object$se.sensitive <- se.sensitive
  return.object$se.misreport <- se.misreport
  return.object$se.outcome <- se.outcome
  return.object$se.outcome.aux <- se.outcome.aux
  return.object$vcov.model <- NULL
  return.object$data <- data
  return.object$se <- TRUE
  
  return.object$w <- exp(w) # Convert log posterior predicted probabilities
  return.object$llik <- llik
  return.object$call <- function.call
  return.object$boot.iter <- boot.iter
  return.object$boot.out <- boot.out
  return.object$boot <- TRUE

  class(return.object) <- "listExperiment"
  return(return.object)

}

#' Predict method for the list experiment
#' 
#' Obtains predictions from a fitted list experiment model of the class \code{listExperiment}.
#' 
#' @param object Object of class "listExperiment"
#' @param newdata An optional data frame from which to calculate predictions.
#' @param treatment.misreport Value of the treatment variable covariate in the misreport sub-model (if included in the model).
#' \describe{
#'     \item{0}{treatment indicator in the misreport sub-model is set to 0 for all individuals (default).}
#'     \item{1}{treatment indicator in the misreport sub-model is set to 1 for all individuals.}
#'     \item{"observed"}{treatment indicator in the misreport sub-model is set to the observed treatment value.}
#' }
#' @param par.control An optional set of control-items sub-model parameters to use in place of those from the fitted model.
#' @param par.sensitive An optional set of sensitive-item sub-model parameters to use in place of those from the fitted model.
#' @param par.misreport An optional set of misreport sub-model parameters to use in place of those from the fitted model.
#' @param ... Additional arguments
#'
#' @details If \code{newdata} is omitted, predictions will be made with
#' the data used for model fitting.
#'
#' @slot z.hat Predicted probability of answering affirmatively to the sensitive item in the list experiment.
#' @slot u.hat Predicted probability of misreporting (assuming respondent holds the sensitive belief).
#' 
#' @references Eady, Gregory. 2016 "The Statistical Analysis of Misreporting on Sensitive Survey Questions."
#' 
#' @examples
#' 
#' data(gender)
#'
#' \dontrun{
#' # Note: substantial computation time
#' model.gender <- listExperiment(y ~ gender + ageGroup + education +
#'                                        motherTongue + region + selfPlacement,
#'                                    data = gender, J = 4,
#'                                    treatment = "treatment", direct = "direct",
#'                                    control.constraint = "none",
#'                                    sensitive.response = 0,
#'                                    misreport.treatment = TRUE)
#' predict(model.gender, treatment.misreport = 0)
#' }
#' 
#' @export
predict.listExperiment <- function(object, newdata = NULL,
                                   treatment.misreport = 0,
                                   par.control = NULL,
                                   par.sensitive = NULL,
                                   par.misreport = NULL,
                                   ...) {

  if(!is.null(par.control)) object$par.control <- par.control
  if(!is.null(par.sensitive)) object$par.sensitive <- par.sensitive
  if(!is.null(par.misreport)) object$par.misreport <- par.misreport

  if(is.null(newdata)) {
    data <- object$data
  } else data <- newdata

  if(as.character(object$formula[[2]]) %in% names(data)) {
    y <- data[, paste(object$formula[[2]])]
  } else stop(paste0("The list experiment response ", as.character(object$formula[[2]]), " not found in data."))

  if(treatment.misreport == "observed") {
    if(object$treatment %in% names(data)) {
      treatment <- data[, paste(object$treatment)]
    } else {
      stop(paste0("Argument treatment.misreport was set to \"observed\", but treatment variable \"", object$treatment, "\" is not in the data."))
    }
  } else {
    treatment <- rep(treatment.misreport, nrow(data))
  }

  if(!is.null(object$direct)) {
    if(object$direct %in% names(data)) {
      d <- data[, paste(object$direct)]
    } else {
      stop(paste0("Direct question variable", object$direct, "\" is not in the data."))
    }
  } else{
    d <- rep(NA, nrow(data))
  }

  if(!is.null(object$outcome)) {
    if(object$outcome %in% names(data)) {
      o <- data[, paste(object$outcome)]
    } else {
      stop(paste0("Outcome variable", object$outcome, "\" is not in the data."))
    }
  } else {
    o <- rep(NA, nrow(data))
  }

  if(all(all.vars(object$formula.sensitive)[-1] %in% names(data))) {
    x.sensitive <- model.matrix(object$formula.sensitive[-2], data = model.frame(~ ., data, na.action = na.pass, xlev = object$xlevels))
  } else {
    stop(paste0("Not all variables used in the sensitive-item sub-model are available in the data"))
  }

  if(!is.null(object$par.misreport)) {
    if(all(all.vars(object$formula.misreport)[-1] %in% names(data))) {
      x.misreport <- model.matrix(object$formula.misreport[-2], data = model.frame(~ ., data, na.action = na.pass, xlev = object$xlevels))
    } else {
      stop(paste0("Not all variables used in the misreport sub-model are available in the data"))
    }
  } else {
    x.misreport <- rep(NA, nrow(data))
  }

  # Prediction for Z*
  z.hat <- as.numeric(plogis(x.sensitive %*% object$par.sensitive))

  # Prediction for U*
  if(object$model.misreport == TRUE) {
    if(object$misreport.treatment == TRUE) {
      u.hat <- as.numeric(plogis(as.matrix(data.frame(x.misreport, treatment)) %*% object$par.misreport))
    } else {
      u.hat <- as.numeric(plogis(as.matrix(data.frame(x.misreport)) %*% object$par.misreport))
    }
  } else u.hat <- NULL

  return(list(z.hat = z.hat, u.hat = u.hat))

}


#' Object summary of the listExperiment class
#' 
#' Summarizes results from a list experiment regression fit using \code{\link{listExperiment}} or \code{\link{bootListExperiment}}.
#' 
#' @param object Object of class "listExperiment".
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments.
#' 
#' @details \code{summary.listExperiment} summarizes the information contained
#' in a listExperiment object for each list experiment regression sub-model.
#' 
#' @references Eady, Gregory. 2016 "The Statistical Analysis of Misreporting on Sensitive Survey Questions."
#' 
#' @examples
#' data(gender)
#'
#' \dontrun{
#' # Note: substantial computation time
#' model.gender <- listExperiment(y ~ gender + ageGroup + education +
#'                                    motherTongue + region + selfPlacement,
#'                                data = gender, J = 4,
#'                                treatment = "treatment", direct = "direct",
#'                                control.constraint = "none",
#'                                sensitive.response = 0,
#'                                misreport.treatment = TRUE)
#' summary(model.gender)
#' }
#' 
#' @export
summary.listExperiment <- function(object, digits = 4, ...) {
  cat("\nList experiment sub-models\n\n")

  cat("Call: ")
  print(object$call)

  if(object$se == TRUE) {
    cat("\nCONTROL ITEMS Pr(Y* = y)\n")
    matrix.control <- cbind(round(object$par.control, digits),
                            round(object$se.control, digits),
                            round(object$par.control/object$se.control, digits),
                            round(2 * pt(abs(object$par.control/object$se.control), object$n - object$df, lower.tail = FALSE), digits))
    colnames(matrix.control) <- c("est.", "se", "z", "p")
    print(formatC(matrix.control, format = "f", digits = digits), quote = FALSE, right = TRUE)
    cat("---\n")

    cat("\nSENSITIVE ITEM Pr(Z* = 1)\n")
    matrix.sensitive <- cbind(round(object$par.sensitive, digits),
                              round(object$se.sensitive, digits),
                              round(object$par.sensitive/object$se.sensitive, digits),
                              round(2 * pt(abs(object$par.sensitive/object$se.sensitive), object$n - object$df, lower.tail = FALSE), digits))
    colnames(matrix.sensitive) <- c("est.", "se", "z", "p")
    print(formatC(matrix.sensitive, format = "f", digits = digits), quote = FALSE, right = TRUE)
    cat("---\n")

    if(object$model.misreport == TRUE) {
      cat("\nMISREPORT Pr(U* = 1)\n")
      matrix.misreport <- cbind(round(object$par.misreport, digits),
                                   round(object$se.misreport, digits),
                                   round(object$par.misreport/object$se.misreport, digits),
                                   round(2 * pt(abs(object$par.misreport/object$se.misreport), object$n - object$df, lower.tail = FALSE), digits))
      colnames(matrix.misreport) <- c("est.", "se", "z", "p")
      print(formatC(matrix.misreport, format = "f", digits = digits), quote = FALSE, right = TRUE)
      cat("---\n")
    }

    if(object$outcome.model != "none") {
      cat("\nOUTCOME\n")
      matrix.outcome <- cbind(round(object$par.outcome, digits),
                              round(object$se.outcome, digits),
                              round(object$par.outcome/object$se.outcome, digits),
                              round(2 * pt(abs(object$par.outcome/object$se.outcome), object$n - object$df, lower.tail = FALSE), digits))
      colnames(matrix.outcome) <- c("est.", "se", "z", "p")
      print(formatC(matrix.outcome, format = "f", digits = digits), quote = FALSE, right = TRUE)
      cat("---")
    }
  } else if(object$se == FALSE) {
    cat("\nCONTROL ITEMS Pr(Y* = y)\n")
    matrix.control <- cbind(round(object$par.control, digits))
    colnames(matrix.control) <- c("est.")
    print(formatC(matrix.control, format = "f", digits = digits), quote = FALSE, right = TRUE)
    cat("---\n")

    cat("\nSENSITIVE ITEM Pr(Z* = 1)\n")
    matrix.sensitive <- cbind(round(object$par.sensitive, digits))
    colnames(matrix.sensitive) <- c("est.")
    print(formatC(matrix.sensitive, format = "f", digits = digits), quote = FALSE, right = TRUE)
    cat("---\n")

    if(object$model.misreport == TRUE) {
       cat("\nMISREPORT Pr(U* = 1)\n")
       matrix.misreport <- cbind(round(object$par.misreport, digits))
       colnames(matrix.misreport) <- c("est.")
       print(formatC(matrix.misreport, format = "f", digits = digits), quote = FALSE, right = TRUE)
       cat("---\n")
    }

    if(object$outcome.model != "none") {
       cat("\nOUTCOME\n")
       matrix.outcome <- cbind(round(object$par.outcome, digits))
       colnames(matrix.outcome) <- c("est.")
       print(formatC(matrix.outcome, format = "f", digits = digits), quote = FALSE, right = TRUE)
       cat("---")
    }
  }

  if(object$boot == TRUE) {
    cat("\nStandard errors calculated by non-parametric bootstrap (", format(object$boot.iter, big.mark = ","), " draws).", sep = "")
  }

  cat("\nObservations:", format(object$n, big.mark = ","))
  # if(nrow(object$data)-object$n != 0)
  cat(" (", format(nrow(object$data)-object$n, big.mark = ","), " of ", format(nrow(object$data), big.mark = ","), " observations removed due to missingness)", sep = "")
  cat("\nLog-likelihood", object$llik)
}

#' Print object summary of listExperiment class
#' 
#' Calls \code{\link{summary.listExperiment}}.
#' 
#' @param x Object of class "listExperiment".
#' @param ... Additional arguments.
#' 
#' @details Prints the object summary of the listExperiment class by calling the
#' \code{\link{summary.listExperiment}} function.
#' 
#' @export
print.listExperiment <- function(x, ...) {
  summary.listExperiment(x, ...)
}






# simPredict <- function(object, var, values, newdata = NULL, treatment.misreport = 0, n.sims = 1000, weight = NULL) {

#   ### Get (new)data
#   if(is.null(newdata)) {
#     data <- object$data
#   } else data <- newdata

#   if(object$treatment %in% names(data)) {
#     treat <- data[, paste(object$treatment)]
#   } else {
#     treat <- NULL
#   }

#   if(treatment.misreport == "observed") {
#     if(!is.null(treat)) {
#       treatment.predict <- treat
#     } else {
#       stop(paste0("treatment.misreport set to \"observed\", but treatment variable \"", object$treatment, "\" not found in data"))
#     }
#   } else treatment.predict <- treatment.misreport

#   if(all(all.vars(object$formula.control)[-1] %in% names(data))) {
#     x.control <- model.matrix(object$formula.control[-2], data = data, na.action = na.pass)
#   } else {
#     x.control <- matrix(NA, nrow = nrow(data), ncol = length(object$par.control))
#   }

#   if(all(all.vars(object$formula.sensitive)[-1] %in% names(data))) {
#     x.sensitive <- model.matrix(object$formula.sensitive[-2], data = data, na.action = na.pass)
#   } else{
#     x.sensitive <- matrix(NA, nrow = nrow(data), ncol = length(object$par.sensitive))
#   }

#   if(!is.null(object$par.misreport) & all(all.vars(object$formula.misreport)[-1] %in% names(data))) {
#     x.misreport <- model.matrix(object$formula.misreport[-2], data = data, na.action = na.pass)
#   } else {
#     x.misreport <- matrix(NA, nrow = nrow(data), ncol = length(object$par.misreport))
#   }

#   if(!is.null(object$par.outcome) & all(all.vars(object$formula.outcome)[-1] %in% names(data))) {
#     x.outcome <- model.matrix(object$formula.outcome[-2], data = data, na.action = na.pass)
#   } else {
#     x.outcome <- matrix(NA, nrow = nrow(data), ncol = length(object$par.outcome))
#   }


#   ### Simulate coefficients
#   coefs <- c(object$par.control, object$par.sens, object$par.misreport)
#   par_sim <- mvtnorm::rmvnorm(n.sims, coefs, object$vcov.mle)

#   # Coefficients for control-items sub-model
#   par_sim_control <- par_sim[, 1:length(object$par.control)]

#   # Coefficients for sensitive-item sub-model
#   par_sim_sensitive <- par_sim[, (length(object$par.control) + 1):(length(coefs) - length(object$par.misreport))]

#   # Coefficients for misreport sub-model
#   par_sim_misreport <- par_sim[, (length(coefs) - length(object$par.misreport) + 1):length(coefs)]

#   O <- data.frame(var = var, value = values,
#                   mean.sensitive = NA, lower.sensitive = NA, upper.sensitive = NA,
#                   mean.diff.sensitive = NA, lower.diff.sensitive = NA, upper.diff.sensitive = NA,
#                   mean.misreport = NA, lower.misreport = NA, upper.misreport = NA,
#                   mean.diff.misreport = NA, lower.diff.misreport = NA, upper.diff.misreport = NA)

#   x.sensitive[, which(colnames(x.sensitive) == var)] <- values[1]
#   x.misreport[, which(colnames(x.misreport) == var)] <- values[1]

#   x.sensitive.1 <- x.sensitive
#   x.misreport.1 <- x.misreport

#   for(i in 1:length(values)) {
#     cat(paste0("\rSimulating for ", var, " = ", values[i], "             "))

#     x.sensitive[, which(colnames(x.sensitive) == var)] <- values[i]
#     x.misreport[, which(colnames(x.misreport) == var)] <- values[i]
#     out_sensitive <- apply(par_sim_sensitive, 1, function(x) mean(plogis(x.sensitive %*% x)))
#     out_sensitive.diff <- apply(par_sim_sensitive, 1, function(x) mean(plogis(x.sensitive %*% x) - plogis(x.sensitive.1 %*% x)))
#     out_misreport <- apply(par_sim_misreport, 1, function(x) mean(plogis(x.misreport %*% x)))
#     out_misreport.diff <- apply(par_sim_misreport, 1, function(x) mean(plogis(x.misreport %*% x) - plogis(x.misreport.1 %*% x)))

#     O$mean.sensitive[i] <- mean(out_sensitive)
#     O$lower.sensitive[i] <- quantile(out_sensitive, 0.05)
#     O$upper.sensitive[i] <- quantile(out_sensitive, 0.95)
#     O$mean.diff.sensitive[i] <- mean(out_sensitive.diff)
#     O$lower.diff.sensitive[i] <- quantile(out_sensitive.diff, 0.05)
#     O$upper.diff.sensitive[i] <- quantile(out_sensitive.diff, 0.95)

#     O$mean.misreport[i] <- mean(out_misreport)
#     O$lower.misreport[i] <- quantile(out_misreport, 0.05)
#     O$upper.misreport[i] <- quantile(out_misreport, 0.95)
#     O$mean.diff.misreport[i] <- mean(out_sensitive.diff)
#     O$lower.diff.misreport[i] <- quantile(out_misreport.diff, 0.05)
#     O$upper.diff.misreport[i] <- quantile(out_misreport.diff, 0.95)

#   }

# ggplot(O, aes(x = value, y = mean.sensitive,
#               ymin = lower.sensitive, ymax = upper.sensitive)) +
#   my.theme() +
#   geom_ribbon(fill = "grey94", color = "grey90", size = 0.25) +
#   geom_line()

# }