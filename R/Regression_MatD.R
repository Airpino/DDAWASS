## Two-component linear regression for discrete distributional data ------

#' Two-component regression of discrete distributional variables
#'
#' Fits the two-component linear model of Irpino and Verde (2015) to discrete
#' distributional data:
#' \deqn{Q_{y_i}(t) = \beta_0 + \sum_j \beta_j \mu_{x_{ij}} + \sum_j \gamma_j Q^c_{x_{ij}}(t) + e_i(t)}
#' where \eqn{Q^c} is the centred quantile function (\eqn{Q - \mu}) and
#' \eqn{\gamma_j \ge 0} so that the predicted quantile functions are
#' non-decreasing. Since centred quantile functions are orthogonal to
#' constants, the least-squares problem splits into an OLS regression on the
#' means (\eqn{\beta}) and a non-negative least-squares problem on the centred
#' parts (\eqn{\gamma}), both solved exactly.
#' @param data a \code{MatD}.
#' @param Yvar the response column.
#' @param Xvars the explanatory columns.
#' @return a named vector of parameters: \code{(AV_Intercept)}, the
#'   \code{AV_} coefficients of the means and the \code{CEN_} coefficients of
#'   the centred distributions.
#' @references Irpino A., Verde R. (2015). Linear regression for numeric
#'   symbolic variables: a least squares approach based on Wasserstein
#'   distance. \emph{Advances in Data Analysis and Classification}, 9(1),
#'   81--106. \doi{10.1007/s11634-015-0197-7}
#' @examples
#' pars <- WD.regression.two.components(Satisfaction, Yvar = "Overall",
#'                                      Xvars = c("Staff", "Price", "Cleanliness"))
#' pars
#' xv <- c("Staff", "Price", "Cleanliness")
#' pred <- WD.regression.two.components.predict(Satisfaction[, xv], pars)
#' WD.regression.GOF(Satisfaction[, "Overall"], pred)[c("RMSE_W", "OMEGA")]
#' @export
WD.regression.two.components <- function(data, Yvar, Xvars) {
  Y <- data[, Yvar]
  X <- data[, Xvars]
  .check_nonempty(Y); .check_nonempty(X)
  n <- nrow(data@M); d <- ncol(X@M)
  xn <- colnames(X@M); yn <- colnames(Y@M)
  mY <- get.MatD.stats(Y, "mean")$mat[, 1]
  mX <- get.MatD.stats(X, "mean")$mat
  # OLS on means
  fit <- stats::lm.fit(cbind(1, mX), mY)
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  # Gram of centred quantile functions
  A <- matrix(0, d, d)
  b <- numeric(d)
  yc <- .col_xc(Y, 1)
  xcs <- lapply(seq_len(d), function(j) .col_xc(X, j))
  for (j in seq_len(d)) {
    for (i in seq_len(n)) {
      b[j] <- b[j] + .wd_dot_c(xcs[[j]]$x[[i]], xcs[[j]]$c[[i]], yc$x[[i]], yc$c[[i]]) -
        mX[i, j] * mY[i]
    }
    for (k in j:d) {
      s <- 0
      for (i in seq_len(n)) {
        s <- s + .wd_dot_c(xcs[[j]]$x[[i]], xcs[[j]]$c[[i]], xcs[[k]]$x[[i]], xcs[[k]]$c[[i]]) -
          mX[i, j] * mX[i, k]
      }
      A[j, k] <- A[k, j] <- s
    }
  }
  gam <- .nnls_gram(A, b)
  pars <- c(beta, gam)
  names(pars) <- c("(AV_Intercept)", paste0("AV_", xn), paste0("CEN_", xn))
  pars
}

# Lawson-Hanson active set for min g'Ag - 2b'g, g >= 0
.nnls_gram <- function(A, b, tol = 1e-12, maxit = 500) {
  d <- length(b)
  A <- A + diag(1e-12 * max(1, max(abs(diag(A)))), d)
  x <- numeric(d)
  P <- rep(FALSE, d)
  wgrad <- b - A %*% x
  it <- 0
  while (any(!P & wgrad > tol) && it < maxit) {
    it <- it + 1
    cand <- which(!P)
    j <- cand[which.max(wgrad[cand])]
    P[j] <- TRUE
    repeat {
      z <- numeric(d)
      z[P] <- solve(A[P, P, drop = FALSE], b[P])
      if (all(z[P] > tol)) {
        x <- z
        break
      }
      neg <- which(P & z <= tol)
      alpha <- min(x[neg] / (x[neg] - z[neg]))
      x <- x + alpha * (z - x)
      P <- P & x > tol
      x[!P] <- 0
    }
    wgrad <- b - A %*% x
  }
  as.numeric(x)
}

#' Predictions of the two-component regression
#'
#' @param data a \code{MatD} with the explanatory variables (in the same
#'   order used for the estimation).
#' @param parameters the output of \code{\link{WD.regression.two.components}}.
#' @return a one-column \code{MatD} of predicted distributions.
#' @examples
#' xv <- c("Staff", "Price")
#' pars <- WD.regression.two.components(Satisfaction, "Overall", xv)
#' pred <- WD.regression.two.components.predict(Satisfaction[, xv], pars)
#' pred[1:3, ]
#' @export
WD.regression.two.components.predict <- function(data, parameters) {
  .check_nonempty(data)
  n <- nrow(data@M); d <- ncol(data@M)
  if (length(parameters) != 2 * d + 1) stop("parameters do not match the number of columns")
  beta <- parameters[1:(d + 1)]
  gam <- parameters[(d + 2):(2 * d + 1)]
  out <- vector("list", n)
  for (i in seq_len(n)) {
    cells <- data@M[i, ]
    grid <- sort(unique(unlist(lapply(cells, function(dd) dd@cdf))))
    grid[length(grid)] <- 1
    q <- rep(beta[1], length(grid))
    for (j in seq_len(d)) {
      dd <- cells[[j]]
      q <- q + beta[j + 1] * dd@m + gam[j] * (compQ(dd, grid) - dd@m)
    }
    out[[i]] <- .compactD(q, grid)
  }
  MatD(out, nrows = n, ncols = 1, rownames = rownames(data@M), varnames = "Predicted")
}

#' Goodness of fit indices for distributional regression
#'
#' Computes the root mean squared Wasserstein error (\code{RMSE_W}), the
#' \code{OMEGA} index and the pseudo-\eqn{R^2} of Irpino and Verde (2015),
#' with the same definitions of \code{HistDAWass::WH.regression.GOF}.
#' @param observed a one-column \code{MatD} of observed distributions.
#' @param predicted a one-column \code{MatD} of predicted distributions.
#' @return a list with \code{RMSE_W}, \code{OMEGA} and \code{PSEUDOR2}.
#' @examples
#' xv <- c("Staff", "Cleanliness")
#' pars <- WD.regression.two.components(Satisfaction, "Overall", xv)
#' pred <- WD.regression.two.components.predict(Satisfaction[, xv], pars)
#' gof <- WD.regression.GOF(Satisfaction[, "Overall"], pred)
#' gof$OMEGA
#' @export
WD.regression.GOF <- function(observed, predicted) {
  n <- nrow(observed@M)
  MO <- WD.vec.mean(observed)
  WassD2 <- vapply(seq_len(n), function(i)
    WassSqDistD(observed@M[i, 1][[1]], predicted@M[i, 1][[1]]), 0)
  TOTSSQ <- sum(vapply(seq_len(n), function(i) WassSqDistD(observed@M[i, 1][[1]], MO), 0))
  SSQR <- sum(vapply(seq_len(n), function(i) WassSqDistD(predicted@M[i, 1][[1]], MO), 0))
  num <- sum(vapply(seq_len(n), function(i) {
    pr <- predicted@M[i, 1][[1]]; (pr@m - MO@m)^2 + pr@s^2 }, 0))
  den <- sum(vapply(seq_len(n), function(i) {
    ob <- observed@M[i, 1][[1]]; (ob@m - MO@m)^2 + ob@s^2 }, 0))
  SSE <- sum(WassD2)
  list(RMSE_W = sqrt(SSE / n), OMEGA = num / den,
       PSEUDOR2 = list(index = max(0, min(SSQR / TOTSSQ, 1 - SSE / TOTSSQ)),
                       details = c(TotSSQ = TOTSSQ, SSQ.R = SSQR, SSQ.E = SSE,
                                   Bias = TOTSSQ - SSQR - SSE,
                                   SSQ.R.rel = SSQR / TOTSSQ, SSQ.E.rel = SSE / TOTSSQ,
                                   SSQ.bias.rel = (TOTSSQ - SSQR - SSE) / TOTSSQ)))
}
