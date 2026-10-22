## Methods for single discrete distributions -----------------------------

#' @importFrom HistDAWass compQ compP rQQ
NULL

setMethod("show", "distributionD", function(object) {
  if (length(object@x) == 0) {
    cat("An empty discrete distribution\n")
    return(invisible(NULL))
  }
  cat("Discrete distribution (", length(object@x), " support values",
      if (!is.na(object@n)) paste0(", n = ", object@n), ")\n", sep = "")
  df <- data.frame(x = object@x, p = round(object@p, 5),
                   cdf = round(object@cdf, 5))
  if (nrow(df) > 12) {
    print(utils::head(df, 6), row.names = FALSE)
    cat("   ... (", nrow(df) - 12, " rows omitted) ...\n", sep = "")
    print(utils::tail(df, 6), row.names = FALSE)
  } else {
    print(df, row.names = FALSE)
  }
  cat(sprintf("mean = %g   std = %g\n", object@m, object@s))
  invisible(NULL)
})

#' Basic statistics of a discrete distribution
#'
#' \code{meanD}, \code{stdD}, \code{varD}, \code{skewD} and \code{kurtD}
#' return the mean, the standard deviation, the variance, the (Pearson)
#' skewness and the (non-excess) kurtosis of a
#' \code{\linkS4class{distributionD}}.
#' @param object a \code{distributionD} object.
#' @return a number.
#' @examples
#' d <- distributionD(x = 1:5, p = c(0.05, 0.1, 0.2, 0.4, 0.25))
#' meanD(d); stdD(d); skewD(d); kurtD(d)
#' @name statsD
#' @aliases meanD
#' @rdname statsD
#' @export
setGeneric("meanD", function(object) standardGeneric("meanD"))
#' @rdname statsD
#' @export
setGeneric("stdD", function(object) standardGeneric("stdD"))
#' @rdname statsD
#' @export
setGeneric("varD", function(object) standardGeneric("varD"))
#' @rdname statsD
#' @export
setGeneric("skewD", function(object) standardGeneric("skewD"))
#' @rdname statsD
#' @export
setGeneric("kurtD", function(object) standardGeneric("kurtD"))

#' @rdname statsD
setMethod("meanD", "distributionD", function(object) object@m)
#' @rdname statsD
setMethod("stdD", "distributionD", function(object) object@s)
#' @rdname statsD
setMethod("varD", "distributionD", function(object) object@s^2)
#' @rdname statsD
setMethod("skewD", "distributionD", function(object) {
  if (object@s == 0) return(0)
  sum(((object@x - object@m) / object@s)^3 * object@p)
})
#' @rdname statsD
setMethod("kurtD", "distributionD", function(object) {
  if (object@s == 0) return(NA_real_)
  sum(((object@x - object@m) / object@s)^4 * object@p)
})

#' Quantiles and cdf of a discrete distribution
#'
#' Methods of the \pkg{HistDAWass} generics \code{compQ} and \code{compP}
#' for \code{\linkS4class{distributionD}} objects. The quantile function is
#' the left-continuous inverse of the cdf, \eqn{Q(t) = \min\{x: F(x) \ge t\}},
#' with \eqn{Q(0) = \min(x)}.
#' @param object a \code{distributionD}.
#' @param p a vector of probability levels in \eqn{[0, 1]}.
#' @param q a vector of values.
#' @return a numeric vector.
#' @examples
#' d <- distributionD(x = 1:5, p = c(0.05, 0.1, 0.2, 0.4, 0.25))
#' compQ(d, c(0, 0.25, 0.5, 0.75, 1))
#' compP(d, c(2, 3.5))
#' @name compQ-distributionD
#' @aliases compQ,distributionD-method compP,distributionD-method
#'   compQ,distributionD,ANY-method compP,distributionD,ANY-method
#' @exportMethod compQ compP
setMethod("compQ", "distributionD", function(object, p) {
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("p must be in [0, 1]")
  .wd_quantile_c(object@x, object@cdf, as.numeric(p))
})

#' @rdname compQ-distributionD
setMethod("compP", "distributionD", function(object, q) {
  vapply(q, function(v) {
    k <- findInterval(v, object@x)
    if (k == 0) 0 else object@cdf[k]
  }, numeric(1))
})

#' L2 Wasserstein distance between discrete distributions
#'
#' Computes the (squared) L2 Wasserstein distance between two
#' \code{\linkS4class{distributionD}} objects, i.e. the L2 distance between
#' their quantile functions. The computation is exact: quantile functions of
#' discrete distributions are step functions and the integral is evaluated on
#' the merged grid of cumulative probabilities.
#'
#' With \code{details = TRUE} the squared distance is decomposed (Irpino and
#' Romano, 2007) into a location, a size and a shape component:
#' \deqn{d_W^2 = (\mu_1 - \mu_2)^2 + (\sigma_1 - \sigma_2)^2 + 2\sigma_1\sigma_2(1 - \rho_{QQ})}
#' @param object1,object2 two \code{distributionD} objects.
#' @param details logical; return the decomposition?
#' @param ... unused.
#' @return the squared distance or, if \code{details = TRUE}, a named vector
#'   with the components \code{SQ.W}, \code{POSITION}, \code{SIZE},
#'   \code{SHAPE} and \code{rQQ}.
#' @references Irpino A., Romano E. (2007). Optimal histogram representation
#'   of large data sets: Fisher vs piecewise linear approximation.
#'   \emph{RNTI} E-9, 99--110.
#' @examples
#' a <- distributionD(x = 1:5, p = c(0.05, 0.1, 0.2, 0.4, 0.25))
#' b <- distributionD(x = 1:5, p = c(0.3, 0.3, 0.2, 0.1, 0.1))
#' WassSqDistD(a, b)
#' WassSqDistD(a, b, details = TRUE)
#' @export
setGeneric("WassSqDistD", function(object1, object2, ...) standardGeneric("WassSqDistD"))

#' @rdname WassSqDistD
setMethod("WassSqDistD", signature(object1 = "distributionD", object2 = "distributionD"),
  function(object1, object2, details = FALSE, ...) {
    d2 <- .wd_sq_c(object1@x, object1@cdf, object2@x, object2@cdf)
    if (!details) return(d2)
    pos <- (object1@m - object2@m)^2
    siz <- (object1@s - object2@s)^2
    sha <- max(d2 - pos - siz, 0)
    rqq <- if (object1@s > 0 && object2@s > 0) 1 - sha / (2 * object1@s * object2@s) else NA_real_
    c(SQ.W = d2, POSITION = pos, SIZE = siz, SHAPE = sha, rQQ = rqq)
  })

#' Quantile-quantile correlation of two discrete distributions
#'
#' Method of the \pkg{HistDAWass} generic \code{rQQ}: the correlation between
#' the quantile functions of two \code{\linkS4class{distributionD}}.
#' @param e1,e2 two \code{distributionD} objects.
#' @return a number in \eqn{[-1, 1]}.
#' @examples
#' a <- distributionD(x = 1:5, p = c(0.05, 0.1, 0.2, 0.4, 0.25))
#' b <- distributionD(x = 1:5, p = c(0.3, 0.3, 0.2, 0.1, 0.1))
#' rQQ(a, b)
#' @name rQQ-distributionD
#' @aliases rQQ,distributionD,distributionD-method
#' @exportMethod rQQ
setMethod("rQQ", signature(e1 = "distributionD", e2 = "distributionD"), function(e1, e2) {
  if (e1@s == 0 || e2@s == 0) return(NA_real_)
  (.wd_dot_c(e1@x, e1@cdf, e2@x, e2@cdf) - e1@m * e2@m) / (e1@s * e2@s)
})

#' Inner product of quantile functions
#'
#' \code{dotpD(a, b)} computes \eqn{\int_0^1 Q_a(t) Q_b(t) dt}, the
#' building block of Wasserstein-based covariances.
#' @param a,b two \code{distributionD} objects.
#' @return a number.
#' @examples
#' a <- distributionD(x = 1:3, p = c(0.2, 0.5, 0.3))
#' b <- distributionD(x = c(2, 4), p = c(0.5, 0.5))
#' dotpD(a, b)
#' @export
dotpD <- function(a, b) .wd_dot_c(a@x, a@cdf, b@x, b@cdf)

## Arithmetic in the Wasserstein (quantile) space ------------------------

#' Arithmetic on discrete distributions
#'
#' Operations are defined on quantile functions, consistently with the
#' L2 Wasserstein geometry: the sum of two distributions is the distribution
#' whose quantile function is \eqn{Q_1 + Q_2}; a distribution plus a
#' number is a shift; a distribution times a number is a rescaling (a
#' negative number also reverses the support, \eqn{Q(t) \to aQ(1-t)}).
#' @param e1,e2 \code{distributionD} objects or numbers.
#' @return a \code{distributionD}.
#' @examples
#' a <- distributionD(x = 1:3, p = c(0.2, 0.5, 0.3))
#' b <- distributionD(x = c(10, 20), p = c(0.5, 0.5))
#' a + b
#' a * 2 + 1
#' @name Arith-distributionD
#' @aliases +,distributionD,distributionD-method +,distributionD,numeric-method
#'   +,numeric,distributionD-method *,distributionD,numeric-method
#'   *,numeric,distributionD-method -,distributionD,numeric-method
#' @exportMethod + - *
setMethod("+", signature(e1 = "distributionD", e2 = "distributionD"), function(e1, e2) {
  grid <- sort(unique(c(e1@cdf, e2@cdf)))
  grid[length(grid)] <- 1
  x <- compQ(e1, grid) + compQ(e2, grid)
  .compactD(x, grid)
})
setMethod("+", signature(e1 = "distributionD", e2 = "numeric"), function(e1, e2) {
  .newD(e1@x + e2, e1@p, e1@n)
})
setMethod("+", signature(e1 = "numeric", e2 = "distributionD"), function(e1, e2) e2 + e1)
setMethod("-", signature(e1 = "distributionD", e2 = "numeric"), function(e1, e2) e1 + (-e2))
setMethod("*", signature(e1 = "distributionD", e2 = "numeric"), function(e1, e2) {
  if (e2 == 0) return(.newD(0, 1, e1@n))
  if (e2 > 0) return(.newD(e1@x * e2, e1@p, e1@n))
  .newD(rev(e1@x * e2), rev(e1@p), e1@n)
})
setMethod("*", signature(e1 = "numeric", e2 = "distributionD"), function(e1, e2) e2 * e1)

#' Plot a discrete distribution
#'
#' @param x a \code{distributionD}.
#' @param y unused.
#' @param type \code{"BAR"} (default, probability masses) or \code{"CDF"}.
#' @param col fill/line colour.
#' @param ... unused.
#' @return a \code{ggplot} object.
#' @examples
#' d <- distributionD(x = 0:10, counts = c(1, 0, 2, 3, 5, 10, 12, 20, 25, 15, 7))
#' plot(d)
#' plot(d, type = "CDF")
#' @aliases plot,distributionD-method
#' @exportMethod plot
setMethod("plot", signature(x = "distributionD"), function(x, y, type = "BAR",
                                                           col = "#2E7D32", ...) {
  df <- data.frame(x = x@x, p = x@p, cdf = x@cdf)
  if (toupper(type) == "CDF") {
    g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$cdf)) +
      ggplot2::geom_step(direction = "hv", colour = col, linewidth = 0.8) +
      ggplot2::geom_point(colour = col) +
      ggplot2::labs(x = NULL, y = "F(x)")
  } else {
    g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$p)) +
      ggplot2::geom_col(fill = col, width = 0.8 * .min_gap(x@x)) +
      ggplot2::labs(x = NULL, y = "p(x)")
  }
  g + ggplot2::theme_minimal()
})

.min_gap <- function(x) if (length(x) < 2) 1 else min(diff(x))

#' @importFrom rlang .data
NULL
