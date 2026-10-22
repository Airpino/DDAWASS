#' @useDynLib DDAWASS, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import methods
NULL

## ---------------------------------------------------------------------------
## distributionD: a discrete distribution with finite support
## ---------------------------------------------------------------------------

#' Class distributionD
#'
#' Class \code{distributionD} describes a discrete distribution with a finite
#' support (e.g. the frequency distribution of a Likert/0-10 rating, of a
#' count, or of any variable observed on a finite set of numeric values).
#' It is the discrete counterpart of the \code{distributionH} class of
#' \pkg{HistDAWass}.
#'
#' @slot x numeric vector: the support values, sorted in increasing order and
#'   without duplicates.
#' @slot p numeric vector: the probability masses (relative frequencies)
#'   associated with \code{x}; they are positive and sum to one.
#' @slot cdf numeric vector: the cumulative distribution function evaluated
#'   at \code{x} (\code{cumsum(p)}, last value equal to 1).
#' @slot n numeric: the number of observations the distribution was computed
#'   from (\code{NA} when unknown).
#' @slot m numeric: the mean of the distribution.
#' @slot s numeric: the standard deviation of the distribution.
#' @seealso \code{\link{distributionD}} for the constructor.
#' @name distributionD-class
#' @rdname distributionD-class
#' @exportClass distributionD
setClass(
  Class = "distributionD",
  representation = representation(
    x = "numeric", p = "numeric", cdf = "numeric",
    n = "numeric", m = "numeric", s = "numeric"
  ),
  prototype = prototype(
    x = numeric(0), p = numeric(0), cdf = numeric(0),
    n = NA_real_, m = NA_real_, s = NA_real_
  ),
  validity = function(object) {
    if (length(object@x) == 0) return(TRUE) # empty distribution allowed
    if (length(object@x) != length(object@p) ||
        length(object@x) != length(object@cdf)) {
      return("x, p and cdf must have the same length")
    }
    if (is.unsorted(object@x, strictly = TRUE)) {
      return("x must be strictly increasing")
    }
    if (any(object@p < 0)) return("p must be non-negative")
    if (abs(sum(object@p) - 1) > 1e-8) return("p must sum to 1")
    if (abs(object@cdf[length(object@cdf)] - 1) > 1e-8) {
      return("the last value of cdf must be 1")
    }
    TRUE
  }
)

#' Constructor of a discrete distribution
#'
#' Builds a \code{\linkS4class{distributionD}} object. Three ways to call it:
#' \itemize{
#'   \item \code{distributionD(x, p)}: support values and probabilities;
#'   \item \code{distributionD(x, counts = k)}: support values and absolute
#'     frequencies;
#'   \item \code{distributionD(obs)}: a vector of raw observations, which is
#'     tabulated.
#' }
#' Duplicated support values are merged, values with zero mass are dropped
#' and the probabilities are normalised to sum to one.
#'
#' @param x numeric vector of support values, or of raw observations when
#'   both \code{p} and \code{counts} are \code{NULL}.
#' @param p (optional) numeric vector of probabilities/relative frequencies.
#' @param counts (optional) numeric vector of absolute frequencies.
#' @param keep.zeros logical; keep support values with zero mass? Default
#'   \code{FALSE}.
#' @return A \code{distributionD} object.
#' @examples
#' # from support and probabilities
#' d1 <- distributionD(x = 1:5, p = c(0.1, 0.2, 0.4, 0.2, 0.1))
#' # from counts
#' d2 <- distributionD(x = 0:10, counts = c(1, 0, 2, 3, 5, 10, 12, 20, 25, 15, 7))
#' # from raw observations
#' d3 <- distributionD(c(3, 4, 4, 5, 5, 5, 2, 4))
#' d1; meanD(d2); stdD(d3)
#' @export
distributionD <- function(x = numeric(0), p = NULL, counts = NULL,
                          keep.zeros = FALSE) {
  if (length(x) == 0) return(new("distributionD"))
  x <- as.numeric(x)
  n <- NA_real_
  if (is.null(p) && is.null(counts)) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(new("distributionD"))
    tab <- table(x)
    counts <- as.numeric(tab)
    x <- as.numeric(names(tab))
  }
  if (!is.null(counts)) {
    if (length(counts) != length(x)) stop("x and counts must have the same length")
    if (any(counts < 0, na.rm = TRUE)) stop("counts must be non-negative")
    counts[is.na(counts)] <- 0
    n <- sum(counts)
    if (n <= 0) stop("counts sum to zero")
    p <- counts / n
  } else {
    if (length(p) != length(x)) stop("x and p must have the same length")
    if (any(p < 0, na.rm = TRUE)) stop("p must be non-negative")
    p[is.na(p)] <- 0
    if (sum(p) <= 0) stop("p sums to zero")
    p <- p / sum(p)
  }
  ok <- !is.na(x)
  x <- x[ok]
  p <- p[ok]
  # merge duplicates and sort
  if (anyDuplicated(x)) {
    agg <- tapply(p, x, sum)
    x <- as.numeric(names(agg))
    p <- as.numeric(agg)
  }
  o <- order(x)
  x <- x[o]
  p <- p[o]
  if (!keep.zeros) {
    keep <- p > 0
    x <- x[keep]
    p <- p[keep]
  }
  .newD(x, p, n)
}

# internal fast constructor (x sorted, p normalised)
.newD <- function(x, p, n = NA_real_) {
  p <- p / sum(p)
  cdf <- cumsum(p)
  cdf[length(cdf)] <- 1
  m <- sum(x * p)
  v <- sum((x - m)^2 * p)
  new("distributionD", x = x, p = p, cdf = cdf, n = as.numeric(n),
      m = m, s = sqrt(max(v, 0)))
}

# compact a distribution whose support may have (nearly) equal consecutive values
.compactD <- function(x, cdf, tol = 1e-10) {
  p <- diff(c(0, cdf))
  keep <- p > tol
  x <- x[keep]
  p <- p[keep]
  if (length(x) > 1) {
    grp <- cumsum(c(TRUE, diff(x) > tol * max(1, abs(x))))
    if (any(duplicated(grp))) {
      x <- as.numeric(tapply(x, grp, function(v) v[1]))
      p <- as.numeric(tapply(p, grp, sum))
    }
  }
  .newD(x, p)
}

## ---------------------------------------------------------------------------
## MatD: a matrix of discrete distributions
## ---------------------------------------------------------------------------

#' Class MatD
#'
#' Class \code{MatD} is a matrix whose cells are
#' \code{\linkS4class{distributionD}} objects: rows are statistical units
#' (individuals, groups, territories, ...) and columns are discrete
#' distributional variables. It mirrors the \code{MatH} class of
#' \pkg{HistDAWass}.
#'
#' @slot M a matrix of lists, each cell containing a \code{distributionD}.
#' @name MatD-class
#' @rdname MatD-class
#' @exportClass MatD
setClass(Class = "MatD", representation = representation(M = "matrix"))

#' Constructor of a matrix of discrete distributions
#'
#' @param x a list of \code{distributionD} objects.
#' @param nrows number of rows (units).
#' @param ncols number of columns (variables).
#' @param rownames (optional) names of the rows.
#' @param varnames (optional) names of the columns.
#' @param by.row logical; if \code{TRUE} the list \code{x} fills the matrix by
#'   rows, otherwise (default) by columns.
#' @return A \code{MatD} object.
#' @examples
#' d <- list(distributionD(c(1, 2, 2, 3)), distributionD(c(2, 3, 3, 3)),
#'           distributionD(c(1, 1, 2)),    distributionD(c(3, 3, 2)))
#' MatD(d, nrows = 2, ncols = 2, rownames = c("A", "B"), varnames = c("Q1", "Q2"))
#' @export
MatD <- function(x = list(), nrows = 1, ncols = 1, rownames = NULL,
                 varnames = NULL, by.row = FALSE) {
  M <- matrix(list(new("distributionD")), nrows, ncols)
  if (length(x) > 0) {
    if (length(x) != nrows * ncols) {
      stop("length(x) must be equal to nrows * ncols")
    }
    idx <- matrix(seq_len(nrows * ncols), nrows, ncols, byrow = by.row)
    for (i in seq_len(nrows)) for (j in seq_len(ncols)) {
      el <- x[[idx[i, j]]]
      if (!is(el, "distributionD")) stop("all the elements of x must be distributionD objects")
      M[i, j][[1]] <- el
    }
  }
  rownames(M) <- if (is.null(rownames)) paste0("I", seq_len(nrows)) else as.character(rownames)
  colnames(M) <- if (is.null(varnames)) paste0("X", seq_len(ncols)) else as.character(varnames)
  new("MatD", M = M)
}
