## Projection of discrete distributions on a given support ---------------

#' Project discrete distributions on a fixed support
#'
#' \code{simplifyD} moves the mass of each support value of a distribution to
#' the nearest value of a given \code{support} (e.g. the categories
#' \code{1:5} of a Likert scale): the value \eqn{s_k} receives the masses of
#' all the points in \eqn{((s_{k-1} + s_k)/2, (s_k + s_{k+1})/2]}; masses
#' below the first midpoint go to \eqn{s_1} and masses above the last one to
#' \eqn{s_K}. Ties (points exactly at a midpoint) go to the lower value.
#'
#' The nearest-point rule gives the distribution supported on \code{support}
#' closest to the original one in the L2 Wasserstein distance (indeed, in any
#' Wasserstein distance): it is the Wasserstein projection on the set of
#' distributions supported on \code{support}. It is typically used to express
#' a barycenter (whose support is generally made of non-integer values) on the
#' original categories of the variable, see \code{\link{WD.vec.mean.simplified}}.
#' The mean of the projected distribution may differ from the original one,
#' by at most half of the largest gap between consecutive support values.
#'
#' For a barycenter \eqn{\bar\mu} of \eqn{\mu_1, \dots, \mu_n}, since
#' \eqn{\sum_i w_i d_W^2(\mu_i, \nu) = \sum_i w_i d_W^2(\mu_i, \bar\mu) +
#' d_W^2(\bar\mu, \nu)}, the projection of \eqn{\bar\mu} on the support is
#' also a barycenter of the data constrained to have that support. A proof
#' is given in the file \code{barycenter_projection.pdf} installed with the
#' package: \code{system.file("doc-theory", "barycenter_projection.pdf",
#' package = "DDAWASS")}.
#'
#' @param object a \code{distributionD} or a \code{MatD}.
#' @param support the target support values. For a \code{MatD} it can be a
#'   list with one vector per column; if \code{NULL} (default, \code{MatD}
#'   only) each column is projected on the union of the support values
#'   observed in that column.
#' @param keep.zeros logical; keep the support values that receive zero mass?
#' @return an object of the same class as \code{object}.
#' @examples
#' b <- distributionD(x = c(1, 1.1, 1.35, 1.6, 1.8, 2.5, 2.6, 3.4, 4.9),
#'                    p = c(1, 2, 4, 8, 3, 5, 10, 6, 1))
#' s <- simplifyD(b, support = 1:5)
#' s
#' WassSqDistD(b, s)   # the smallest possible among distributions on 1:5
#' @export
simplifyD <- function(object, support = NULL, keep.zeros = TRUE) {
  if (is(object, "distributionD")) {
    if (is.null(support)) stop("support must be given for a distributionD")
    return(.project_support(object, support, keep.zeros))
  }
  if (!is(object, "MatD")) stop("object must be a distributionD or a MatD")
  p <- ncol(object@M)
  for (j in seq_len(p)) {
    sup <- if (is.list(support)) support[[j]] else if (!is.null(support)) support else
      sort(unique(unlist(lapply(object@M[, j], function(d) d@x))))
    for (i in seq_len(nrow(object@M))) {
      d <- object@M[i, j][[1]]
      if (length(d@x) > 0) object@M[i, j][[1]] <- .project_support(d, sup, keep.zeros)
    }
  }
  object
}

.project_support <- function(d, support, keep.zeros) {
  s <- sort(unique(as.numeric(support)))
  if (length(s) == 1) return(.newD(s, 1, d@n))
  mids <- (s[-1] + s[-length(s)]) / 2
  k <- findInterval(d@x, mids, left.open = TRUE) + 1
  p <- vapply(seq_along(s), function(h) sum(d@p[k == h]), 0)
  out <- distributionD(s, p = p, keep.zeros = keep.zeros)
  out@n <- d@n
  out
}

#' Barycenter expressed on the original support values
#'
#' Computes the exact Wasserstein barycenter of each column of a \code{MatD}
#' (see \code{\link{WD.vec.mean}}) and then projects it on the support of the
#' variable with \code{\link{simplifyD}}, so that the result has the same
#' categories of the data (e.g. the 1-7 scores of a Likert item) and new
#' probabilities. The result is the distribution on the given support closest
#' to the exact barycenter, and it is also the barycenter of the data
#' constrained to that support (see \code{\link{simplifyD}}). It is an
#' approximation, useful for presentation and
#' interpretation: the exact barycenter, returned by \code{WD.vec.mean}, is
#' the one used in all the computations of the package (variances,
#' covariances, k-means prototypes, ...).
#' @param object a \code{MatD}.
#' @param w (optional) weights of the rows.
#' @param support (optional) the target support (a vector, or a list with
#'   one vector per column); by default the union of the support values
#'   observed in each column.
#' @return a \code{distributionD} if \code{object} has one column, otherwise
#'   a one-row \code{MatD}.
#' @examples
#' exact <- WD.vec.mean(Satisfaction[, "Staff"])
#' approx <- WD.vec.mean.simplified(Satisfaction[, "Staff"], support = 0:10)
#' approx
#' WassSqDistD(exact, approx)
#' @export
WD.vec.mean.simplified <- function(object, w = 1, support = NULL) {
  p <- ncol(object@M)
  if (is.null(support)) {
    support <- lapply(seq_len(p), function(j)
      sort(unique(unlist(lapply(object@M[, j], function(d) d@x)))))
  }
  m <- WD.vec.mean(object, w)
  if (p == 1) {
    return(simplifyD(m, if (is.list(support)) support[[1]] else support))
  }
  simplifyD(m, support)
}
