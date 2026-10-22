## Basic Wasserstein statistics of discrete distributional variables ----

#' Barycenter (Frechet mean) of discrete distributional variables
#'
#' Computes, for each column of a \code{\linkS4class{MatD}}, the (weighted)
#' L2 Wasserstein barycenter of its distributions, i.e. the distribution
#' whose quantile function is the weighted average of the quantile functions.
#' The barycenter of discrete distributions is itself discrete; its support
#' is in general not a subset of the original supports.
#' @param object a \code{MatD}.
#' @param w (optional) weights of the rows (default: equal weights).
#' @return a \code{distributionD} if \code{object} has one column, otherwise
#'   a one-row \code{MatD}.
#' @examples
#' WD.vec.mean(Satisfaction[, "Staff"])
#' WD.vec.mean(Satisfaction)
#' @export
WD.vec.mean <- function(object, w = 1) {
  .check_nonempty(object)
  n <- nrow(object@M); p <- ncol(object@M)
  w <- .norm_w(w, n)
  out <- lapply(seq_len(p), function(j) {
    xc <- .col_xc(object, j)
    b <- .wd_barycenter_c(xc$x, xc$c, w)
    .compactD(b$x, b$cdf)
  })
  if (p == 1) return(out[[1]])
  MatD(out, nrows = 1, ncols = p, rownames = "Mean", varnames = colnames(object@M))
}

#' Wasserstein variance-covariance and correlation matrices
#'
#' \code{WD.var.covar} computes the covariance matrix of the discrete
#' distributional variables of a \code{MatD}, following the L2 Wasserstein
#' (quantile) approach of Irpino and Verde (2015):
#' \deqn{C_{jk} = \sum_i w_i \int_0^1 (Q_{ij}(t) - \bar{Q}_j(t))(Q_{ik}(t) - \bar{Q}_k(t)) dt}
#' where \eqn{\bar{Q}_j} is the quantile function of the barycenter of
#' variable \eqn{j}. The diagonal contains the Wasserstein (Frechet)
#' variances. \code{WD.correlation} returns the corresponding correlation
#' matrix.
#' @param object a \code{MatD}.
#' @param w (optional) weights of the rows.
#' @return a \code{p x p} matrix.
#' @references Irpino A., Verde R. (2015). Basic statistics for distributional
#'   symbolic variables: a new metric-based approach. \emph{Advances in Data
#'   Analysis and Classification}, 9(2), 143--175.
#'   \doi{10.1007/s11634-014-0176-4}
#' @examples
#' WD.var.covar(Satisfaction)
#' round(WD.correlation(Satisfaction), 3)
#' @export
WD.var.covar <- function(object, w = 1) {
  .check_nonempty(object)
  n <- nrow(object@M); p <- ncol(object@M)
  w <- .norm_w(w, n)
  bary <- lapply(seq_len(p), function(j) {
    xc <- .col_xc(object, j)
    .wd_barycenter_c(xc$x, xc$c, w)
  })
  C <- matrix(0, p, p, dimnames = list(colnames(object@M), colnames(object@M)))
  for (j in seq_len(p)) {
    xj <- .col_xc(object, j)
    for (k in j:p) {
      xk <- if (k == j) xj else .col_xc(object, k)
      s <- 0
      for (i in seq_len(n)) {
        s <- s + w[i] * .wd_dot_c(xj$x[[i]], xj$c[[i]], xk$x[[i]], xk$c[[i]])
      }
      s <- s - .wd_dot_c(bary[[j]]$x, bary[[j]]$cdf, bary[[k]]$x, bary[[k]]$cdf)
      C[j, k] <- C[k, j] <- s
    }
  }
  C
}

#' @rdname WD.var.covar
#' @export
WD.correlation <- function(object, w = 1) {
  C <- WD.var.covar(object, w)
  s <- sqrt(diag(C))
  R <- C / outer(s, s)
  diag(R) <- 1
  R
}

#' Total sum of squares of a discrete distributional variable
#'
#' \code{WD.SSQ} returns the (weighted) sum of squared Wasserstein distances
#' of the distributions of each column from their barycenter (the deviance).
#' @param object a \code{MatD}.
#' @param w (optional) weights of the rows (not normalised).
#' @return a named numeric vector (one value per column).
#' @examples
#' WD.SSQ(Satisfaction[, 1:3])
#' @export
WD.SSQ <- function(object, w = 1) {
  n <- nrow(object@M)
  ww <- if (length(w) == 1) rep(1, n) else w
  v <- diag(WD.var.covar(object, ww))
  v * sum(ww)
}

## Distances ---------------------------------------------------------------

#' Matrix of L2 Wasserstein distances between the rows of a MatD
#'
#' The distance between two rows (units) is
#' \eqn{d(i, i') = \sqrt{\sum_j \lambda_j d_W^2(x_{ij}, x_{i'j})}}.
#' The computation is done in C++ with an exact linear-time algorithm.
#' @param object a \code{MatD}.
#' @param simplify ignored (kept for compatibility with \code{WH_MAT_DIST}).
#' @param weights (optional) weights \eqn{\lambda_j} of the variables.
#' @param squared logical: return squared distances? Default \code{FALSE}.
#' @return a symmetric \code{n x n} matrix.
#' @examples
#' D <- WD_MAT_DIST(Satisfaction)
#' round(D[1:5, 1:5], 3)
#' @export
WD_MAT_DIST <- function(object, simplify = FALSE, weights = 1, squared = FALSE) {
  .check_nonempty(object)
  p <- ncol(object@M)
  w <- if (length(weights) == 1) rep(weights, p) else weights
  if (length(w) != p) stop("weights must have one value per variable")
  xc <- .all_xc(object)
  D <- .wd_distmat_c(xc$x, xc$c, as.numeric(w))
  D[D < 0] <- 0
  if (!squared) D <- sqrt(D)
  dimnames(D) <- list(rownames(object@M), rownames(object@M))
  D
}

## Clustering ---------------------------------------------------------------

#' Hierarchical clustering of discrete distributional data
#'
#' Agglomerative clustering (\code{\link[stats]{hclust}}) based on the L2
#' Wasserstein distance between rows.
#' @param x a \code{MatD}.
#' @param simplify ignored (compatibility with \code{WH_hclust}).
#' @param method agglomeration method passed to \code{hclust} (default
#'   \code{"complete"}; \code{"ward.D2"} is the natural choice for the
#'   Wasserstein metric).
#' @param weights (optional) weights of the variables.
#' @return an object of class \code{hclust}.
#' @examples
#' hc <- WD_hclust(Satisfaction, method = "ward.D2")
#' plot(hc)
#' @export
WD_hclust <- function(x, simplify = FALSE, method = "complete", weights = 1) {
  D <- WD_MAT_DIST(x, weights = weights)
  stats::hclust(stats::as.dist(D), method = method)
}

#' K-means clustering of discrete distributional data
#'
#' Dynamic clustering (k-means) of the rows of a \code{\linkS4class{MatD}}
#' using the L2 Wasserstein distance; prototypes are the Wasserstein
#' barycenters of the clusters. Initial prototypes are chosen with the
#' k-means++ seeding, the best of \code{rep} runs is returned.
#' @param x a \code{MatD}.
#' @param k number of clusters.
#' @param rep number of repetitions (random starts).
#' @param simplify ignored (compatibility with \code{WH_kmeans}).
#' @param qua ignored (compatibility with \code{WH_kmeans}).
#' @param standardize logical; if \code{TRUE} each variable is weighted by
#'   the inverse of its Wasserstein variance.
#' @param weights (optional) weights of the variables (overridden by
#'   \code{standardize}).
#' @param iter.max maximum number of iterations.
#' @param verbose print progress?
#' @return a list with \code{solution} (a list with \code{IDX}, the cluster
#'   labels, \code{cardinality}, \code{centers}, a \code{MatD} of prototypes,
#'   and \code{Crit}, the within sum of squares), \code{quality} (the ratio
#'   between the between and the total sum of squares), \code{TSS} and
#'   \code{all.crit}, the criteria of all the repetitions.
#' @examples
#' set.seed(1)
#' km <- WD_kmeans(Satisfaction, k = 3, rep = 5)
#' km$solution$IDX
#' km$quality
#' @export
WD_kmeans <- function(x, k, rep = 5, simplify = FALSE, qua = 10,
                      standardize = FALSE, weights = 1, iter.max = 100,
                      verbose = FALSE) {
  .check_nonempty(x)
  n <- nrow(x@M); p <- ncol(x@M)
  k <- as.integer(k)
  if (k < 2 || k >= n) stop("k must be between 2 and n - 1")
  w <- if (length(weights) == 1) rep(weights, p) else weights
  if (standardize) w <- 1 / diag(WD.var.covar(x))
  xc <- .all_xc(x)

  # total sum of squares
  bary <- lapply(seq_len(p), function(j) .wd_barycenter_c(xc$x[[j]], xc$c[[j]], rep(1 / n, n)))
  TSS <- sum(.wd_dist_to_protos_c(xc$x, xc$c,
                                  lapply(bary, function(b) list(b$x)),
                                  lapply(bary, function(b) list(b$cdf)), w))

  protos_of <- function(ids) {
    px <- pc <- vector("list", p)
    for (j in seq_len(p)) {
      px[[j]] <- pc[[j]] <- vector("list", k)
      for (h in seq_len(k)) {
        sel <- which(ids == h)
        if (length(sel) == 0) sel <- sample.int(n, 1)
        b <- .wd_barycenter_c(xc$x[[j]][sel], xc$c[[j]][sel], rep(1 / length(sel), length(sel)))
        px[[j]][[h]] <- b$x
        pc[[j]][[h]] <- b$cdf
      }
    }
    list(x = px, c = pc)
  }

  best <- NULL
  all.crit <- numeric(rep)
  for (r in seq_len(rep)) {
    # k-means++ seeding
    seeds <- sample.int(n, 1)
    dmin <- .wd_dist_to_protos_c(xc$x, xc$c,
                                 lapply(xc$x, function(col) col[seeds]),
                                 lapply(xc$c, function(col) col[seeds]), w)[, 1]
    while (length(seeds) < k) {
      prob <- dmin / sum(dmin)
      if (!is.finite(sum(prob)) || sum(prob) == 0) prob <- rep(1 / n, n)
      s <- sample.int(n, 1, prob = prob)
      seeds <- c(seeds, s)
      dnew <- .wd_dist_to_protos_c(xc$x, xc$c,
                                   lapply(xc$x, function(col) col[s]),
                                   lapply(xc$c, function(col) col[s]), w)[, 1]
      dmin <- pmin(dmin, dnew)
    }
    prot <- list(x = lapply(xc$x, function(col) col[seeds]),
                 c = lapply(xc$c, function(col) col[seeds]))
    ids <- rep(0L, n)
    crit_old <- Inf
    for (it in seq_len(iter.max)) {
      D <- .wd_dist_to_protos_c(xc$x, xc$c, prot$x, prot$c, w)
      ids <- max.col(-D, ties.method = "first")
      # empty clusters: move the farthest point
      for (h in seq_len(k)) {
        if (!any(ids == h)) {
          far <- which.max(D[cbind(seq_len(n), ids)])
          ids[far] <- h
        }
      }
      prot <- protos_of(ids)
      D <- .wd_dist_to_protos_c(xc$x, xc$c, prot$x, prot$c, w)
      crit <- sum(D[cbind(seq_len(n), ids)])
      if (verbose) message(sprintf("rep %d it %d crit %.6g", r, it, crit))
      if (crit_old - crit <= 1e-10 * max(1, crit)) break
      crit_old <- crit
    }
    all.crit[r] <- crit
    if (is.null(best) || crit < best$crit) best <- list(ids = ids, prot = prot, crit = crit)
  }

  cent <- vector("list", k * p)
  cnt <- 0
  for (j in seq_len(p)) for (h in seq_len(k)) {
    cnt <- cnt + 1
    cent[[cnt]] <- .compactD(best$prot$x[[j]][[h]], best$prot$c[[j]][[h]])
  }
  centers <- MatD(cent, nrows = k, ncols = p, rownames = paste0("Cl.", seq_len(k)),
                  varnames = colnames(x@M))
  IDX <- stats::setNames(best$ids, rownames(x@M))
  list(solution = list(IDX = IDX, cardinality = table(IDX), centers = centers,
                       Crit = best$crit),
       quality = 1 - best$crit / TSS, TSS = TSS, all.crit = all.crit)
}

## Reliability ---------------------------------------------------------------

#' Cronbach's alpha for discrete distributional variables
#'
#' Computes Cronbach's alpha using the Wasserstein covariance matrix of a set
#' of distributional items (\code{GLOB}), and its decomposition into the part
#' due to the covariance of the means (\code{AVE}) and the part due to the
#' covariance of the centred distributions (\code{VAR}). The
#' \code{alpha.if.deleted} matrix reports the values obtained removing each
#' item in turn.
#' @param object a \code{MatD} with at least three columns.
#' @return a list with \code{alpha} and \code{alpha.if.deleted}.
#' @examples
#' WD.cronbach(Satisfaction)
#' @export
WD.cronbach <- function(object) {
  p <- ncol(object@M)
  if (p < 3) stop("at least three items are needed")
  n <- nrow(object@M)
  C <- WD.var.covar(object)
  M <- get.MatD.stats(object, "mean")$mat
  Cm <- stats::cov(M) * (n - 1) / n
  Cd <- C - Cm
  alpha <- function(S) {
    q <- ncol(S)
    q / (q - 1) * (1 - sum(diag(S)) / sum(S))
  }
  a <- c(GLOB = alpha(C), AVE = alpha(Cm), VAR = alpha(Cd))
  del <- t(vapply(seq_len(p), function(j) {
    c(GLOB = alpha(C[-j, -j]), AVE = alpha(Cm[-j, -j]), VAR = alpha(Cd[-j, -j]))
  }, numeric(3)))
  rownames(del) <- colnames(object@M)
  list(alpha = a, alpha.if.deleted = del)
}
