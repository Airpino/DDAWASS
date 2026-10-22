## Principal components of discrete distributional variables -------------

#' Matrix of quantiles of a MatD
#'
#' For each row and each selected variable, the quantiles at levels
#' \code{seq(pp[1], pp[2], length.out = quantiles + 1)} are computed and
#' arranged in a numeric matrix (one block of \code{quantiles + 1} columns
#' per variable).
#' @param object a \code{MatD}.
#' @param quantiles number of quantile intervals.
#' @param vars columns to use (default all).
#' @param pp range of probability levels (default \code{c(0, 1)}).
#' @return a numeric matrix.
#' @examples
#' WD.quantile.matrix(Satisfaction, quantiles = 4, vars = 1:2)[1:5, ]
#' @export
WD.quantile.matrix <- function(object, quantiles = 10, vars = seq_len(ncol(object@M)),
                               pp = c(0, 1)) {
  vars <- .var_index(object, vars)
  .check_nonempty(object[, vars])
  p <- seq(pp[1], pp[2], length.out = quantiles + 1)
  lab <- c("Min", sprintf("Q(%.2f)", p[-c(1, length(p))]), "Max")
  if (pp[1] > 0) lab[1] <- sprintf("Q(%.2f)", pp[1])
  if (pp[2] < 1) lab[length(lab)] <- sprintf("Q(%.2f)", pp[2])
  vn <- colnames(object@M)[vars]
  blocks <- lapply(seq_along(vars), function(v) {
    j <- vars[v]
    Q <- t(vapply(object@M[, j], function(d) compQ(d, p), numeric(length(p))))
    if (length(p) == 1) Q <- t(Q)
    colnames(Q) <- paste(vn[v], lab, sep = ".")
    Q
  })
  Q <- do.call(cbind, blocks)
  rownames(Q) <- rownames(object@M)
  Q
}

#' Principal component analysis of a single discrete distributional variable
#'
#' Performs the quantile-based PCA of one distributional variable (Verde and
#' Irpino, 2010): a classical non-standardised PCA (with
#' \code{\link[FactoMineR]{PCA}}) of the matrix of quantiles, with the means
#' and the standard deviations of the distributions as supplementary
#' variables.
#' @param data a \code{MatD}.
#' @param var the column to analyse.
#' @param quantiles number of quantile intervals (default 10).
#' @param plots logical; draw the \pkg{FactoMineR} plots?
#' @param ncp number of components kept.
#' @return the \code{PCA} object of \pkg{FactoMineR} with the additional
#'   element \code{quantile_matrix}.
#' @examples
#' res <- WD.1d.PCA(Satisfaction, var = "Staff", plots = FALSE)
#' res$eig[1:3, ]
#' WD.plot.Spanish.fan(res)
#' @export
WD.1d.PCA <- function(data, var = 1, quantiles = 10, plots = TRUE, ncp = 5) {
  if (length(var) > 1) stop("var must be a single column")
  var <- .var_index(data, var)
  n <- nrow(data@M)
  if (n < 2) stop("at least two rows are needed")
  Q <- WD.quantile.matrix(data, quantiles, var)
  st <- data@M[, var]
  MF <- cbind(Q / sqrt(quantiles + 1),
              Mean = vapply(st, function(d) d@m, 0),
              Std = vapply(st, function(d) d@s, 0))
  MF <- .jitter_constant(MF)
  res <- FactoMineR::PCA(MF, quanti.sup = c(ncol(MF) - 1, ncol(MF)),
                         scale.unit = FALSE, graph = plots, ncp = ncp)
  res$quantile_matrix <- Q
  res
}

#' Multiple factor analysis of discrete distributional variables
#'
#' Quantile-based multiple factor analysis (MFA) of several distributional
#' variables: each variable is represented by a group of quantile columns
#' and the groups are balanced as in \code{\link[FactoMineR]{MFA}}.
#' @param data a \code{MatD}.
#' @param list.of.vars the columns to analyse (at least two).
#' @param quantiles number of quantile intervals (default 10).
#' @param ncp number of components kept.
#' @param pp range of probability levels.
#' @param plots logical; draw the \pkg{FactoMineR} plots?
#' @references Crocetta C., Irpino A., Antonucci L., Marin C. (2025). A
#'   performance indicator and its decomposition according to the impacts of
#'   different aspects based on distributional data. \emph{Quality & Quantity},
#'   59(3), 2091--2110. \doi{10.1007/s11135-024-01956-0}
#' @return the \code{MFA} object with the additional element
#'   \code{quantile_matrix}.
#' @examples
#' res <- WD.MultiplePCA(Satisfaction, 1:4, quantiles = 5, plots = FALSE)
#' res$eig[1:3, ]
#' @export
WD.MultiplePCA <- function(data, list.of.vars, quantiles = 10, ncp = 5,
                           pp = c(0, 1), plots = TRUE) {
  if (length(list.of.vars) < 2) stop("at least two variables are needed")
  list.of.vars <- .var_index(data, list.of.vars)
  Q <- WD.quantile.matrix(data, quantiles, list.of.vars, pp)
  MF <- .jitter_constant(Q / sqrt(quantiles + 1))
  vn <- colnames(data@M)[list.of.vars]
  res <- FactoMineR::MFA(as.data.frame(MF), group = rep(quantiles + 1, length(list.of.vars)),
                         type = rep("c", length(list.of.vars)), ncp = ncp,
                         name.group = vn, graph = plots)
  res$quantile_matrix <- Q
  res
}

# MFA/PCA fail with constant columns (frequent with discrete data: e.g. all
# minima equal to 0): add a negligible noise
.jitter_constant <- function(M) {
  for (cc in seq_len(ncol(M))) {
    if (stats::sd(M[, cc]) < 1e-10) M[, cc] <- M[, cc] + stats::rnorm(nrow(M), 0, 1e-8)
  }
  M
}

#' Spanish-fan plot of a quantile PCA
#'
#' Draws the plot of the quantile variables of a \code{\link{WD.1d.PCA}}
#' as consecutive coloured triangles (the "Spanish fan"): colours go from
#' red (extreme quantiles) to yellow (median).
#' @param res the output of \code{WD.1d.PCA}.
#' @param axes the two axes to plot.
#' @return a \code{ggplot} object.
#' @examples
#' res <- WD.1d.PCA(Satisfaction, var = "Price", plots = FALSE)
#' WD.plot.Spanish.fan(res, axes = c(1, 2))
#' @export
WD.plot.Spanish.fan <- function(res, axes = c(1, 2)) {
  CV <- res$var$coord[, axes, drop = FALSE]
  nq <- nrow(CV)
  centr <- (nq + 1) / 2
  poly <- do.call(rbind, lapply(2:nq, function(t) {
    g <- 1 - abs(t - 0.5 - centr) / (centr - 1)
    data.frame(id = t, x = c(0, CV[t - 1, 1], CV[t, 1]), y = c(0, CV[t - 1, 2], CV[t, 2]),
               fill = grDevices::rgb(1, max(min(g, 1), 0), 0))
  }))
  labs <- data.frame(x = CV[, 1], y = CV[, 2], lab = rownames(CV))
  labs <- labs[unique(round(seq(1, nq, length.out = min(nq, 5)))), ]
  ev <- res$eig[axes, 2]
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_polygon(data = poly, ggplot2::aes(x = .data$x, y = .data$y,
                                                     group = .data$id, fill = .data$fill),
                          colour = "black", alpha = 0.8) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(data = labs, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                       size = 3, vjust = -0.5) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = sprintf("Axis %d (%.1f%%)", axes[1], ev[1]),
                  y = sprintf("Axis %d (%.1f%%)", axes[2], ev[2]),
                  title = "Spanish-fan plot of the quantiles") +
    ggplot2::theme_minimal()
}

# convert variable names or indices into column indices
.var_index <- function(object, vars) {
  if (is.character(vars)) {
    idx <- match(vars, colnames(object@M))
    if (anyNA(idx)) stop("unknown variable(s): ", paste(vars[is.na(idx)], collapse = ", "))
    return(idx)
  }
  vars
}
