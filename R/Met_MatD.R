## Methods for MatD ------------------------------------------------------

#' Accessors of a MatD object
#'
#' @param object a \code{\linkS4class{MatD}} object.
#' @return \code{get.MatD.nrows} and \code{get.MatD.ncols} return integers,
#'   \code{get.MatD.rownames} and \code{get.MatD.varnames} character
#'   vectors, \code{get.cell.MatD} a \code{distributionD}.
#' @examples
#' get.MatD.nrows(Satisfaction)
#' get.MatD.varnames(Satisfaction)
#' get.cell.MatD(Satisfaction, 1, 1)
#' @name MatD-accessors
#' @export
get.MatD.nrows <- function(object) nrow(object@M)
#' @rdname MatD-accessors
#' @export
get.MatD.ncols <- function(object) ncol(object@M)
#' @rdname MatD-accessors
#' @export
get.MatD.rownames <- function(object) rownames(object@M)
#' @rdname MatD-accessors
#' @export
get.MatD.varnames <- function(object) colnames(object@M)
#' @rdname MatD-accessors
#' @param r,c row and column index.
#' @export
get.cell.MatD <- function(object, r, c) object@M[r, c][[1]]

#' @rdname MatD-accessors
#' @param value a \code{distributionD}.
#' @export
set.cell.MatD <- function(object, r, c, value) {
  if (!is(value, "distributionD")) stop("value must be a distributionD")
  object@M[r, c][[1]] <- value
  object
}

setMethod("show", "MatD", function(object) {
  nr <- nrow(object@M); nc <- ncol(object@M)
  cat("a matrix of discrete distributions\n", nr, " rows x ", nc, " columns\n", sep = "")
  showr <- min(nr, 10); showc <- min(nc, 6)
  out <- matrix("", showr, showc,
                dimnames = list(rownames(object@M)[seq_len(showr)],
                                colnames(object@M)[seq_len(showc)]))
  for (i in seq_len(showr)) for (j in seq_len(showc)) {
    d <- object@M[i, j][[1]]
    out[i, j] <- if (length(d@x) == 0) "(empty)" else
      sprintf("[m=%.3g s=%.3g k=%d]", d@m, d@s, length(d@x))
  }
  print(noquote(out))
  if (nr > showr || nc > showc) cat("... (only the first rows/columns shown)\n")
  cat("k = number of support values\n")
  invisible(NULL)
})

#' Extract a sub-matrix of a MatD
#' @param x a \code{MatD}.
#' @param i,j row and column indices (numbers, names or logicals).
#' @param ... unused.
#' @param drop ignored: the result is always a \code{MatD}.
#' @return a \code{MatD}.
#' @examples
#' Satisfaction[1:3, c("Staff", "Price")]
#' @aliases [,MatD,ANY,ANY,ANY-method
#' @exportMethod [
setMethod("[", signature(x = "MatD", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, ..., drop = FALSE) {
    if (missing(i)) i <- seq_len(nrow(x@M))
    if (missing(j)) j <- seq_len(ncol(x@M))
    new("MatD", M = x@M[i, j, drop = FALSE])
  })

#' Bind MatD objects by rows or by columns
#' @param a,b two \code{MatD} objects.
#' @return a \code{MatD}.
#' @examples
#' WD.bind.col(Satisfaction[, 1:2], Satisfaction[, 3])
#' @name WD.bind
#' @export
WD.bind.col <- function(a, b) {
  if (nrow(a@M) != nrow(b@M)) stop("different number of rows")
  new("MatD", M = cbind(a@M, b@M))
}
#' @rdname WD.bind
#' @export
WD.bind.row <- function(a, b) {
  if (ncol(a@M) != ncol(b@M)) stop("different number of columns")
  new("MatD", M = rbind(a@M, b@M))
}

#' Statistics of the distributions in a MatD
#'
#' Returns a matrix with a summary statistic computed for each cell.
#' @param object a \code{MatD}.
#' @param stat one of \code{"mean"}, \code{"median"}, \code{"min"},
#'   \code{"max"}, \code{"std"}, \code{"var"}, \code{"skewness"},
#'   \code{"kurtosis"}, \code{"quantile"}, \code{"n"} (number of
#'   observations) or \code{"k"} (number of support values).
#' @param prob probability level used when \code{stat = "quantile"}.
#' @return a list with elements \code{stat} and \code{mat}.
#' @examples
#' get.MatD.stats(Satisfaction)$mat[1:5, ]
#' get.MatD.stats(Satisfaction, stat = "quantile", prob = 0.1)$mat[1:5, ]
#' @export
get.MatD.stats <- function(object, stat = "mean", prob = 0.5) {
  f <- switch(stat,
    mean = function(d) d@m,
    std = function(d) d@s,
    var = function(d) d@s^2,
    skewness = skewD,
    kurtosis = kurtD,
    median = function(d) compQ(d, 0.5),
    min = function(d) d@x[1],
    max = function(d) d@x[length(d@x)],
    quantile = function(d) compQ(d, prob),
    n = function(d) d@n,
    k = function(d) length(d@x),
    stop("unknown stat")
  )
  mat <- apply(object@M, c(1, 2), function(cell) {
    d <- cell[[1]]
    if (length(d@x) == 0) NA_real_ else f(d)
  })
  mat <- matrix(mat, nrow(object@M), ncol(object@M),
                dimnames = dimnames(object@M))
  list(stat = stat, mat = mat)
}

#' Plot a MatD
#'
#' Draws a grid of bar charts (rows = units, columns = variables).
#' @param x a \code{MatD}.
#' @param y unused.
#' @param type \code{"BAR"} (default) or \code{"CDF"}.
#' @param col fill colour.
#' @param ... unused.
#' @return a \code{ggplot} object.
#' @examples
#' plot(Satisfaction[1:4, 1:4])
#' @aliases plot,MatD-method
setMethod("plot", signature(x = "MatD"), function(x, y, type = "BAR", col = "#2E7D32", ...) {
  df <- MatD2long(x)
  g <- ggplot2::ggplot(df)
  if (toupper(type) == "CDF") {
    g <- g + ggplot2::geom_step(ggplot2::aes(x = .data$x, y = .data$cdf),
                                colour = col, direction = "hv")
  } else {
    g <- g + ggplot2::geom_col(ggplot2::aes(x = .data$x, y = .data$p), fill = col)
  }
  g + ggplot2::facet_grid(rows = ggplot2::vars(.data$row),
                          cols = ggplot2::vars(.data$variable), scales = "free_x") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
})

## internal helpers ------------------------------------------------------

# list of supports and cdfs for column j
.col_xc <- function(object, j) {
  cells <- object@M[, j]
  list(x = lapply(cells, function(d) d@x), c = lapply(cells, function(d) d@cdf))
}

.all_xc <- function(object) {
  p <- ncol(object@M)
  L <- lapply(seq_len(p), function(j) .col_xc(object, j))
  list(x = lapply(L, `[[`, "x"), c = lapply(L, `[[`, "c"))
}

.check_nonempty <- function(object) {
  e <- vapply(object@M, function(d) length(d@x) == 0, logical(1))
  if (any(e)) stop("the MatD contains empty distributions (", sum(e), " cells)")
  invisible(TRUE)
}

.norm_w <- function(w, n) {
  if (length(w) == 1) return(rep(1 / n, n))
  if (length(w) != n) stop("wrong length of the weights")
  if (any(w < 0)) stop("weights must be non-negative")
  w / sum(w)
}
