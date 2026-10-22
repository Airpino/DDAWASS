## Building MatD objects from data and converting to/from other formats ---

#' Build a MatD from raw (micro) data
#'
#' Each row of \code{data} is an observation (e.g. a respondent) belonging to
#' a group defined by one or more categorical columns (\code{unit}); for each
#' group and each variable in \code{vars} the frequency distribution of the
#' observed values is computed.
#' @param data a \code{data.frame}.
#' @param unit name(s) or index(es) of the grouping column(s). With several
#'   columns, groups are their combinations (labels joined by \code{sep});
#'   empty combinations are dropped.
#' @param vars names or indices of the discrete numeric variables.
#' @param support (optional) a list (one element per variable, or a single
#'   vector for all) with the admissible values (e.g. \code{1:7} for a
#'   7-point Likert scale): values with zero frequency are then kept in the
#'   support, so that all the distributions of a variable share the same
#'   categories.
#' @param min.n minimum number of observations for a group to be kept
#'   (groups with fewer rows are dropped with a message).
#' @param sep separator used to build the labels of combined groups.
#' @return a \code{MatD} with one row per group. The number of observations
#'   of each distribution is stored in its \code{n} slot (see
#'   \code{get.MatD.stats(x, "n")}).
#' @examples
#' set.seed(2)
#' raw <- data.frame(shop = rep(c("A", "B", "C"), each = 50),
#'                   gender = sample(c("F", "M"), 150, replace = TRUE),
#'                   q1 = sample(1:5, 150, replace = TRUE),
#'                   q2 = sample(0:10, 150, replace = TRUE))
#' data2MatD(raw, unit = "shop", vars = c("q1", "q2"))
#' data2MatD(raw, unit = c("shop", "gender"), vars = c("q1", "q2"),
#'           support = list(1:5, 0:10), min.n = 20)
#' @export
data2MatD <- function(data, unit, vars, support = NULL, min.n = 1, sep = "_") {
  data <- as.data.frame(data)
  if (is.numeric(unit)) unit <- colnames(data)[unit]
  if (!all(unit %in% colnames(data))) stop("unit column(s) not found")
  g <- if (length(unit) == 1) factor(data[[unit]]) else
    interaction(data[unit], sep = sep, drop = TRUE, lex.order = TRUE)
  g <- droplevels(g)
  vn <- if (is.numeric(vars)) colnames(data)[vars] else vars
  if (!all(vn %in% colnames(data))) stop("some vars not found")
  if (min.n > 1) {
    sz <- table(g)
    small <- names(sz)[sz < min.n]
    if (length(small) > 0) {
      message(length(small), " group(s) with fewer than ", min.n,
              " observations dropped: ", paste(small, collapse = ", "))
      keep <- !(g %in% small)
      data <- data[keep, , drop = FALSE]
      g <- droplevels(g[keep])
    }
  }
  lev <- levels(g)
  out <- vector("list", length(lev) * length(vn))
  cnt <- 0
  for (j in seq_along(vn)) {
    sup <- if (is.list(support)) support[[j]] else support
    for (u in lev) {
      cnt <- cnt + 1
      v <- as.numeric(data[[vn[j]]][g == u])
      v <- v[!is.na(v)]
      out[[cnt]] <- if (length(v) == 0) new("distributionD") else if (is.null(sup)) {
        distributionD(v)
      } else {
        if (any(!v %in% sup)) warning("values of ", vn[j], " outside the support were ignored")
        distributionD(sup, counts = vapply(sup, function(s) sum(v == s), 0), keep.zeros = TRUE)
      }
    }
  }
  MatD(out, nrows = length(lev), ncols = length(vn), rownames = lev, varnames = vn)
}

#' Build a MatD from a long frequency table
#'
#' @param data a \code{data.frame} in long format, with one row per
#'   (unit, variable, value) and a column of frequencies.
#' @param unit,variable,value,freq names of the columns with the row
#'   labels, the variable names, the support values and the (absolute or
#'   relative) frequencies.
#' @return a \code{MatD}.
#' @examples
#' tab <- expand.grid(unit = c("u1", "u2"), variable = c("A", "B"), value = 1:3)
#' tab$freq <- c(5, 1, 2, 2, 3, 3, 2, 2, 2, 6, 6, 6)
#' long2MatD(tab)
#' @export
long2MatD <- function(data, unit = "unit", variable = "variable", value = "value",
                      freq = "freq") {
  data <- as.data.frame(data)
  u <- factor(data[[unit]], levels = unique(data[[unit]]))
  v <- factor(data[[variable]], levels = unique(data[[variable]]))
  out <- vector("list", nlevels(u) * nlevels(v))
  cnt <- 0
  for (j in levels(v)) for (i in levels(u)) {
    cnt <- cnt + 1
    sel <- u == i & v == j
    out[[cnt]] <- if (!any(sel) || sum(data[[freq]][sel], na.rm = TRUE) <= 0) new("distributionD") else
      distributionD(data[[value]][sel], counts = data[[freq]][sel])
  }
  MatD(out, nrows = nlevels(u), ncols = nlevels(v), rownames = levels(u), varnames = levels(v))
}

#' Convert a MatD into a long data.frame
#' @param object a \code{MatD}.
#' @return a \code{data.frame} with columns \code{row}, \code{variable}
#'   (factors, in the order of the \code{MatD}), \code{x}, \code{p} and
#'   \code{cdf}.
#' @examples
#' head(MatD2long(Satisfaction))
#' @export
MatD2long <- function(object) {
  rn <- rownames(object@M); cn <- colnames(object@M)
  L <- vector("list", length(object@M))
  cnt <- 0
  for (j in seq_along(cn)) for (i in seq_along(rn)) {
    cnt <- cnt + 1
    d <- object@M[i, j][[1]]
    if (length(d@x) == 0) next
    L[[cnt]] <- data.frame(row = rn[i], variable = cn[j], x = d@x, p = d@p, cdf = d@cdf)
  }
  df <- do.call(rbind, L)
  df$row <- factor(df$row, levels = rn)
  df$variable <- factor(df$variable, levels = cn)
  df
}

#' Convert legacy nested tibbles into a MatD (and back)
#'
#' The research code that preceded this package stored discrete
#' distributional data as a tibble whose cells are \code{data.frame}s with
#' the support in the first column and a \code{freq} column (optionally
#' \code{n} and \code{cdf}). \code{tibble2MatD} converts such objects into a
#' \code{MatD}; \code{MatD2tibble} does the reverse.
#' @param tib a tibble/data.frame of list-columns.
#' @param labels (optional) name or index of a character column holding the
#'   row labels; it is excluded from the variables.
#' @return a \code{MatD} (\code{tibble2MatD}) or a tibble (\code{MatD2tibble}).
#' @examples
#' tb <- MatD2tibble(Satisfaction[1:3, 1:2])
#' tb
#' tibble2MatD(tb, labels = "ID_name")
#' @export
tibble2MatD <- function(tib, labels = NULL) {
  tib <- as.data.frame(tib, stringsAsFactors = FALSE)
  rn <- NULL
  if (!is.null(labels)) {
    rn <- as.character(tib[[labels]])
    tib <- tib[, setdiff(colnames(tib), if (is.numeric(labels)) colnames(tib)[labels] else labels),
               drop = FALSE]
  }
  vars <- colnames(tib)[vapply(tib, is.list, logical(1))]
  n <- nrow(tib)
  out <- vector("list", n * length(vars))
  cnt <- 0
  for (j in vars) for (i in seq_len(n)) {
    cnt <- cnt + 1
    cell <- as.data.frame(tib[[j]][[i]])
    w <- if ("freq" %in% names(cell)) cell$freq else if ("n" %in% names(cell)) cell$n else
      diff(c(0, cell$cdf))
    out[[cnt]] <- distributionD(cell[[1]], p = w)
    if ("n" %in% names(cell)) out[[cnt]]@n <- sum(cell$n)
  }
  MatD(out, nrows = n, ncols = length(vars), rownames = rn, varnames = vars)
}

#' @rdname tibble2MatD
#' @param object a \code{MatD}.
#' @export
MatD2tibble <- function(object) {
  cols <- lapply(seq_len(ncol(object@M)), function(j) {
    lapply(object@M[, j], function(d) {
      data.frame(x = d@x, n = if (is.na(d@n)) NA_real_ else d@p * d@n, freq = d@p, cdf = d@cdf)
    })
  })
  names(cols) <- colnames(object@M)
  tib <- tibble::as_tibble(cols)
  tib$ID_name <- rownames(object@M)
  tib
}

#' Discretize a matrix of histograms (MatH) into a MatD
#'
#' Each histogram of a \code{HistDAWass::MatH} is transformed into a discrete
#' distribution on \code{n_levels} classes: the masses are the probabilities
#' of the classes defined by the breaks. This is the transformation used to
#' compare continuous and discrete distributional data and to draw Green Eyes
#' plots and distributional heatmaps of histogram data.
#' @param object a \code{MatH}.
#' @param n_levels number of classes.
#' @param scheme how the breaks are defined: \code{"variable"} (default,
#'   equal-width classes between the minimum and the maximum of each
#'   variable), \code{"global"} (the same classes for all variables) or
#'   \code{"standardized"} (each variable is standardised with its
#'   Wasserstein mean and standard deviation, and cut into equal-width
#'   classes between -3 and 3, with two open extreme classes).
#' @param breaks (optional) a list with a vector of breaks for each variable
#'   (overrides \code{scheme}).
#' @param values support values of the discrete distributions: \code{"index"}
#'   (default, the class number \code{1, ..., n_levels}) or \code{"midpoint"}
#'   (the class centre in the original scale).
#' @param return.breaks logical; if \code{TRUE} returns a list with the
#'   \code{MatD} and the breaks used.
#' @return a \code{MatD} (or a list).
#' @examples
#' library(HistDAWass)
#' B <- MatH2MatD(BLOOD, n_levels = 10)
#' B
#' @export
MatH2MatD <- function(object, n_levels = 10, scheme = c("variable", "global", "standardized"),
                      breaks = NULL, values = c("index", "midpoint"), return.breaks = FALSE) {
  if (!is(object, "MatH")) stop("object must be a MatH")
  scheme <- match.arg(scheme)
  values <- match.arg(values)
  M <- object@M
  nr <- nrow(M); nc <- ncol(M)
  mins <- vapply(seq_len(nc), function(j) min(vapply(M[, j], function(d) d@x[1], 0)), 0)
  maxs <- vapply(seq_len(nc), function(j) max(vapply(M[, j], function(d) d@x[length(d@x)], 0)), 0)
  ctr <- rep(0, nc); scl <- rep(1, nc)
  if (is.null(breaks)) {
    breaks <- switch(scheme,
      variable = lapply(seq_len(nc), function(j) seq(mins[j], maxs[j], length.out = n_levels + 1)),
      global = rep(list(seq(min(mins), max(maxs), length.out = n_levels + 1)), nc),
      standardized = {
        for (j in seq_len(nc)) {
          ctr[j] <- HistDAWass::WH.vec.mean(object[, j])@m
          scl[j] <- sqrt(as.numeric(HistDAWass::WH.var.covar(object[, j])))
        }
        rep(list(c(-Inf, seq(-3, 3, length.out = n_levels - 1), Inf)), nc)
      })
  } else {
    if (!is.list(breaks)) breaks <- rep(list(breaks), nc)
  }
  out <- vector("list", nr * nc)
  cnt <- 0
  for (j in seq_len(nc)) {
    br <- breaks[[j]]
    L <- length(br) - 1
    brs <- br * scl[j] + ctr[j]   # breaks in the original scale
    sup <- if (values == "index") seq_len(L) else {
      fin <- brs
      fin[!is.finite(fin)] <- c(mins[j], maxs[j])[(fin[!is.finite(fin)] > 0) + 1]
      (fin[-1] + fin[-length(fin)]) / 2
    }
    for (i in seq_len(nr)) {
      cnt <- cnt + 1
      d <- M[i, j][[1]]
      cdf <- vapply(brs, function(q) {
        if (q == -Inf) 0 else if (q == Inf) 1 else HistDAWass::compP(d, q)
      }, 0)
      cdf[1] <- 0; cdf[length(cdf)] <- 1
      out[[cnt]] <- distributionD(sup, p = pmax(diff(cdf), 0))
    }
  }
  res <- MatD(out, nrows = nr, ncols = nc, rownames = rownames(M), varnames = colnames(M))
  if (return.breaks) {
    names(breaks) <- colnames(M)
    return(list(MatD = res, breaks = lapply(seq_len(nc), function(j) breaks[[j]] * scl[j] + ctr[j]),
                scheme = scheme))
  }
  res
}
