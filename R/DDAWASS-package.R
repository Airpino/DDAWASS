#' DDAWASS: Discrete Distribution Data Analysis with Wasserstein distance
#'
#' Classes and methods for data whose cells are discrete distributions with
#' finite support (frequency distributions of ratings, scores, counts, ...),
#' analysed in the L2 Wasserstein (quantile) metric space. The package is the
#' discrete counterpart of \pkg{HistDAWass}: it provides the
#' \code{\linkS4class{distributionD}} and \code{\linkS4class{MatD}} classes,
#' exact and fast (C++) Wasserstein distances, barycenters, variance,
#' covariance and correlation, k-means and hierarchical clustering,
#' quantile-based PCA and MFA, two-component regression and Cronbach's alpha.
#'
#' Main functions:
#' \itemize{
#'   \item building data: \code{\link{distributionD}}, \code{\link{MatD}},
#'     \code{\link{data2MatD}}, \code{\link{long2MatD}},
#'     \code{\link{MatH2MatD}}, \code{\link{tibble2MatD}};
#'   \item statistics: \code{\link{WD.vec.mean}}, \code{\link{WD.vec.mean.simplified}},
#'     \code{\link{simplifyD}}, \code{\link{WD.var.covar}},
#'     \code{\link{WD.correlation}}, \code{\link{get.MatD.stats}};
#'   \item distances and clustering: \code{\link{WassSqDistD}},
#'     \code{\link{WD_MAT_DIST}}, \code{\link{WD_hclust}},
#'     \code{\link{WD_kmeans}};
#'   \item factorial methods and models: \code{\link{WD.1d.PCA}},
#'     \code{\link{WD.MultiplePCA}}, \code{\link{WD.regression.two.components}},
#'     \code{\link{WD.cronbach}}.
#' }
#' Graphical tools (Green Eyes plots, distributional heatmaps) are in the
#' companion package \pkg{DistrViz}.
#' @keywords internal
"_PACKAGE"

#' Simulated hotel satisfaction survey
#'
#' A \code{\linkS4class{MatD}} with 30 fictitious hotels (rows) and the
#' distributions of the 0-10 scores given by their guests to seven items and
#' to the overall satisfaction (columns). Data are simulated (see
#' \code{data-raw/Satisfaction.R} in the source package);
#' \code{Satisfaction_raw} contains the simulated individual answers.
#' @format \code{Satisfaction}: a \code{MatD} with 30 rows and 8 columns
#'   (\code{Staff}, \code{Cleanliness}, \code{Comfort}, \code{Location},
#'   \code{Price}, \code{Food}, \code{WiFi}, \code{Overall}).
#'   \code{Satisfaction_raw}: a \code{data.frame} with one row per guest.
#' @examples
#' Satisfaction
#' plot(Satisfaction[1:3, 1:3])
#' @name Satisfaction
#' @aliases Satisfaction_raw
"Satisfaction"
