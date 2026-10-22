brute_w2 <- function(a, b, N = 200000) {
  t <- (seq_len(N) - 0.5) / N
  mean((compQ(a, t) - compQ(b, t))^2)
}

test_that("constructor merges, sorts and normalises", {
  d <- distributionD(c(3, 1, 3, 2), p = c(1, 1, 1, 1))
  expect_equal(d@x, c(1, 2, 3))
  expect_equal(d@p, c(0.25, 0.25, 0.5))
  expect_equal(d@cdf[3], 1)
  r <- distributionD(c(1, 1, 2, 5, NA))
  expect_equal(r@n, 4)
  expect_equal(meanD(r), mean(c(1, 1, 2, 5)))
  expect_equal(stdD(r)^2, mean((c(1, 1, 2, 5) - 2.25)^2))
})

test_that("Wasserstein distance is exact and decomposes", {
  set.seed(1)
  a <- distributionD(0:10, counts = rpois(11, 10))
  b <- distributionD(c(1.5, 3, 7), p = c(0.2, 0.5, 0.3))
  expect_equal(WassSqDistD(a, b), brute_w2(a, b), tolerance = 1e-4)
  expect_equal(WassSqDistD(a, a), 0)
  expect_equal(WassSqDistD(a, b), WassSqDistD(b, a))
  det <- WassSqDistD(a, b, details = TRUE)
  expect_equal(unname(det["SQ.W"]), unname(det["POSITION"] + det["SIZE"] + det["SHAPE"]))
  expect_equal(unname(det["rQQ"]), rQQ(a, b), tolerance = 1e-10)
})

test_that("quantiles and cdf are consistent", {
  d <- distributionD(1:4, p = c(0.1, 0.4, 0.3, 0.2))
  expect_equal(compQ(d, c(0, 0.1, 0.10001, 0.5, 0.8, 1)), c(1, 1, 2, 2, 3, 4))
  expect_equal(compP(d, c(0, 1, 2.5, 4, 10)), c(0, 0.1, 0.5, 1, 1))
})

test_that("barycenter variance equals mean squared distance", {
  X <- Satisfaction[, "Staff"]
  m <- WD.vec.mean(X)
  v <- mean(vapply(X@M[, 1], function(d) WassSqDistD(d, m), 0))
  expect_equal(as.numeric(WD.var.covar(X)), v, tolerance = 1e-10)
  expect_equal(meanD(m), mean(get.MatD.stats(X)$mat))
})

test_that("covariance matrix is symmetric PSD and correlation in [-1, 1]", {
  C <- WD.var.covar(Satisfaction)
  expect_true(isSymmetric(C))
  expect_true(min(eigen(C, symmetric = TRUE)$values) > -1e-10)
  R <- WD.correlation(Satisfaction)
  expect_true(all(abs(R) <= 1 + 1e-10))
})

test_that("distance matrix agrees with pairwise distances", {
  X <- Satisfaction[1:5, 1:3]
  D <- WD_MAT_DIST(X, squared = TRUE)
  s <- sum(vapply(1:3, function(j) WassSqDistD(X@M[2, j][[1]], X@M[4, j][[1]]), 0))
  expect_equal(D[2, 4], s)
  expect_equal(diag(D), rep(0, 5), ignore_attr = TRUE)
})

test_that("clustering works", {
  set.seed(3)
  km <- WD_kmeans(Satisfaction, 3, rep = 2)
  expect_equal(length(unique(km$solution$IDX)), 3)
  expect_true(km$quality > 0 && km$quality < 1)
  hc <- WD_hclust(Satisfaction)
  expect_s3_class(hc, "hclust")
})

test_that("regression gives non-negative shape coefficients and a good fit", {
  xv <- c("Staff", "Cleanliness", "Comfort", "Price")
  pars <- WD.regression.two.components(Satisfaction, "Overall", xv)
  expect_true(all(pars[grep("^CEN_", names(pars))] >= 0))
  pred <- WD.regression.two.components.predict(Satisfaction[, xv], pars)
  gof <- WD.regression.GOF(Satisfaction[, "Overall"], pred)
  expect_true(gof$OMEGA > 0.8)
})

test_that("NNLS solver matches a brute-force solution", {
  A <- matrix(c(2, 1, 1, 2), 2)
  b <- c(-1, 3)
  expect_equal(DDAWASS:::.nnls_gram(A, b), c(0, 1.5), tolerance = 1e-10)
})

test_that("converters round-trip", {
  tb <- MatD2tibble(Satisfaction[1:4, 1:3])
  back <- tibble2MatD(tb, labels = "ID_name")
  expect_equal(get.MatD.stats(back)$mat, get.MatD.stats(Satisfaction[1:4, 1:3])$mat)
  M2 <- data2MatD(Satisfaction_raw, "hotel", c("Staff", "Price"))
  expect_equal(get.MatD.stats(M2)$mat, get.MatD.stats(Satisfaction[, c("Staff", "Price")])$mat)
  lg <- MatD2long(Satisfaction[1:2, 1:2])
  names(lg) <- c("unit", "variable", "value", "freq", "cdf")
  expect_equal(get.MatD.stats(long2MatD(lg))$mat, get.MatD.stats(Satisfaction[1:2, 1:2])$mat,
               ignore_attr = TRUE)
})

test_that("MatH discretization preserves the mass", {
  B <- MatH2MatD(HistDAWass::BLOOD, n_levels = 8)
  expect_equal(get.MatD.nrows(B), 14)
  expect_true(all(get.MatD.stats(B, "max")$mat <= 8))
})

test_that("variables can be selected by name in PCA functions", {
  Q1 <- WD.quantile.matrix(Satisfaction, 4, vars = c("Staff", "Price"))
  Q2 <- WD.quantile.matrix(Satisfaction, 4, vars = c(1, 5))
  expect_equal(Q1, Q2)
  expect_false(any(grepl("^NA", colnames(Q1))))
  res <- WD.MultiplePCA(Satisfaction, c("Staff", "Price", "Food"), quantiles = 4, plots = FALSE)
  expect_equal(rownames(res$group$coord), c("Staff", "Price", "Food"))
  data2 <- data2MatD(Satisfaction_raw, c("hotel"), c("Staff", "Price"), support = 0:10, min.n = 100)
  expect_true(get.MatD.nrows(data2) < 30)
})

test_that("simplifyD assigns masses with the midpoint rule", {
  b <- distributionD(x = c(1, 1.1, 1.35, 1.5, 1.6, 2.5, 2.51, 4.9, 6),
                     p = c(1, 2, 4, 8, 3, 5, 10, 6, 1))
  s <- simplifyD(b, support = 1:5)
  expect_equal(s@x, 1:5)
  # 1: (-Inf, 1.5] -> 1+2+4+8 ; 2: (1.5, 2.5] -> 3+5 ; 3: (2.5,3.5] -> 10 ; 5: > 4.5 -> 6+1
  expect_equal(s@p, c(15, 8, 10, 0, 7) / 40)
  # it is the closest distribution on the support
  alt <- distributionD(1:5, p = c(16, 7, 10, 0, 7))
  expect_lt(WassSqDistD(b, s), WassSqDistD(b, alt))
})

test_that("simplified barycenter lives on the variable support", {
  ex <- WD.vec.mean(Satisfaction[, "Staff"])
  ap <- WD.vec.mean.simplified(Satisfaction[, "Staff"], support = 0:10)
  expect_true(all(ap@x %in% 0:10))
  expect_equal(sum(ap@p), 1)
  expect_lt(abs(meanD(ap) - meanD(ex)), 0.5)
  M2 <- WD.vec.mean.simplified(Satisfaction[, 1:3])
  expect_s4_class(M2, "MatD")
  S <- simplifyD(Satisfaction[1:3, 1:2], support = list(0:10, c(0, 5, 10)))
  expect_equal(get.cell.MatD(S, 1, 2)@x, c(0, 5, 10))
})
