// Core L2 Wasserstein computations for discrete (finite-support) distributions.
//
// A discrete distribution is described by its sorted support `x` and its
// cumulative distribution `cdf` (non-decreasing, last value equal to 1).
// Its quantile function is the step function Q(t) = x[k] for
// t in (cdf[k-1], cdf[k]].  All the quantities below are integrals of
// products of two such step functions, computed exactly by merging the two
// cdf grids in linear time (two-pointer sweep).

#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

static const double EPS_CDF = 1e-12;

// generic sweep: returns integral_0^1 g(Qa(t), Qb(t)) dt
// mode 0: (a-b)^2 ; mode 1: a*b
static double sweep(const double* xa, const double* ca, int na,
                    const double* xb, const double* cb, int nb, int mode) {
  int i = 0, j = 0;
  double t = 0.0, acc = 0.0;
  while (i < na && j < nb) {
    double nxt = std::min(ca[i], cb[j]);
    double w = nxt - t;
    if (w > 0) {
      if (mode == 0) {
        double d = xa[i] - xb[j];
        acc += w * d * d;
      } else {
        acc += w * xa[i] * xb[j];
      }
      t = nxt;
    }
    bool adv_a = ca[i] <= nxt + EPS_CDF;
    bool adv_b = cb[j] <= nxt + EPS_CDF;
    if (adv_a) i++;
    if (adv_b) j++;
  }
  return acc;
}

// Squared L2 Wasserstein distance between two discrete distributions (C++)
// @param xa,xb sorted support points
// @param ca,cb cumulative probabilities (last value = 1)
// @return the squared distance
// @keywords internal
// [[Rcpp::export(name = ".wd_sq_c")]]
double wd_sq_c(NumericVector xa, NumericVector ca,
               NumericVector xb, NumericVector cb) {
  return sweep(xa.begin(), ca.begin(), xa.size(),
               xb.begin(), cb.begin(), xb.size(), 0);
}

// Inner product of the quantile functions of two discrete distributions (C++)
// @inheritParams wd_sq_c
// @return integral of Qa(t) Qb(t) dt over (0,1)
// @keywords internal
// [[Rcpp::export(name = ".wd_dot_c")]]
double wd_dot_c(NumericVector xa, NumericVector ca,
                NumericVector xb, NumericVector cb) {
  return sweep(xa.begin(), ca.begin(), xa.size(),
               xb.begin(), cb.begin(), xb.size(), 1);
}

// helper: extract x and cdf vectors of cell (i) of column list
struct Cells {
  std::vector<NumericVector> x, c;
};

static Cells unpack(List X, List C) {
  Cells out;
  int n = X.size();
  out.x.reserve(n);
  out.c.reserve(n);
  for (int i = 0; i < n; i++) {
    out.x.push_back(as<NumericVector>(X[i]));
    out.c.push_back(as<NumericVector>(C[i]));
  }
  return out;
}

// Matrix of squared distances between rows of a MatD (C++)
// @param X list (one element per variable) of lists (one per row) of supports
// @param C list (one element per variable) of lists (one per row) of cdfs
// @param w weights of the variables
// @return a n x n matrix of weighted sums of squared distances
// @keywords internal
// [[Rcpp::export(name = ".wd_distmat_c")]]
NumericMatrix wd_distmat_c(List X, List C, NumericVector w) {
  int p = X.size();
  std::vector<Cells> cols;
  for (int j = 0; j < p; j++) cols.push_back(unpack(X[j], C[j]));
  int n = cols[0].x.size();
  NumericMatrix D(n, n);
  for (int i1 = 0; i1 < n - 1; i1++) {
    for (int i2 = i1 + 1; i2 < n; i2++) {
      double acc = 0.0;
      for (int j = 0; j < p; j++) {
        const NumericVector& xa = cols[j].x[i1];
        const NumericVector& ca = cols[j].c[i1];
        const NumericVector& xb = cols[j].x[i2];
        const NumericVector& cb = cols[j].c[i2];
        acc += w[j] * sweep(xa.begin(), ca.begin(), xa.size(),
                            xb.begin(), cb.begin(), xb.size(), 0);
      }
      D(i1, i2) = acc;
      D(i2, i1) = acc;
    }
  }
  return D;
}

// Squared distances of rows to a set of prototypes (C++)
// @param X,C data supports and cdfs (see .wd_distmat_c)
// @param XP,CP prototypes supports and cdfs, same layout
// @param w weights of the variables
// @return a n x k matrix
// @keywords internal
// [[Rcpp::export(name = ".wd_dist_to_protos_c")]]
NumericMatrix wd_dist_to_protos_c(List X, List C, List XP, List CP,
                                  NumericVector w) {
  int p = X.size();
  std::vector<Cells> cols, prot;
  for (int j = 0; j < p; j++) {
    cols.push_back(unpack(X[j], C[j]));
    prot.push_back(unpack(XP[j], CP[j]));
  }
  int n = cols[0].x.size();
  int k = prot[0].x.size();
  NumericMatrix D(n, k);
  for (int i = 0; i < n; i++) {
    for (int h = 0; h < k; h++) {
      double acc = 0.0;
      for (int j = 0; j < p; j++) {
        const NumericVector& xa = cols[j].x[i];
        const NumericVector& ca = cols[j].c[i];
        const NumericVector& xb = prot[j].x[h];
        const NumericVector& cb = prot[j].c[h];
        acc += w[j] * sweep(xa.begin(), ca.begin(), xa.size(),
                            xb.begin(), cb.begin(), xb.size(), 0);
      }
      D(i, h) = acc;
    }
  }
  return D;
}

// Gram matrix of quantile inner products between two lists of distributions
// @param XA,CA first list of supports/cdfs
// @param XB,CB second list of supports/cdfs
// @return matrix G[a,b] = integral Qa Qb
// @keywords internal
// [[Rcpp::export(name = ".wd_gram_c")]]
NumericMatrix wd_gram_c(List XA, List CA, List XB, List CB) {
  Cells A = unpack(XA, CA), B = unpack(XB, CB);
  int na = A.x.size(), nb = B.x.size();
  NumericMatrix G(na, nb);
  for (int a = 0; a < na; a++)
    for (int b = 0; b < nb; b++)
      G(a, b) = sweep(A.x[a].begin(), A.c[a].begin(), A.x[a].size(),
                      B.x[b].begin(), B.c[b].begin(), B.x[b].size(), 1);
  return G;
}

// Weighted Wasserstein barycenter (Frechet mean) of discrete distributions
// @param X list of supports
// @param C list of cdfs
// @param w weights (summing to one)
// @return list(x, cdf) of the barycenter (support not yet compacted)
// @keywords internal
// [[Rcpp::export(name = ".wd_barycenter_c")]]
List wd_barycenter_c(List X, List C, NumericVector w) {
  Cells D = unpack(X, C);
  int n = D.x.size();
  // merged grid of cdf levels
  std::vector<double> grid;
  for (int i = 0; i < n; i++)
    for (int k = 0; k < D.c[i].size(); k++) grid.push_back(D.c[i][k]);
  std::sort(grid.begin(), grid.end());
  std::vector<double> g;
  for (size_t k = 0; k < grid.size(); k++) {
    if (g.empty() || grid[k] > g.back() + EPS_CDF) g.push_back(grid[k]);
  }
  if (g.empty() || g.back() < 1.0 - EPS_CDF) g.push_back(1.0);
  g.back() = 1.0;
  int m = g.size();
  NumericVector bx(m), bc(m);
  for (int k = 0; k < m; k++) bc[k] = g[k];
  for (int i = 0; i < n; i++) {
    const NumericVector& xi = D.x[i];
    const NumericVector& ci = D.c[i];
    int pos = 0, len = xi.size();
    for (int k = 0; k < m; k++) {
      // quantile at level just below g[k]: first pos with ci[pos] >= g[k]
      while (pos < len - 1 && ci[pos] < g[k] - EPS_CDF) pos++;
      bx[k] += w[i] * xi[pos];
    }
  }
  return List::create(_["x"] = bx, _["cdf"] = bc);
}

// Quantiles of a discrete distribution (C++)
// @param x support, @param cdf cumulative probs, @param p levels
// @return Q(p) with the convention Q(0) = min(x)
// @keywords internal
// [[Rcpp::export(name = ".wd_quantile_c")]]
NumericVector wd_quantile_c(NumericVector x, NumericVector cdf,
                            NumericVector p) {
  int n = p.size(), len = x.size();
  NumericVector q(n);
  for (int k = 0; k < n; k++) {
    int lo = 0, hi = len - 1;
    // first index with cdf >= p (binary search)
    while (lo < hi) {
      int mid = (lo + hi) / 2;
      if (cdf[mid] >= p[k] - EPS_CDF) hi = mid; else lo = mid + 1;
    }
    q[k] = x[lo];
  }
  return q;
}
