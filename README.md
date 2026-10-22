# DDAWASS <img src="https://img.shields.io/badge/lifecycle-experimental-orange.svg" alt="lifecycle: experimental" align="right"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/Airpino/DDAWASS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Airpino/DDAWASS/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/DDAWASS)](https://CRAN.R-project.org/package=DDAWASS)
<!-- badges: end -->

**Discrete Distribution Data Analysis with Wasserstein distance**

DDAWASS is an R package for the analysis of data tables whose cells are
**discrete distributions** with a finite support: for example the
distributions of the scores given on a Likert scale by groups of respondents,
of grades, of ratings or of counts. Each row is a statistical unit (a group, a
territory, a company, ...), each column a distributional variable, and
distributions are compared through the **L2 Wasserstein distance** between
their quantile functions.

It is the discrete counterpart of
[HistDAWass](https://CRAN.R-project.org/package=HistDAWass) (histogram data)
and provides:

- S4 classes `distributionD` (a discrete distribution) and `MatD` (a matrix
  of discrete distributions), built from raw microdata, frequency tables,
  histograms or legacy nested tibbles;
- exact Wasserstein distances (computed in C++), with the decomposition into
  position, size and shape components;
- barycenters (Fréchet means), Wasserstein variance, covariance and
  correlation, and a barycenter projected on the original categories of the
  scale for presentation;
- hierarchical clustering and k-means;
- quantile-based principal component analysis and multiple factor analysis;
- two-component linear regression;
- a distributional Cronbach's alpha.



## Installation

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Airpino/DDAWASS", build_vignettes = TRUE)
```

The package contains C++ code: on Windows,
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) is needed to install
it from source.

## A quick example

The `univer` data of the [CUB](https://CRAN.R-project.org/package=CUB)
package contain the 1-7 evaluations of the orientation services of the
University of Naples Federico II given by 2179 students. We describe each
faculty-by-gender group by the distributions of its scores.

```r
library(DDAWASS)
data(univer, package = "CUB")
univer$Faculty <- factor(sprintf("F%02d", univer$faculty))
univer$Gender  <- factor(univer$gender, levels = 0:1, labels = c("M", "W"))
items <- c("informat", "willingn", "officeho", "compete", "global")

M <- data2MatD(univer, unit = c("Faculty", "Gender"), vars = items,
               support = 1:7, min.n = 10)
M[1:4, 1:3]
#> a matrix of discrete distributions
#> 4 rows x 3 columns
#>       informat             willingn             officeho
#> F01_M [m=5.35 s=1.73 k=7]  [m=6.19 s=1.23 k=7]  [m=5.77 s=1.41 k=7]
#> F01_W [m=5.16 s=1.72 k=7]  [m=6.42 s=0.674 k=7] [m=5.42 s=1.57 k=7]
#> F02_M [m=6.62 s=0.785 k=7] [m=6.65 s=0.748 k=7] [m=6.3 s=1.19 k=7]
#> F02_W [m=6.47 s=1.07 k=7]  [m=6.49 s=0.965 k=7] [m=6.23 s=1.34 k=7]
```

Distance between two groups, decomposed into the contributions of the means
(`POSITION`), of the standard deviations (`SIZE`) and of the shapes:

```r
WassSqDistD(get.cell.MatD(M, "F02_M", "global"),
            get.cell.MatD(M, "F06_M", "global"), details = TRUE)
#>      SQ.W  POSITION      SIZE     SHAPE       rQQ
#> 5.4861184 4.0232281 0.9578624 0.5050279 0.7773405
```

Barycenter of the global satisfaction, expressed on the 1-7 scale, and
Wasserstein correlations between the items:

```r
WD.vec.mean.simplified(M[, "global"], support = 1:7)
#> Discrete distribution (7 support values)
#>  x       p     cdf
#>  1 0.00000 0.00000
#>  2 0.00870 0.00870
#>  3 0.04008 0.04878
#>  4 0.10507 0.15385
#>  5 0.15385 0.30769
#>  6 0.45012 0.75781
#>  7 0.24219 1.00000
#> mean = 5.72317   std = 1.11335

round(WD.correlation(M), 2)
#>          informat willingn officeho compete global
#> informat     1.00     0.64     0.56    0.78   0.75
#> willingn     0.64     1.00     0.51    0.66   0.72
#> officeho     0.56     0.51     1.00    0.53   0.58
#> compete      0.78     0.66     0.53    1.00   0.83
#> global       0.75     0.72     0.58    0.83   1.00
```

Clustering and regression:

```r
set.seed(123)
km <- WD_kmeans(M, k = 3, rep = 20)
km$solution$cardinality
#>  1  2  3
#>  4 12 10

pars <- WD.regression.two.components(M, Yvar = "global",
                                     Xvars = c("informat", "willingn", "officeho", "compete"))
round(pars, 3)
#> (AV_Intercept)    AV_informat    AV_willingn    AV_officeho     AV_compete
#>         -0.217          0.049          0.232          0.156          0.610
#>   CEN_informat   CEN_willingn   CEN_officeho    CEN_compete
#>          0.162          0.240          0.056          0.419
```

## Vignettes

1. *Building discrete distributional data*: from microdata or frequency
   tables to a `MatD`.
2. *Basic statistics*: statistics of the cells, barycenters, variance
   decomposition, covariance, correlation and reliability.
3. *Clustering*: Wasserstein distances, hierarchical clustering and k-means.
4. *Principal component analysis*: quantile PCA, Spanish-fan plot and
   multiple factor analysis.
5. *Linear regression*: the two-component model, predictions, goodness of fit
   and bootstrap.

```r
browseVignettes("DDAWASS")
```

## Main functions

| Task | Functions |
|---|---|
| Build data | `distributionD()`, `MatD()`, `data2MatD()`, `long2MatD()`, `MatH2MatD()`, `tibble2MatD()` |
| Describe | `get.MatD.stats()`, `meanD()`, `stdD()`, `skewD()`, `compQ()`, `compP()` |
| Barycenter and dispersion | `WD.vec.mean()`, `WD.vec.mean.simplified()`, `simplifyD()`, `WD.var.covar()`, `WD.correlation()`, `WD.SSQ()` |
| Distances | `WassSqDistD()`, `WD_MAT_DIST()`, `rQQ()` |
| Clustering | `WD_hclust()`, `WD_kmeans()` |
| Factorial methods | `WD.1d.PCA()`, `WD.MultiplePCA()`, `WD.plot.Spanish.fan()` |
| Regression | `WD.regression.two.components()`, `WD.regression.two.components.predict()`, `WD.regression.GOF()` |
| Reliability | `WD.cronbach()` |

## References

- Crocetta C., Irpino A., Antonucci L., Marin C. (2025). A performance
  indicator and its decomposition according to the impacts of different
  aspects based on distributional data. *Quality & Quantity*, 59(3),
  2091–2110. <https://doi.org/10.1007/s11135-024-01956-0>
- Irpino A., Verde R. (2015). Basic statistics for distributional symbolic
  variables: a new metric-based approach. *Advances in Data Analysis and
  Classification*, 9(2), 143–175. <https://doi.org/10.1007/s11634-014-0176-4>
- Irpino A., Verde R. (2015). Linear regression for numeric symbolic
  variables: a least squares approach based on Wasserstein distance.
  *Advances in Data Analysis and Classification*, 9(1), 81–106.
  <https://doi.org/10.1007/s11634-015-0197-7>

## License

GPL (>= 3)
