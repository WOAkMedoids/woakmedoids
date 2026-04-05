# WOAkMedoids

[![CRAN status](https://www.r-pkg.org/badges/version/WOAkMedoids)](https://CRAN.R-project.org/package=WOAkMedoids)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

An R package that implements the Whale Optimization Algorithm (WOA) for k-medoids clustering, providing tools for effective and efficient cluster analysis in various data sets. Supported distance measures include Dynamic Time Warping (DTW) and Euclidean Distance (ED).

## Installation

Install the released version from CRAN:
```r
install.packages("WOAkMedoids")
```

Or install the development version from GitHub:
```r
devtools::install_github("WOAkMedoids/woakmedoids")
```

## Usage
```r
library(WOAkMedoids)

data(Lightning7)
Lightning7_data <- Lightning7[, -1]  # Remove the first column of classification data

result <- woa_kmedoids(Lightning7_data, ClusNum = 7, distance_method = "dtw", learned_w = 5)
print(result)
```

## Reference

Chenan H. and Tsutsumida N. (2025) A scalable k-medoids clustering via whale optimization algorithm, *Array*, 28, 100599. [https://doi.org/10.1016/j.array.2025.100599](https://doi.org/10.1016/j.array.2025.100599)