#' Lightning7 Data for Testing
#'
#' A dataset containing example data for testing purposes from the UCR Time Series Classification Archive.
#' This dataset is a time series dataset with correct classifications in the first column.
#' There are 7 classes in this dataset.
#' It contains 73 series, each with 319 time points, and the best DTW window length for this dataset is 5.
#'
#' @docType data
#' @keywords datasets
#' @name Lightning7
#' @usage data(Lightning7)
#' @format A data frame with 73 rows and 320 columns. The first column (V1) is a factor vector of correct classifications, and the remaining 319 columns (V2 to V320) are numeric vectors of time series data.
#' @source UCR Time Series Classification Archive
#' @references
#' \itemize{
#'   \item Eads, Damian R., et al. "Genetic algorithms and support vector machines for time series classification." Applications and Science of Neural Networks, Fuzzy Systems, and Evolutionary Computation V. Vol. 4787. International Society for Optics and Photonics, 2002.
#'   \item \url{http://www.timeseriesclassification.com/description.php?Dataset=Lightning7}
#' }
#' @examples
#' data(Lightning7)
#' head(Lightning7)
NULL
