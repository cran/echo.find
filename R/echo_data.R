#' Synthetic expression data for 12 genes.
#'
#' A dataset containing the names and expression values for 12 synthetically
#' generated samples. This example data has time points from 2 to 48 hours with
#' 2 hour resolution and 3 replicates. Random missing data is also included.
#' Synthetic data was created by randomly selecting parameters for the
#' extended harmonic oscillator equation (see journal paper link in vignette
#' for the equation), then adding random uniform noise to each expression.
#'
#' Note the data format: its first column first column has gene labels/names, and
#' all other columns have expression data. This expression data is ordered by
#' time point then by replicate, and has evenly spaced time points. Any missing
#' data has cells left blank.
#'
#' @format A data frame with 12 rows and 73 variables (column 1: sample labels, columns to 2 to 73: numerical values for gene expression in the forsmat CTX.Y (time point X, replicate Y)).
#'
"expressions"
