# Extended Oscillations Function Source - Package work
# By Hannah De los Santos
# ECHO v 1.61
# Code description: Contains funcitons necessary for extended harmonic oscillator work, in order to have less confusion between scripts.

#' Function to calculate the results for all genes using the extended circadian harmonic oscillator (ECHO) method.
#'
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @param begin first time point for dataset
#' @param end last time point for dataset
#' @param resol resolution of time points
#' @param num_reps number of replicates
#' @param low lower limit when looking for rhythms, in hours. May be unused if finding rhythms of any length within timecouse (run_all_per is TRUE).
#' @param high upper limit when looking for rhythms, in hours. May be unused if finding rhythms of any length within timecouse (run_all_per is TRUE).
#' @param run_all_per boolean which indicates whether or not rhythms of any length within timecourse should be searched for.
#' @param paired if replicate data, whether the replicates are related (paired) or not (unpaired)
#' @param rem_unexpr boolean indicating whether genes with less than rem_unexpr_amt percent expression should not be considered
#' @param rem_unexpr_amt percentage of expression for which genes should not be considered if rem_unexpr is TRUE
#' @param rem_unexpr_amt_below cutoff for expression
#' @param is_normal boolean that indicates whether data should be normalized or not
#' @param is_de_linear_trend boolean that indicates whether linear trends should be removed from data or not
#' @param is_smooth boolean that indicates whether data should be smoothed or not
#' @return results, a data frame which contains:
#'   \item{Gene Name}{gene name}
#'   \item{Convergence}{did the fit converge, or descriptor of type of data (constant, unexpressed, etc.)}
#'   \item{Iterations}{number of iterations}
#'   \item{Amplitude.Change.Coefficient}{Amplitude change coefficient value for fit}
#'   \item{Oscillation Type}{Type of oscillation (damped, driven, etc.)}
#'   \item{Initial.Amplitude}{Initial amplitude value for fit}
#'   \item{Radian.Frequency}{Radian frequency for fit}
#'   \item{Period}{Period for fit (in time units)}
#'   \item{Phase Shift}{Phase shift for fit (radians)}
#'   \item{Hours Shifted}{Phase shift for fit (hours)}
#'   \item{Equilibrium Value}{Equilibrium shift for fit}
#'   \item{Slope}{Slope value of original data, if linear baseline is removed}
#'   \item{Tau}{Kendall's tau between original and fitted values}
#'   \item{P-value}{P-value calculated based on Kendall's tau}
#'   \item{BH Adj P-Value}{Benjamini-Hochberg adjusted p-values}
#'   \item{BY Adj P-Value}{Benjamini-Yekutieli adjusted p-values}
#'   \item{Original TPX.Y}{Original values for gene expression at time point X, replicate Y}
#'   \item{Fitted TPX}{Fitted values for gene expression at time point X}
#' @export
#' @importFrom stats p.adjust
#' @examples
#' # for more elaboration, please see the vignette
#' # "expressions" is the example echo.find data frame
#' \donttest{ # long example - commented out
#' echo_find(genes = expressions, begin = 2, end = 48, resol = 2,
#'   num_reps = 3, low = 20, high = 26, run_all_per = FALSE,
#'   paired = FALSE, rem_unexpr = FALSE, rem_unexpr_amt = 70, rem_unexpr_amt_below=0,
#'   is_normal = FALSE, is_de_linear_trend = FALSE, is_smooth = FALSE)
#' }
echo_find <- function(genes, begin, end, resol, num_reps, low = 1, high = 2, run_all_per, paired, rem_unexpr, rem_unexpr_amt = 70, rem_unexpr_amt_below=0, is_normal, is_de_linear_trend, is_smooth){
    # creating times sequence used for the genes
    timen <- seq(begin,end,resol) # time points for cicadian rhythms
    tied <- paired
    print("here")
    if (num_reps == 1){ # one replicate, default to true paired-ness
      tied <- TRUE
    }

  # if genes only has one row, add another constant row to alleviate automatic vectorization
  if (nrow(genes)==1){
    # creating a constant row and adding it to genes
    add_row <- data.frame(matrix(0L, 1, ncol(genes)))
    add_row[1,1] <- "not considering"
    colnames(add_row) <- colnames(genes)
    genes <- rbind(genes,add_row)
    add_one <- TRUE # marker for appropriate displays for progress bar
  } else{
    add_one <- FALSE # marker for appropriate displays for progress bar
  }

  # run all genes:

  #rem_unexpr # indicator for removing unexpressed genes
  rem_unexpr_amt <- (rem_unexpr_amt)/100 # threshold for removing unexpressed genes, converted to a decimal
  # if yes, check for genes that are unexpressed before preprocessing
  if (rem_unexpr){
    rem_unexpr_vect <- genes_unexpressed_all(rem_unexpr_amt, abs(rem_unexpr_amt_below), genes)
  } else{
    rem_unexpr_vect <- rep(FALSE,nrow(genes))
  }

  # normalize and store original data
  if (is_normal){
    norm_list <- normalize_all(genes)
    genes <- norm_list$dat
  }

  # remove baseline
  if (is_de_linear_trend){
    res_list <- de_linear_trend_all(timen,num_reps,tied,genes)
    genes <- res_list$res_df # expressions with removed baseline
    beta <- res_list$beta # slopes
  } else {
    beta <- rep(NA, nrow(genes))
  }

  # getting average data, for more than one replicate
  avg_genes <- avg_all_rep(num_reps,genes,timen)

  # smooth the data, if requested
  if (is_smooth){
    is_weighted <- TRUE # only offer weighted smoothing

    # create smoothed matrix, reassign genes
    if (tied){ # if paired replicates
      genes <- smoothing_all_tied(is_weighted, num_reps, genes,tied,timen)
    }
    else{ # if unpaired replicates
      genes <- smoothing_all_untied(is_weighted, num_reps, avg_genes, genes,timen)
    }
  }

  if (is_smooth){ # average should correspond to smoothed data
    # not unsmoothed data
    avg_genes <- avg_all_rep(num_reps,genes,timen)
  }

  # figuring out whether a range is wanted, adjusting accordingly
  if (run_all_per){ # empty low input, adjust to time series
    if (resol >= 1){
      low <- 2*pi/resol
    }
    else{ # if the begining is <1, smallest period available is 1 (because hours)
      low <- 2*pi/1
    }
    high <- 2*pi/(resol*length(timen))
  } else{ # periods are specified
    low <- 2*pi/as.numeric(low)
    high <- 2*pi/as.numeric(high)
  }

  # if more than one replicate or requested, an exact distribution is needed
  # create exact distribution for pvalues
  jtklist <- jtkdist(length(timen), reps = num_reps)
  jtk.alt <- list() # preallocate pvalue distribution for missing data

  # where we put the result
  # preallocate results matrix
  total_results <- data.frame(matrix(0,nrow = nrow(genes), ncol = 13+length(rep(timen, each = num_reps))+length(timen)))

  # renaming columns of the final results
  colnames(total_results) <- c("Gene Name","Convergence","Iterations","Amplitude.Change.Coefficient","Oscillation Type","Initial.Amplitude","Radian.Frequency","Period","Phase Shift","Hours Shifted","Equilibrium Value", "Tau", "P-Value", paste(rep("Original TP",length(rep(timen, each = num_reps))),rep(timen, each = num_reps),rep(".",length(rep(timen, each = num_reps))),rep(c(1:num_reps), length(timen)),sep=""), paste(rep("Fitted TP",length(timen)),timen,sep=""))

  # now go through one by one and get results
  for (i in 1:nrow(genes)){
    res <- calculate_param(i, timen, resol, num_reps, tied = tied, is_smooth = is_smooth, is_weighted = is_weighted,low = low,high = high,rem_unexpr = rem_unexpr, rem_unexpr_amt = rem_unexpr_amt, jtklist, genes, rem_unexpr_vect, avg_genes, jtk.alt)
    total_results[i,] <- res$results
    if (any(is.na(unlist(genes[i,-1])))){
      jtk.alt <- res$jtk.alt
    }
  }

  # remove the fake row I added if there is only one gene
  if (add_one){
    total_results <- total_results[-nrow(total_results),]
  }

  # add slope
  total_results <- cbind(total_results[,c(1:11)],`Slope` = beta, total_results[,c(12:ncol(total_results))])

  adjusted_p_val_us <- p.adjust(unlist(total_results$`P-Value`), method = "BH") # benjamini-hochberg adjust p-values
  total_results <- cbind(total_results[,c(1:14)],`BH Adj P-Value` = adjusted_p_val_us, total_results[,c(15:ncol(total_results))]) # assign to data frame

  # adding the benjamini-hochberg-yekutieli p-value adjustment
  total_results <- cbind(total_results[,c(1:15)],`BY Adj P-Value` = p.adjust(unlist(total_results$`P-Value`), method = "BY"), total_results[,c(16:ncol(total_results))])

  return(total_results)
}
