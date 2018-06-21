# Extended Oscillations Function Source - Package work
# By Hannah De los Santos
# ECHO v 1.61
# Code description: Contains funcitons necessary for extended harmonic oscillator work, in order to have less confusion between scripts - unecessary.

#' Formula for damped oscillator with phase and equilibrium shift.
#'
#' @param a Amplitude
#' @param gam Forcing coefficient (amount of damping/driving)
#' @param omega Radial frequency
#' @param phi Phase Shift (radians)
#' @param y_shift Equilibrium shift
#' @param t time
#' @return result of inputs into formula
#' @keywords internal
#' @noRd
alt_form <- function(a,gam,omega,phi,y_shift,t){
  return (a*exp(-1*gam*(t)/2)*cos((omega*t)+phi)+y_shift)
}

#' Function to calculate the average of replicates for all given expressions. Used primarily for cases with multiple replicates.
#'
#' @param num_reps number of replicates
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @param timen vector of time points for gene expression
#' @return a matrix of means of gene expressions for replicates
#' @keywords internal
#' @noRd
avg_all_rep <- function(num_reps,genes,timen){
  # originally based on heat map code, but it will work fine here

  #get matrix of just the relative expression over time
  hm_mat <- as.matrix(genes[,2:ncol(genes)])

  #if there are replicates, average the relative expression for each replicate
  mtx_reps <- list() # to store actual replicate matrix
  mtx_count <- list() # to store how many are NA
  for (i in 1:num_reps){
    mtx_reps[[i]] <- hm_mat[, seq(i,ncol(hm_mat), by=num_reps)]
    mtx_count[[i]] <- is.na(mtx_reps[[i]])
    mtx_reps[[i]][is.na(mtx_reps[[i]])] <- 0
  }
  repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat))+num_reps # to store how many we should divide by
  hm_mat <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat)) # to store the final result
  for (i in 1:num_reps){
    hm_mat <- hm_mat + mtx_reps[[i]] # sum the replicates
    repmtx <- repmtx - mtx_count[[i]] # how many replicates are available for each time point
  }
  repmtx[repmtx==0] <- NA # to avoid division by 0 and induce NAs if there are no time points available
  hm_mat <- hm_mat/repmtx

  return(hm_mat)
}

#' Function to calculate the variance of replicates at a certain time point.
#'
#' @param x time point
#' @param current_gene row number of gene being examined
#' @param num_reps number of replicates
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @return the variance of replicates at a certain time point
#' @importFrom stats var
#' @keywords internal
#' @noRd
calc_var <- function(x, current_gene, num_reps, genes){
  eps <- 1e-7 # slight adjustment for 0 variance case
  std2 <- var(unlist(genes[current_gene,c(x:(num_reps-1+x))]), na.rm = TRUE) # calc variance
  return(1/(std2+eps))
}

#' Function to calculate the weights for replicate fitting (these weights are the inverse of the variance at each time point).
#'
#' @param current_gene row number of gene being examined
#' @param num_reps number of replicates
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @return the variance of replicates at a certain time point
#' @keywords internal
#' @noRd
calc_weights <- function(current_gene, num_reps, genes){
  return(sapply(seq(2,ncol(genes), by = num_reps), function(x) calc_var(x,current_gene,num_reps, genes)))
}

#' Shewchuk algorithms for adaptive precision summation used in jtkdist
#' http://www.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
#' @param a first number to sum
#' @param b first number to sum
#' @return two different sums
#' @author Jonathan Shewchuk
#' @keywords internal
#' @noRd
fast.two.sum <- function(a,b) { # known abs(a) >= abs(b)
  x <- a+b
  bv <- x-a
  y <- b-bv
  if(y==0) return(x)
  c(x,y)
}

#' Shewchuk algorithms for adaptive precision summation used in jtkdist
#' http://www.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
#' @param a first number to sum
#' @param b first number to sum
#' @return two different sums
#' @author Jonathan Shewchuk
#' @keywords internal
#' @noRd
two.sum <- function(a,b) { # unknown order
  x <- a+b
  bv <- x-a
  av <- x-bv
  br <- b-bv
  ar <- a-av
  y <- ar+br
  if(y==0) return(x)
  c(x,y)
}

#' Shewchuk algorithms for adaptive precision summation used in jtkdist
#' http://www.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
#' @param g sum to expand
#' @return two pieces of expanded sum
#' @author Jonathan Shewchuk
#' @keywords internal
#' @noRd
expansion.sum <- function(g) {
  g <- g[order(abs(g))]
  z <- fast.two.sum(g[2],g[1])
  q <- z[1]
  h <- NULL
  if(length(z)!=1) h <- z[2]
  n <- length(g)
  if(n==2) return(c(h,q))
  for(i in 3:n) {
    z <- two.sum(q,g[i])
    q <- z[1]
    if(length(z)!=1) h <- c(h,z[2])
  }
  c(h,q) # strongly non-overlapping values
}

#'  details: http://www.jstor.org/pss/2347656
#'
#' @param  timepoints number of total amount of timepoints
#' @param  reps number of replicates
#' @param  normal boolean indicating whether a normal distribution is desired
#' @param  alt boolean indicating whether there are missing time points
#' @param jtk.alt contains the exact p-value distribution for replicate data with missing data
#' @return exact distribution of JTK, whether that be based on missing values or not
#' @author E.F. Harding
#' @keywords internal
#' @noRd
jtkdist <- function(timepoints,reps=1,normal=FALSE,alt=FALSE, jtk.alt) {


  if(length(reps)==timepoints) {
    tim <- reps # support for unbalanced replication
  } else {
    tim <- rep(reps[1],timepoints) # balanced replication
  }

  maxnlp <- lfactorial(sum(tim))-sum(lfactorial(tim)) #maximum possible negative log p-value
  limit <- log(.Machine$double.xmax) #largest representable nlp
  normal <- normal | (maxnlp>limit-1) #switch to normal approximation if maxnlp is too large

  if(alt) { # if there are missing time points
    lab <- paste(sort(tim), collapse=",") # create label
    if(lab %in% names(jtk.alt)) { # if this distribution has already been created
      alt.id <- match(lab, names(jtk.alt))
      jtkaltlist <- list("jtk.alt"=jtk.alt, "alt.id"=length(jtk.alt))
      return(jtkaltlist) # no need to create a new distribution, use the the old label
    }
    nn <- sum(tim) # number of time points
    M <- (nn^2-sum(tim^2))/2 # number of combinations
    jtk.alt[[lab]] <- list() # save in the global space for efficiency
    jtk.alt[[lab]]$MAX <- M
    if(normal) { # create normal distribution is desired or necessary
      var <- (nn^2*(2*nn+3) -
                sum(tim^2*(2*tim+3)))/72 # variance
      jtk.alt[[lab]]$SDV <- sqrt(var) # standard deviation
      jtk.alt[[lab]]$EXV <- M/2 # expected value2
      return(length(jtk.alt))
    }
  } else { # exact distribution
    jtk.grp.size <- tim # sizes of each replicate group
    jtk.num.grps <- length(tim) # timepoints = number of groups
    jtk.num.vals <- nn <- sum(tim) # number of data values (independent of period and lag)
    jtk.max <- M <- (nn^2-sum(tim^2))/2 # maximum possible jtk statistic
    jtk.grps <- rep(1:length(tim), ti=tim) # group labels
    jtk.dims <- c(nn*(nn-1)/2,1) # dimensions for distribution

    if(normal) { # unused, practically
      JTK.VAR <- (nn^2*(2*nn+3) -
                    sum(tim^2*(2*tim+3)))/72 # variance of jtk
      JTK.SDV <- sqrt(JTK.VAR) # standard deviation of jtk
      JTK.EXV <- jtk.max/2 # expected value of jtk
      JTK.EXACT <- FALSE
      return(invisible(0)) # omit calculation of exact distribution
    }
  }
  MM <- floor(M/2) # mode of this possibly alternative jtk distribution
  cf <- as.list(rep(1,MM+1)) # initial lower half cumulative frequency distribution

  size <- tim # sizes of each group of known replicate values
  size <- size[order(size)]  # ascending order for fastest calculation
  k <- length(tim) # number of groups of known replicate values

  # start building distribution
  N <- size[k]
  if(k>2) for(i in (k-1):2) {
    N <- c(size[i]+N[1],N)
  }
  for(i in 1:(k-1)) { # count permutations using the Harding algorithm
    m <- size[i]
    n <- N[i]

    if(n < MM) {
      P <- min(m+n,MM)
      for(t in (n+1):P) { # zero-based offset t
        for(u in 1+MM:t) { # one-based descending index u
          cf[[u]] <- expansion.sum( # Shewchuck algorithm
            c(cf[[u]],-cf[[u-t]]))
        }
      }
    }
    Q <- min(m,MM)
    for(s in 1:Q) { # zero-based offset s
      for(u in 1+s:MM) { # one-based ascending index u
        cf[[u]] <- expansion.sum( # Shewchuck algorithm
          c(cf[[u]],cf[[u-s]]))
      }
    }
  }
  cf <- sapply(cf,sum)

  # cf now contains the lower-half cumulative frequency distribution;
  # append the symmetric upper-half cumulative distribution to cf

  if(M %% 2) {
    cf <- c(cf,2*cf[MM+1]-c(cf[MM:1],0)) # if M is odd (mode is duplicated)
  } else {
    cf <- c(cf,cf[MM+1]+cf[MM]-c(cf[MM:2-1],0)) # if M is even (unique mode is in lower half)
  }
  jtkcf <- rev(cf) # upper-tail cumulative frequencies for all integer jtk
  ajtkcf <- (jtkcf[-length(cf)]+jtkcf[-1])/2 # interpolated cumulative frequency values for all half-integer jtk

  id <- 1+0:(2*M) # one-based indices for all jtk values
  cf <- id # container for the jtk frequency distribution
  cf[!!id%%2] <- jtkcf # odd indices for integer jtk
  cf[!id%%2] <- ajtkcf # even indices for half-integer jtk
  cp <- cf/jtkcf[1] # all upper-tail p-values

  if(alt) { # if missing values
    jtk.alt[[lab]]$CP <- cp
    jtkaltlist <- list("jtk.alt"=jtk.alt, "alt.id"=length(jtk.alt))
    return(jtkaltlist)
  }
  jtk.cp <- cp
  jtk.exact <- TRUE

  # list of relevant statistics
  jtklist <- list("jtk.grp.size"=jtk.grp.size,"jtk.num.grps"=jtk.num.grps,"jtk.num.vals"=jtk.num.vals,"jtk.max"=jtk.max,"jtk.grps"=jtk.grps,"jtk.dims"=jtk.dims,"jtk.cp"=jtk.cp,"jtk.exact"=jtk.exact)

  return(jtklist)
}

#' Function to calculate pvalues for replicates using the exact JTK distribution.
#'
#' @param ref_waveform answer to damped oscillator
#' @param times vector of time points
#' @param num_reps number of replicates
#' @param current_gene row number of gene being examined
#' @param jtklist exact distribution for JTK
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @param jtk.alt contains the exact p-value distribution for replicate data with missing data
#' @return list containing:
#' \item{pval}{pvalue}
#' \item{tau}{Kendall's tau}
#' \item{S}{Kendall's S statistic}
#' @keywords internal
#' @noRd
calc_tau_and_p_jtk <- function(ref_waveform,times,num_reps,current_gene,jtklist, genes, jtk.alt){
  # forming the correlation matrix for the experimental values
  rep_numerator <- sign(outer(unlist(genes[current_gene,-1]),unlist(genes[current_gene,-1]),"-"))
  rep_numerator <- rep_numerator[lower.tri(rep_numerator)]

  # forming the correlation matrix for the fitted values
  ref_waveform_rep <- rep(ref_waveform,each=num_reps)
  ref_numerator <- sign(outer(ref_waveform_rep,ref_waveform_rep,"-"))
  ref_numerator <- ref_numerator[lower.tri(ref_numerator)]
  S <- sum(rep_numerator*ref_numerator,na.rm = TRUE) # Kendall's S statistic

  tim <- rep(times, each = num_reps)[!is.na(genes[current_gene,-1])] # time points examined (not including points with missing values)
  tim <- as.integer(table(tim)) # getting the number of occurences for each time point, in vector form
  nn <- sum(tim) # sum number of occurences

  alt <- any(is.na(unlist(genes[current_gene,-1]))) # are there any missing values?
  if(alt) { # if missing values
    tab <- table(jtklist$jtk.grps[is.finite(unlist(genes[current_gene,-1]))])
    jtkaltlist <- jtkdist(length(tab),
                          as.integer(tab),
                          normal=FALSE,
                          alt=alt,
                          jtk.alt = jtk.alt)
    alt.id <- jtkaltlist$alt.id
    jtk.alt <- jtkaltlist$jtk.alt
  }
  M <- switch(1+alt, # maximum possible S score for this distribution
              jtklist$jtk.max, # no missing values
              jtk.alt[[alt.id]]$MAX) # if missing values

  if(!S){return (list("pval"=1,"S"=0,"tau"=0))} # if not calculated, must have a pvalue of 1
  jtk <- (abs(S)+M)/2 # two-tailed jtk statistic

  if(jtklist$jtk.exact) { # if we are using the exact distribution
    jtki <- 1+2*jtk # index into the exact upper-tail distribution
    p <- switch(1+alt, # deciding between missing values (2) or not (1)
                2*jtklist$jtk.cp[jtki],
                2*jtk.alt[[alt.id]]$CP[jtki])
  }
  # if we are using the normal distribution (unused, practically)
  # else {
  #   p <- switch(1+alt,
  #               2*pnorm(-(jtk-1/2),
  #                       -JTK.EXV,JTK.SDV),
  #               2*pnorm(-(jtk-1/2),
  #                       -JTK.ALT[[alt.id]]$EXV,
  #                       JTK.ALT[[alt.id]]$SDV))
  # }

  # include tau = S/M for this distribution
  if (alt){ # if missing value, return the alternate distribution
    return(list("pval"=p,"S"=S,"tau"=S/M, "jtk.alt"=jtk.alt))
  }
  # if no missing values, we're relying on the larger complete distribution
  return(list("pval"=p,"S"=S,"tau"=S/M))
}

#' Function to calculate the parameters for the extended harmonic oscillator equation for a specific gene.
#'
#' @param current_gene row number of current gene we want to calculate parameters for
#' @param times time points for dataset
#' @param resol resolution of time points
#' @param num_reps number of replicates
#' @param tied if replicate data, whether the replicates are related (paired) or not (unpaired)
#' @param is_smooth boolean that indicates whether data should be smoothed or not
#' @param is_weighted if there is smoothing, is it weighted (1,2,1) smoothing, or unweighted smoothing (1,1,1)
#' @param low the highest frequency we are looking for, in radians (lowest period)
#' @param high the lowest frequency we are looking for, in radians (highest period)
#' @param rem_unexpr boolean indicating whether genes with less than rem_unexpr_amt percent expression should not be considered
#' @param rem_unexpr_amt percentage of expression for which genes should not be considered
#' @param jtklist contains the exact p-value distribution for replicate data
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @param rem_unexpr_vect vector of booleans indicating whether genes with less than rem_unexpr_amt percent expression should not be considered
#' @param avg_genes matrix of average expressions over all replicates for each gene
#' @param jtk.alt contains the exact p-value distribution for replicate data with missing data
#' @return res, list which contains:
#'   \item{results}{a data frame which contains:}
#'      \item{gene}{gene name}
#'      \item{conv}{did the fit converge, or descriptor of type of data (constant, unexpressed, etc.)}
#'      \item{iter}{number of iterations}
#'      \item{gamma}{forcing coefficient value for fit}
#'      \item{type_gam}{Type of oscillation (damped, driven, etc.)}
#'      \item{amplitude}{Amplitude value for fit}
#'      \item{omega}{Radial frequency for fit}
#'      \item{period}{Period for fit (in time units)}
#'      \item{phase.shift}{Phase shift for fit (radians)}
#'      \item{hours.shift}{Phase shift for fit (hours)}
#'      \item{tau}{Kendall's tau between original and fitted values}
#'      \item{y_shift}{Equilibrium shift for fit}
#'      \item{pval}{P-value calculated based on Kendall's tau}
#'      \item{original.values}{original values for gene}
#'      \item{fitted.values}{fitted values for gene}
#'  \item{jtk.alt}{contains the exact p-value distribution for replicate data with missing data}
#' @import minpack.lm
#' @importFrom stats coef
#' @keywords internal
#' @noRd
calculate_param <- function(current_gene,times,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,jtklist=list(),genes,rem_unexpr_vect, avg_genes, jtk.alt){

  gene_n <- as.character(genes[current_gene,1]) # gene name
  # first we need to check whether or not the gene is just a straight line
  if (!is_deviating(current_gene,genes)){ # one replicate
    if (num_reps == 1){
      results <- data.frame(gene = gene_n, conv = "No Deviation", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
      return (results)
    }
    else{ # multiple replicates
      results <- data.frame(gene = gene_n, conv = "No Deviation", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
      return (results)
    }
  }

  # then we need to check if 70% are expressed (if desired)
  if (rem_unexpr){
    if (rem_unexpr_vect[current_gene]){
      if (num_reps == 1){ # one replicate

        results <- data.frame(gene = gene_n, conv = "Unexpressed", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
        return (results)

      } else{ # multiple replicates

        results <- data.frame(gene = gene_n, conv = "Unexpressed", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
        return (results)
      }
    }
  }

  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- rbind(as.numeric(as.character(t(genes[current_gene,c(2:ncol(genes))])))) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }

  tryCatch({ # throw exception upon error
    # calculate the amount of peaks
    peaks <- c(); # vector of peak values
    peaks_time <- c(); # vector of peak times
    counting <- 1; # counter
    if (resol <= 1){
      # go through gene values and find maximum as compared to 8 surrounding values
      # finding peaks for first 4 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[1:9], na.rm = TRUE)) != -Inf){
      #   for (i in 1:4){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[1:9], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }

      for(i in 5:(length(y_val)-4)){
        # deal with complete missingness
        if (suppressWarnings(max(y_val[i-4],y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3],y_val[i+4], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
          next
        }
        # otherwise continue as normal
        if (y_val[i] == max(y_val[i-4],y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3],y_val[i+4], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }

      # finding peaks for last 4 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[(length(y_val)-8):length(y_val)], na.rm = TRUE)) != -Inf){
      #   for (i in (length(y_val)-3):length(y_val)){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[(length(y_val)-8):length(y_val)], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
    } else if (resol <=2){
      # go through gene values and find maximum as compared to six surrounding values
      # finding peaks for first 3 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[1:7], na.rm = TRUE)) != -Inf){
      #   for (i in 1:3){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[1:7], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
      for(i in 4:(length(y_val)-3)){
        # deal with complete missingness
        if (suppressWarnings(max(y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
          next
        }
        # otherwise continue as normal
        if (y_val[i] == max(y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
      # finding peaks for last 3 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[(length(y_val)-6):length(y_val)], na.rm = TRUE)) != -Inf){
      #   for (i in (length(y_val)-2):length(y_val)){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[(length(y_val)-6):length(y_val)], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
    } else if (resol <= 4){
      # finding peaks for first 2 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[1:5], na.rm = TRUE)) != -Inf){
      #   for (i in 1:2){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[1:5], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
      # go through gene values and find maximum as compared to four surrounding values
      for(i in 3:(length(y_val)-2)){
        # to deal with complete missingness
        if(suppressWarnings(max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE))== -Inf  | is.na(y_val[i])){
          next
        }
        if (y_val[i] == max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }

      # finding peaks for last 3 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[(length(y_val)-4):length(y_val)], na.rm = TRUE)) != -Inf){
      #   for (i in (length(y_val)-1):length(y_val)){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[(length(y_val)-4):length(y_val)], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }

    } else{
      # finding peaks for first point
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[1:3], na.rm = TRUE)) != -Inf){
      #   for (i in 1){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[1:3], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }

      # go through gene values and find maximum as compared to two surrounding values
      for(i in 2:(length(y_val)-1)){
        # to deal with complete missingness
        if(suppressWarnings(max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE))== -Inf  | is.na(y_val[i])){
          next
        }
        if (y_val[i] == max(y_val[i-1],y_val[i],y_val[i+1], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }

      # finding peaks for last 3 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[(length(y_val)-2):length(y_val)], na.rm = TRUE)) != -Inf){
      #   for (i in length(y_val)){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[(length(y_val)-2):length(y_val)], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
    }

    # calculate starting amplitude, y_shift
    y0 <- mean(y_val,na.rm = TRUE) # intial value for the equilibrium shift
    if (y0 < 10^-10 && y0 > -10^-10){
      y0 <- 10^-8 # avoiding 0 mean, which makes gradient singular
    }
    x0 <- min(times) # the x start parameter
    a0 <- max(y_val,na.rm = TRUE) - y0 # initial guess for amplitude

    # intial value for gamma
    if (length(peaks)==0){ # if there are no peaks, we account for that
      gam0 <- 0
    } else if (which.max(peaks)==1){ # if the highest peak is first, then damping is likely
      if (length(peaks)>1){
        # trying to figure out gamma based on logarithmic decrement
        n <- peaks_time[2]-peaks_time[1]
        log_dec <- (1/n)*log(abs(peaks[1]/peaks[2]))
        gam0 <- 1/(sqrt(1+((2*pi/log_dec)^2)))
      } else{ # use a guess
        gam0 <- .01
      }
    } else{ # otherwise driving is likely
      if (length(peaks)>1){
        # trying to figure out gamma based on logarithmic decrement
        n <- peaks_time[2]-peaks_time[1]
        log_dec <- (1/n)*log(abs(peaks[2]/peaks[1]))
        gam0 <- -1*1/(sqrt(1+((2*pi/log_dec)^2)))
      } else{ # use a guess
        gam0 <- -.01
      }
    }

    # let frequency depend on amount of peaks = (length(times)*resol/(no of peaks+1 [accounts for phase shift])
    if (length(peaks) == 0){
      if (high == -Inf || low == Inf){
        w0 <- 2*pi/(length(times)*resol/2)
      } else{
        # want to get their actual integer period values
        highfix <- (high/2/pi)^-1
        lowfix <- (low/2/pi)^-1
        w0 <- 2*pi/(length(times)*resol/((highfix+lowfix)/2))
      }
    } else if (length(peaks) == 1){ # phase shift cases only one peak to appear
      w0 <- 2*pi/(length(times)*resol/(length(peaks)+1))
    } else{
      w0 <- 2*pi/(length(times)*resol/(length(peaks)))
    }

    # can't be outside the specified parameters
    if (w0 > low){
      w0 <- low
    } else if (w0 < high){
      w0 <- high
    }


    # we estimate our phase shift on the second and third nonmissing points for accuracy
    # if you have less than 3 points nonmissing, I have no hope for you
    second <- which(!is.na(y_val))[2]
    third <- which(!is.na(y_val))[3]
    min_i <- 0;
    for (i in 0:11){
      # if the phase shift at both the second and third gene values are the smallest value available for the fitted value
      if ((abs(alt_form(a0,gam0,w0,(i*pi/6),y0,times[second])-y_val[second]))<(abs(alt_form(a0,gam0,w0,(min_i*pi/6),y0,times[second])-y_val[second]))){
        if ((abs(alt_form(a0,gam0,w0,(i*pi/6),y0,times[third])-y_val[third]))<(abs(alt_form(a0,gam0,w0,(min_i*pi/6),y0,times[third])-y_val[third]))){
          min_i <- i
        }
      }
    }
    phi0 <- min_i*pi/6 # intial value for phase shift

    if (num_reps == 1){ # one replicate
      # put the times into a data frame
      temp <- data.frame(y=t(y_val),t=times)
      temp <- temp[!is.na(temp$y),] # remove any missing data points

      # fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0),
                                               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                               upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                               control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                        ftol=1e-6, ptol=1e-6, gtol=1e-6)))
    } else{ # multiple replicates
      #put the times and data point into a data frame
      weights <- calc_weights(current_gene,num_reps, genes)
      temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),
                         t=cbind(rep(times,each = num_reps)),
                         w=cbind(rep(weights,each = num_reps)))
      temp <- temp[!is.na(temp$y),] # remove any missing data points

      #fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0),
                                               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                               upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                               control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                        ftol=1e-6, ptol=1e-6, gtol=1e-6),
                                               weights = temp$w))
    }

    did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
    num_iter <- oscillator.fit$convInfo$finIter # amount of iterations

    parameters <- coef(oscillator.fit) #extract parameter estimates

    # alt_form parameters:
    # the parameters go in the order of: gam,a,omega,phi,y_shift
    gam <- parameters[1]
    a <- parameters[2]
    omega <- parameters[3]
    phi <- parameters[4]
    y_shift <- parameters[5]

    # calculating whether (over)damped, (over)driven, harmonic
    if (gam < -.15){
      type_gam <- "Overexpressed"
    } else if (gam <= -.01){
      type_gam <- "Driven"
    } else if (gam <= .01){
      type_gam <- "Harmonic"
    } else if (gam <= .15){
      type_gam <- "Damped"
    } else{
      type_gam <- "Repressed"
    }

    # calculating the phase shift in terms of period (omega inverse of period)
    if (phi > 0){ # shift to the left
      frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
      dist_peak <- frac_part*(2*pi/omega) # distance from first peak
      phase_hours <- (2*pi/omega)-dist_peak
    } else { # shift to the right
      frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
      dist_peak <- frac_part*(2*pi/omega) # distance from first peak
      phase_hours <- (2*pi/omega)+abs(dist_peak)
      phase_hours <- abs(dist_peak)
    }
    # should be no negative shifts
    # if (phase_hours < 0){ # only output positive shifts
    #   phase_hours <- phase_hours + (2*pi/omega)
    # }

    # calculate p-value
    ref_wave <- (alt_form(a,gam,omega,phi,y_shift,times)) # fitted values
    taulist <- calc_tau_and_p_jtk(ref_wave,times,num_reps,current_gene,jtklist,genes, jtk.alt)
    tau <- taulist$tau # Kendall's tau
    pval <- taulist$pval # pvalue
    # if there are any missing genes, we want to add the new distribution
    if (any(is.na(unlist(genes[current_gene,-1])))){
      jtk.alt <- taulist$jtk.alt
    }

    # list of parameters and other resulting values
    if (num_reps == 1){
      results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
      results <- cbind(results, y_val, rbind(ref_wave))
    } else{
      results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
    }
    res <- list("results"=results)
    if (any(is.na(unlist(genes[current_gene,-1])))){
      res$jtk.alt <- jtk.alt
    }
    return (res)
  }, error = function(e){ # if there's failure in convergence

    if (num_reps == 1){
      results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
      results <- cbind(results, y_val, rbind(rep(NA,length(times))))
    } else{
      results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
    }

    res <- list("results"=results)
    if (any(is.na(unlist(genes[current_gene,-1])))){
      res$jtk.alt <- jtk.alt
    }
    return (results)
  })
}

#' Function for determining whether gene values are unexpressed (less than rem_unexpr_amt percent expressed (i.e., not 0)) for full matrix.
#'
#' @param rem_unexpr_amt percentage of expression for which genes should not be considered
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @return boolean if there is rem_unexpr_amt percent expression
#' @keywords internal
#' @noRd
genes_unexpressed_all <- function(rem_unexpr_amt, genes){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])
  # get how many genes are expressed for each gene
  tot_expressed <- rowSums(all_reps != 0,na.rm = TRUE)

  # return false if amount is less than threshold
  return(tot_expressed <= (ncol(all_reps)*rem_unexpr_amt))
}


#' Function for determining whether gene values are deviating or constant.
#' @param current_gene row number of current gene we want to calculate parameters for
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @return boolean if there is deviation within the gene values
#' @importFrom stats sd
#' @keywords internal
#' @noRd
is_deviating <- function(current_gene,genes){
  stdev <- sd(genes[current_gene,-1], na.rm = TRUE)
  return (stdev != 0)
}

#' Function to calculate a rolling average (3-wide window) for a set of data (smoothing) for paired (tied) replicates, which smooths each replicate separately.
#'
#' @param is_weighted boolean if there is smoothing, is it weighted (1,2,1) smoothing,  or unweighted smoothing (1,1,1)
#' @param num_reps number of replicates
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @param tied whether replicates are paired (tied) or unpaired (untied)
#' @param timen vector of time points for gene expression
#' @return smoothed expression data
#' @keywords internal
#' @noRd
smoothing_all_tied <- function(is_weighted, num_reps, genes, tied, timen){
  # originally based on heat map code, but it will work fine here

  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])

  if (!is_weighted){

    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are NA
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]
      mtx_count[[i]] <- is.na(center_reps[[i]])
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }

    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+3 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 2 at edges
  } else{ # weighted averaging

    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are NA
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]*2
      mtx_count[[i]] <- is.na(center_reps[[i]])*2
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }

    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+4 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges

  }
  for (i in 1:num_reps){
    # sum the replicates
    left <- cbind(matrix(0,nrow(all_reps),1),center_reps[[i]][,-ncol(center_reps[[i]])]/2) # left shifted matrix
    right <- cbind(center_reps[[i]][,-1]/2,matrix(0,nrow(genes),1)) # right shifted matrix
    center_reps[[i]] <- left + center_reps[[i]] + right

    # figure out how many replicates are actually available for each time point
    left_na <- cbind(matrix(0,nrow(all_reps),1),mtx_count[[i]][,-ncol(mtx_count[[i]])]) # left shifted matrix
    right_na <- cbind(mtx_count[[i]][,-1],matrix(0,nrow(genes),1)) # right shifted matrix
    repmtx_l[[i]] <- repmtx - left_na - mtx_count[[i]] - right_na
    # to avoid division by 0 and induce NAs if there are no time points available
    repmtx_l[[i]][repmtx_l[[i]]==0] <- NA
  }

  dat <- genes
  for (x in 0:(num_reps-1)){
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  }
  dat[is.na(genes)] <- NA # do not impute missing values
  return(dat)
}

#' Function to calculate a rolling average (3-wide window) for a set of data (smoothing) for unpaired (untied) replicates, which smooths each replicate based on the average of all replicates.
#'
#' @param  is_weighted boolean if there is smoothing, is it weighted (1,2,1) smoothing, or unweighted smoothing (1,1,1)
#' @param  num_reps number of replicates
#' @param avg_genes matrix of average expressions over all replicates for each gene
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @param timen vector of time points for gene expression
#' @return smoothed expression data
#' @keywords internal
#' @noRd
smoothing_all_untied <- function(is_weighted, num_reps, avg_genes, genes,timen){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])

  if (!is_weighted){

    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are

    side_reps <- avg_genes
    mtx_side_count <- is.na(side_reps)
    side_reps[is.na(side_reps)] <- 0
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]
      mtx_count[[i]] <- is.na(center_reps[[i]])
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }

    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+3 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges
  } else{ # weighted averaging

    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are

    side_reps <- avg_genes
    mtx_side_count <- is.na(side_reps)
    side_reps[is.na(side_reps)] <- 0

    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]*2
      mtx_count[[i]] <- is.na(center_reps[[i]])*2
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }

    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+4 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges

  }
  for (i in 1:num_reps){
    # sum the replicates
    left <- cbind(matrix(0,nrow(genes),1),side_reps[,-ncol(side_reps)]) # left shifted matrix
    right <- cbind(side_reps[,-1],matrix(0,nrow(genes),1)) # right shifted matrix
    center_reps[[i]] <- left + center_reps[[i]] + right

    # figure out how many replicates are actually available for each time point
    left_na <- cbind(matrix(0,nrow(all_reps),1),mtx_side_count[,-ncol(mtx_side_count)]) # left shifted matrix
    right_na <- cbind(mtx_side_count[,-1],matrix(0,nrow(genes),1)) # right shifted matrix
    repmtx_l[[i]] <- repmtx - left_na - mtx_count[[i]] - right_na
    # to avoid division by 0 and induce NAs if there are no time points available
    repmtx_l[[i]][repmtx_l[[i]]==0] <- NA
  }

  dat <- genes # assigning to dataframe to return
  for (x in 0:(num_reps-1)){
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  }
  dat[is.na(genes)] <- NA # do not impute missing values
  return(dat)
}

#' Function to normalize expressions in a matricized manner, by row.
#'
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @return normalized expression data
#' @keywords internal
#' @noRd
normalize_all <- function(genes){

  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])

  # vector of row means
  all_row_mean <- rowMeans(all_reps, na.rm = TRUE)
  # vector of standard deviations
  all_row_stdev <- sqrt(rowSums((all_reps - rowMeans(all_reps,na.rm = TRUE))^2, na.rm = TRUE)/(dim(all_reps)[2] - 1))

  # get the matrix of normalized expressions
  all_reps_normal <- (all_reps - all_row_mean)/all_row_stdev
  # if standard deviation is 0, imposes NA, so this expression shouldn't be considered anyway
  # and is now constant
  all_reps_normal[is.na(all_reps_normal)] <- 0

  # create dataframe with normalized expressions
  dat <- genes
  dat[,-1] <- all_reps_normal
  dat[is.na(genes)] <- NA # do not impute missing values

  return(list("dat"=dat, "means"=all_row_mean, "stdevs"=all_row_stdev))
}

#' Function to remove linear trend for each expression over all data.
#'
#' @param timen vector of time points for all expressions
#' @param num_reps number of replicates
#' @param tied whether replicates are paired (tied) or unpaired (untied)
#' @param genes data frame of genes with the following specifications: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.
#' @return a data frame containing the detrended expressions
#' @keywords internal
#' @noRd
de_linear_trend_all <- function(timen,num_reps,tied,genes){
  all_rep <- as.matrix(genes[,-1]) # y values for linear fit
  if (!tied){ # if they're not paired, we just fit an aggregate data model

    # x values for linear fit
    xrow <- rep(timen,each=num_reps)
    xmtx <- matrix(rep(xrow,each=nrow(all_rep)),nrow = nrow(all_rep))

    # covariance
    cov <- rowSums((all_rep-rowMeans(all_rep,na.rm = TRUE))*(xmtx-rowMeans(xmtx)),na.rm = TRUE)
    # variance
    var <- rowSums((xmtx - rowMeans(xmtx))^2,na.rm = TRUE)

    # fitted coefficients - a+bx
    beta <- cov/var
    alph <- rowMeans(all_rep,na.rm = TRUE)-(beta*rowMeans(xmtx))

    # remove linear trend
    df <- all_rep-(alph+(beta*xmtx)) # linear fit
  } else { # we have to do the models separately for each replicate
    # preallocate matrix where we put results
    df <- matrix(NA, nrow = dim(all_rep)[1], ncol = dim(all_rep)[2])

    # x values for linear fit
    xmtx <- matrix(rep(timen,each=nrow(all_rep)),nrow = nrow(all_rep))
    for (i in 1:num_reps){
      each_rep <- all_rep[,seq(i,ncol(all_rep),by=num_reps)]

      # covariance
      cov <- rowSums((each_rep-rowMeans(each_rep,na.rm = TRUE))*(xmtx-rowMeans(xmtx)),na.rm = TRUE)
      # variance
      var <- rowSums((xmtx - rowMeans(xmtx))^2,na.rm = TRUE)
      # fitted coefficients - a+bx
      beta <- cov/var
      alph <- rowMeans(each_rep,na.rm = TRUE)-(beta*rowMeans(xmtx))

      # remove linear trend
      df[,seq(i,ncol(all_rep),by=num_reps)] <- each_rep -(alph+(beta*xmtx)) # linear fit
    }
  }
  # get the data frame correctly set up for returning
  res_df <- genes
  res_df[,-1] <- df

  return (res_df)
}
