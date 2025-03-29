
compute_auc <- function(data) {
  # Ensure the dataframe has the required columns
  if (!all(c("time", "conc") %in% names(data))) {
    stop("Data must contain 'time' and 'conc' columns.")
  }

  # Sort data by time to ensure correct integration
  data <- data[order(data$time), ]

  # Compute AUC using the trapezoidal rule
  auc <- sum(diff(data$time) * (head(data$conc, -1) + tail(data$conc, -1)) / 2)

  return(auc)
}

compute_cmax <- function(data) {
  # Ensure the dataframe has the required columns
  if (!all(c("time", "conc") %in% names(data))) {
    stop("Data must contain 'time' and 'conc' columns.")
  }

  # Compute Cmax
  cmax <- max(data$conc)

  return(cmax)
}

getAucCmaxDf <- function(df){
  # Get unique values of id, drug, and per
  unique_drugs <- unique(df$drug)
  unique_ids <- unique(df$id)
  unique_pers <- unique(df$per)

  # Initialize result storage
  results <- data.frame(id = numeric(), drug = character(), per = numeric(), auc = numeric() , cmax = numeric())

  # Loop over id, drug, and per
  for (drug in unique_drugs) {
    for (id in unique_ids) {
      for (per in unique_pers) {
        # Subset the dataframe
        subset_data <- df[df$id == id & df$drug == drug & df$per == per, ]

        # Check if subset is non-empty
        if (nrow(subset_data) > 1) {
          auc_value <- compute_auc(subset_data)
          cmax_value <- compute_cmax(subset_data)
          results <- rbind.data.frame(results, data.frame(id = id, drug = drug, per = per, auc = auc_value, cmax = cmax_value))
        }
      }
    }
  }
  return(results)
}


get_lmer <- function(nSubjects,pkParameterData,n_trials){
  #set up progress bar
  res <- lapply(1:n_trials,this_samp,pkParameterData,nSubjects)
  return(res)
}

this_samp <- function(trialNumber,pkParameterData,nSubjects){
  #create subset with n random individuals by drug for parallel design

  pkDf <- pkParameterData
  par_df <- NULL
  drugNames <- unique(pkDf$drug)
  for (drugName in drugNames){
    drugPkDf <- pkDf[pkDf$drug == drugName & pkDf$per == "1" ,]
    idsThisDrug <- sample(x = unique(drugPkDf$id),size = nSubjects,replace = FALSE)
    par_df <- rbind.data.frame(par_df,drugPkDf[drugPkDf$id %in% idsThisDrug,])
    pkDf <- pkDf[!(pkDf$id %in% idsThisDrug),]
  }

  #create subset with the same nSubjects random individuals by drug, period
  #for replicate (all periods) and cross-over (period 1 only)
  n_curves <- as.data.frame(table(pkParameterData$per, pkParameterData$drug))
  colnames(n_curves) <- c("per", "drug", "N")
  n_curves <- n_curves$N[1]
  to_keep <- sample(n_curves,nSubjects) #replace=F by default
  rep_df <- data.frame()  # Initialize empty data frame
  group_levels <- unique(pkParameterData[c("drug", "per")])  # Get unique group combinations
  for (i in seq_len(nrow(group_levels))) {
    # Extract current group
    subset_data <- pkParameterData[pkParameterData$drug == group_levels$drug[i] & pkParameterData$per == group_levels$per[i], ]
    # Select rows based on `to_keep`
    selected_rows <- subset_data[to_keep, , drop = FALSE]
    # Append to result
    rep_df <- rbind(rep_df, selected_rows)
  }

  cross_df <- rep_df %>% filter(per=="1" )

  lmPar <- list(rep(NA,2)) #for parallel
  lmCross <- list(rep(NA,2)) #for crossover
  lmerRep <- list(rep(NA,2)) #for replicate

  lmPar[[1]] <- suppressWarnings(nlme::lme(fixed = log10(auc) ~  drug, random = ~ 1 | per, data = par_df))
  lmPar[[2]] <- suppressWarnings(nlme::lme(fixed = log10(cmax) ~  drug, random = ~ 1 | per, data = par_df))
  lmCross[[1]] <- suppressWarnings(nlme::lme(fixed = log10(auc) ~  drug, random = ~ 1 | id, data = cross_df))
  lmCross[[2]] <- suppressWarnings(nlme::lme(fixed = log10(cmax) ~  drug, random = ~ 1 | id, data = cross_df))

  ctrl <- lme4::lmerControl(optimizer ="Nelder_Mead")
  periods <- unique(pkParameterData$per)
  n_pers <- length(periods) #each period is a replicate for that drug
  if(n_pers>1) {
    lmerRep[[1]] <- suppressWarnings(lmerTest::lmer(log10(auc) ~  per + drug + (1 | id), data = rep_df, control=ctrl))
    lmerRep[[2]] <- suppressWarnings(lmerTest::lmer(log10(cmax) ~  per + drug + (1 | id), data = rep_df, control=ctrl))
  }
  return(list(par = lmPar, cross = lmCross, rep = lmerRep))
}

#' Generates bioequivalence statistics
#'
#' Generates BE statistics for parallel, cross-over and full replication designs.
#' From simulated data with 1000 subjects per drug, the function first calculates
#' area under the concentration-time curve (AUC) and maximum concentration (Cmax)
#' for each simulated profiles. Periods are defined as >1 for each replicated
#' administration (e.g., when intra-occasion varibility is simulated.)
#' Linear regression (\code{\link{lm}}) is used on \code{log10(auc)} and \code{log10(cmax)}
#' for parallel and cross-over designs, and mixed effect modeling (\code{\link{lmer}})
#' is used for the replicated design.
#'
#' @title Calculation of BE statistics
#' @param data A dataframe or tibble  with these columns:
#' subject id, time, concentrations, period, drug.
#' @param n_trials Integer coding the number of random trials per sample size.
#' @param subj_min Integer defining the minimum number of subjects per trial
#' @param subj_max Integer defining the maximum number of subjects per trial
#' @param subj_step Integer defining the step size between \code{subj_min} and
#' \code{subj_max}, so that the program will test sample sizes
#' \emph{N = \{subj_min, subj_min + subj_step, subj_min + 2*subj_step,...,subj_max\}}.
#' @param seed Random number seed.
#' @return A list of lists: \[\[1:n_samp\]\]\[\[1:\code{n_trials}\]\]\[\[$par, $cross,
#' $rep\]\]\[\[1,2]\], where n_samp is the number of sample sizes tested, e.g.,
#' (\code{subj_max} - \code{subj-min}) / \code{subj_step}.  The first item in the
#' design is auc, and the second is cmax.  For example the comparison of auc using
#' the third trial of the second sample size for the replicated design would be
#' x\[\[2\]\]\[\[3\]\]$rep\[\[1\]\].
#' @seealso \code{\link{extractBE}}, \code{\link{plotBE}}
#' @author Michael Neely
#' @author Abdullah Hamadeh
#' @export

calc_be <- function(data,
                    n_trials = 50,
                    subj_min = 5, subj_max = 50, subj_step = 5,
                    seed = -17, ci_lcut=.8, ci_ucut=1.25){

  pkParameterData <- getAucCmaxDf(data)

  #extract data
  periods <- unique(pkParameterData$per)
  n_pers <- length(periods) #each period is a replicate for that drug
  nSubjects <- seq(from = subj_min, to = subj_max, by = subj_step) #step through sample sizes
  n_samp <- length(nSubjects)
  set.seed(seed)



  res_list <- lapply(nSubjects, get_lmer, pkParameterData, n_trials) #purrr::map(nSubjects,get_lmer,n_trials,pkParameterData)
  return(res_list)

}
