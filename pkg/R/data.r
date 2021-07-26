#' Dataset containing variables about the game addiction of young subjects
#'
#' An original dataset example \code{OD} to test \code{seqimpute.R} and its aside functions.
#'
#' \itemize{
#'   \item 500 sequences (i.e. 500 rows)
#'   \item 4 time measurements (i.e. 4 columns)
#'   \item The variables can be either 'no', 'yes' or NA
#' }
#'
#' @docType data
#' @keywords datasets
#' @name OD
#' @usage data(OD)
#' @format A data frame of factor variables with 500 sequences and 4 columns.
"OD"





#' Dataset containing 3 fixed covariates about the game addiction of young subjects
#'
#' These covariates are respectively '\code{Gender}' (\code{male}/\code{female}), '\code{Age}' (at \code{T1}, in years) and '\code{Track}' (\code{school}/\code{apprenticeship}).
#'
#' \itemize{
#'   \item 500 sequences
#'   \item 3 columns
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CO
#' @usage data(CO)
#' @format A data frame of fixed covariates with 500 sequences and 3 columns.
"CO"



#' Dataset containing 1 time-dependant covariate about the game addiction of young subjects
#'
#' This time-dependant covariate is the '\code{Gambling}' (\code{no}/\code{gambler}/\code{problematic gambler}) and contains thus the same number of columns as the original data frame \code{OD}: 4 columns.
#'
#' \itemize{
#'   \item 500 sequences
#'   \item 4 columns
#' }
#'
#' @docType data
#' @keywords datasets
#' @name COt
#' @usage data(COt)
#' @format A data frame of time-dependant covariates with 500 sequences and 4 columns.
"COt"




#' Dataset specifically designed to test specially located gap imputation
#'
#'
#' \itemize{
#'   \item 1001 sequences
#'   \item 25 columns
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SLGData
#' @usage data(SLGData)
#' @formatA data frame of factor variables with 500 sequences and 25 columns.
"SLGData"


NULL
