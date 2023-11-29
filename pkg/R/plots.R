
#' Extract all the sequences without missing value.
#' 
#' @param seqdata a state sequence object built with the TraMineR package
#' @return Returns either a data frame or a state sequence object, depending
#' the type of data that was provided to the function
#' 
#' @export
seqcomplete <- function(seqdata){
  if(!inherits(seqdata,"stslist")){
    stop("seqdata is not a sequence object.")
  }
  rowsNA <- rowSums(seqdata==attr(seqdata,"nr"))
  return(seqdata[rowsNA==0,])
}

#' Extract all the sequences with at least one missing value
#' 
#' @param seqdata a state sequence object built with the TraMineR package
#' @return Returns either a data frame or a state sequence object, depending
#' the type of data that was provided to the function
#' 
#' @export
seqwithmiss <- function(seqdata){
  if(!inherits(seqdata,"stslist")){
    stop("seqdata is not a sequence object.")
  }
  rowsNA <- rowSums(seqdata==attr(seqdata,"nr"))
  return(seqdata[rowsNA!=0,])
}


#' Plot of the most common patterns of missing data.
#' 
#' @param seqdata a state sequence object built with the TraMineR package
#' @param ... parameters to be passed to the seqfplot function 
#' 
#' @export
seqmissfplot <- function(seqdata,...){
  seqmiss <- seqwithmiss(seqdata)
  misspatterns <- matrix(NA,nrow(seqmiss),ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  
  misspatterns[seqmiss==attr(seqmiss,"nr")] <- "missing"
  misspatterns[seqmiss!=attr(seqmiss,"nr")] <- "observed"
  seqtest <- suppressMessages(TraMineR::seqdef(misspatterns,alphabet = c("observed","missing"),cpal=c("blue","red"),xtstep=attr(seqmiss,"xtstep")))
  
  TraMineR::seqfplot(seqtest,...)
  
}


#' Plot all the patterns of missing data.
#' 
#' @param seqdata a state sequence object built with the TraMineR package
#' @param ... parameters to be passed to the seqIplot function 
#' 
#' @export
seqmissIplot <- function(seqdata,...){
  seqmiss <- seqwithmiss(seqdata)
  misspatterns <- matrix(NA,nrow(seqmiss),ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  
  misspatterns[seqmiss==attr(seqmiss,"nr")] <- "missing"
  misspatterns[seqmiss!=attr(seqmiss,"nr")] <- "observed"
  seqtest <- suppressMessages(TraMineR::seqdef(misspatterns,alphabet = c("observed","missing"),cpal=c("blue","red"),xtstep=attr(seqmiss,"xtstep")))
  
  TraMineR::seqIplot(seqtest,...)
  
}