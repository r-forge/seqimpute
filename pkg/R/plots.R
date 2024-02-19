#' Extract only the completed datasets from the results obtained with seqimpute function. 
#' Therefore, the original dataset with missing values is discarded,
#' together with the two first columns (".id" and ".imp")".
#' 
#' 
#' @param impdata a dataframe obtained with the seqimpute() function
#' @return Returns a data frame that has the form for the application of
#' Halpin's clustering strategy
#' 
#' @export
onlyimputed <- function(impdata){
  if(".imp"%in%colnames(impdata)&".id"%in%colnames(impdata)){
    to_remove <- c(".imp", ".id")
    impdata <- impdata[ impdata$.imp!=0, !(names(impdata) %in% to_remove)]
    return(impdata)
  }
  
}

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
#' @return Returns a state sequence object containing only
#' the sequences with at least one missing value
#' 
#' @export
seqwithmiss <- function(seqdata){
  if(!inherits(seqdata,"stslist")){
    stop("seqdata is not a sequence object.")
  }
  if(is.na(attr(seqdata,"nr"))){
    tmp <- seqdata
    for(i in 1:ncol(seqdata)){
      tmp[,i] <- as.character(seqdata[,i])
    }
    rowsNA <- rowSums(is.na(tmp))
  }else{
    rowsNA <- rowSums(seqdata==attr(seqdata,"nr"))
  }
  return(seqdata[rowsNA!=0,])
}


#' Plot of the most common patterns of missing data.
#' 
#' @param seqdata a state sequence object built with the TraMineR package
#' @param with.complete a logical stating if complete trajectories should be included or not in the plot
#'
#' @param ... parameters to be passed to the seqfplot function 
#' 
#' @export
seqmissfplot <- function(seqdata, with.complete=TRUE,...){
  if(with.complete==TRUE){
    seqmiss <- seqdata
  }else{
    seqmiss <- seqwithmiss(seqdata)
  }
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
#' @param with.complete a logical stating if complete trajectories should be included or not in the plot
#' @param ... parameters to be passed to the seqIplot function 
#' 
#' @export
seqmissIplot <- function(seqdata, with.complete=TRUE,...){
  if(with.complete==TRUE){
    seqmiss <- seqdata
  }else{
    seqmiss <- seqwithmiss(seqdata)
  }
  if(nrow(seqmiss)==0){
    stop("The provided sequence object has no missing values.")
    
  }
  misspatterns <- matrix(NA,nrow(seqmiss),ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  if(!is.na(attr(seqmiss,"nr"))){
    misspatterns[seqmiss==attr(seqmiss,"nr")] <- "missing"
    misspatterns[seqmiss!=attr(seqmiss,"nr")] <- "observed"
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns,alphabet = c("observed","missing"),cpal=c("blue","red"),xtstep=attr(seqmiss,"xtstep")))
    
    TraMineR::seqIplot(seqtest,...)
  }else{
    misspatterns[is.na(seqmiss)] <- "missing"
    misspatterns[!is.na(seqmiss)] <- "observed"
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns,alphabet = c("observed","missing"),cpal=c("blue","red"),xtstep=attr(seqmiss,"xtstep")))
    
    TraMineR::seqIplot(seqtest,...)
  }
  
}

#' Function built on the seqimplic function of the TraMineRextras package.
#' Visualization and identification of the states 
#' that best characterize sequence with missing data vs the sequences without missing data at each position (time point). 
#' See the seqimplic helps to more details on how it works.
#' 
#' @param seqdata a state sequence object built with the TraMineR package
#' @param ... parameters to be passed to the seqimplic() function 
#' 
#' @export
seqmissimplic <- function(seqdata,...){
  tt <- rep("missing",nrow(seqdata))
  tt[rowSums(seqdata==attr(seqdata,"nr"))==0] <- "observed"
  imp <- TraMineRextras::seqimplic(seqdata,tt)
  return(imp) 
}


#' Function that adds the clustering result to an imputed dataset
#' obtained with seqimpute
#' 
#' @param imputed a dataframe obtained with the seqimpute function
#' @param clustering clustering made on the stacked on multiple imputed datasets
#' 
#' @export
addcluster <- function(imputed, clustering){
  imputed$cluster <- NA
  if(!".imp"%in%colnames(imputed)){
    stop("The provided dataset does not have any columns .imp")
  }
  if(!0%in%imputed$.imp){
    imputed$cluster <- clustering
  }else{
    tt <- which(rowSums(is.na(imputed[imputed$.imp==0,]))==0)
    imputed[tt,] <- clustering[tt]
    imputed[imputed$.imp!=0,"cluster"] <- clustering
  }
  return(imputed)
}