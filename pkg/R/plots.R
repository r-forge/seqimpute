
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
  if(nrow(seqmiss)==0){
    stop("The provided sequence object has no missing values. Note that the 
            default behavior of the seqdef() function is to delete missing values 
            at the end of the sequences. If it is not what you want, please consider
            modifying it.")
    
  }
  misspatterns <- matrix(NA,nrow(seqmiss),ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  if(!is.na(attr(seqmiss,"nr"))){
    misspatterns[seqmiss==attr(seqmiss,"nr")] <- "missing"
    misspatterns[seqmiss!=attr(seqmiss,"nr")] <- "not missing"
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns,alphabet = c("not missing","missing"),cpal=c("blue","red"),xtstep=attr(seqmiss,"xtstep")))
    
    TraMineR::seqIplot(seqtest,...)
  }else{
    misspatterns[is.na(seqmiss)] <- "missing"
    misspatterns[!is.na(seqmiss)] <- "not missing"
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns,alphabet = c("not missing","missing"),cpal=c("blue","red"),xtstep=attr(seqmiss,"xtstep")))
    
    TraMineR::seqIplot(seqtest,...)
  }
  
}

#' Function built on the seqimplic function of the TraMineRextras package.
#' Visualization and identification of the states 
#' that best characterize sequence with missing data vs the sequences without missing data at each position (time point). 
#' See the seqimplic helps to more details on how it works.
#' 
#' @param seqdata a state sequence object built with the TraMineR package
#' @param ... parameters to be passed to the seqIplot function 
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