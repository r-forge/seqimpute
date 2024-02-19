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