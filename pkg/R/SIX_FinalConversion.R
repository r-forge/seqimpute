# X. Final conversions ---------------------------------------------------------

################################################################################
# X - Final conversions
################################################################################
FinalResultConvert <- function(RESULT, ODClass, ODlevels, rownamesDataset, nrowsDataset, nr, nc, rowsNA, include, mi, mice.return){

 
  # ## Adjustment of the rendered form of RESULT
  if (include == FALSE) {    # case include == 1 (not including initial data
    # set OD)
    if (mi == 1) {       # case mi == 1, we render only the single matrix
      # ODi (without the numerotation variables aside)
      # Getting rid of the part of RESULT containing the matrix OD and the
      # column of variables '1' on the left-hand side
      RESULT <- RESULT[(nr+1):(2*nr),2:(nc+1)]
    } else {
      # Otherwise (i.e. if mi > 1), we just get rid of the part of RESULT
      # containing the matrix OD
      # (and we keep the variables numbering the imputations aside)
      RESULT <- RESULT[(nr+1):((mi+1)*nr),]
    }
    
  }
  # Else (meaning that we are in the case include == 2 (including initial
  # dataset OD)), we simply don't do any change in the form of RESULT (which
  # has already been constructed to fit the shape option '2')
  
  # Transformation of the columns of RESULT into numeric
  # So that it originally could fit the file "Simulation_Ascona_3"
  # of Andre
  RESULT <- apply(RESULT,2,as.numeric)
  

  # Transforming RESULT in a data frame
  RESULT <- as.data.frame(RESULT)
  
  # In case of a factor dataset OD:
  # RE-RECODING RESULT to go from "1", "2", etc. to "words"
  #*************************************
  if (ODClass == "factor") {
    if (mi == 1 & include == FALSE) {
      RESULT <- as.data.frame( sapply(RESULT, mapvalues,
                                      from = as.character(as.vector(1:length(ODlevels))),
                                      to = ODlevels, warn_missing=FALSE) )
    } else {
      # Taking account ot the special notation of RESULT that has an extra
      # column on the left of RESULT (as soon as mi > 1 or in any case if
      # include == 2)
      RESULT[,2:ncol(RESULT)] <- as.data.frame( sapply(RESULT[,2:ncol(RESULT)],
                                                       mapvalues,
                                                       from = as.character(as.vector(1:length(ODlevels))),
                                                       to = ODlevels,warn_missing=FALSE) )
    }
  }
  

  
  #### We put again the rows having only NA's discarder at the beginning
  if(length(rowsNA)>0){
    if (include == FALSE) {
      if (mi == 1) {
        for(i in 1:length(rowsNA)){
          if(rowsNA[i]==1){
            RESULT <- rbind(rep(NA,ncol(RESULT)),RESULT)
          }else if(rowsNA[i]==nrowsDataset){
            RESULT <- rbind(RESULT,rep(NA,ncol(RESULT)))
          }else{
            RESULT <- rbind(RESULT[1:(rowsNA[i]-1),],rep(NA,ncol(RESULT)),RESULT[rowsNA[i]:nrow(RESULT),])
          }
        }
      }else{
        for(j in 1:mi){
          for(i in 1:length(rowsNA)){
            if(j==1 & rowsNA[i]==1){
              RESULT <- rbind(c(j,rep(NA,ncol(RESULT))),RESULT)
            }else if(j==mi & rowsNA[length(rowsNA)]==nrowsDataset){
              RESULT <- rbind(RESULT,c(mi,rep(NA,ncol(RESULT))))
            }else{
              RESULT <- rbind(RESULT[1:(nrowsDataset*(j-1)+rowsNA[i]-1),],c(j,rep(NA,ncol(RESULT)-1)),RESULT[(nrowsDataset*(j-1)+rowsNA[i]):nrow(RESULT),])
            }
          }
        }
      }
    }else{
      for(j in 1:(mi+1)){
        for(i in 1:length(rowsNA)){
          if(j==1 & rowsNA[i]==1){
            RESULT <- rbind(c(j-1,rep(NA,ncol(RESULT))),RESULT)
          }else if(j==(mi+1) & rowsNA[length(rowsNA)]==nrowsDataset){
            RESULT <- rbind(RESULT,c(mi-1,rep(NA,ncol(RESULT))))
          }else{
            RESULT <- rbind(RESULT[1:(nrowsDataset*(j-1)+rowsNA[i]-1),],c(j-1,rep(NA,ncol(RESULT))),RESULT[(nrowsDataset*(j-1)+rowsNA[i]):nrow(RESULT),])
          }
        }
      }
    }
  }
  
  if(include==TRUE){
    RESULT <- cbind(RESULT[,1],rep(rownamesDataset,mi+1),RESULT[,2:ncol(RESULT)])
    colnames(RESULT)[1] <- ".imp"
    colnames(RESULT)[2] <- ".id"
  }else if(include==FALSE & mi>1){
    RESULT <- cbind(RESULT[,1],rep(rownamesDataset,mi),RESULT[,2:ncol(RESULT)])
    colnames(RESULT)[1] <- ".imp"
    colnames(RESULT)[2] <- ".id"
  }else{
    rownames(RESULT)<-rownamesDataset
  }
  
  if(mice.return==TRUE){
    RESULT <- as.mids(RESULT)
  }
  
  
  # Returning the matrix composed of every imputations
  return(RESULT)
}