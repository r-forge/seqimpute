#' Generation of missing data under the form of gaps, which
#' is the typical form of missing data with longitudinal data.
#' A missing completely at random (MCAR) mechanism is used.
#' 
#' @param data a data frame containing sequences of a multinomial variable without missing
#' data, on which missing data are simulated
#' @param pstart probability to start a missing data
#' @param propdata proportion of the observations on which missing data will be simulated
#' 
#' @export
make_missing_MCAR<- function(data,pstart=0.1,propdata=0.6){
  sizehalf <- round(propdata*nrow(data))
  rowsmiss <- sample(1:nrow(data),size=sizehalf,replace=FALSE)
  matrix_missing <- matrix(NA,nrow(data),ncol(data))
  for(i in 1:length(rowsmiss)){
    nmis <- ncol(data)
    while(nmis>floor(0.75*ncol(data))){
      for(j in 1:ncol(data)){
        if(j==1){
          matrix_missing[rowsmiss[i],j] <- sample(x=c(0,1),size=1,p=c(pstart,1-pstart))
        }else{
          if(matrix_missing[rowsmiss[i],j-1]==1){
            matrix_missing[rowsmiss[i],j] <- sample(x=c(0,1),size=1,p=c(pstart,1-pstart))
          }else{
            matrix_missing[rowsmiss[i],j] <- sample(x=c(0,1),size=1,p=c(0.66,0.34))
          }
        }
      }
      nmis <- sum(matrix_missing[rowsmiss[i],]==0)
    }
  }
  data[matrix_missing==0] <- NA
  return(data)
}


#' Generation of missing data under the form of gaps, which
#' is the typical form of missing data with longitudinal data.
#' A missing completely at random (MAR) mechanism is used.
#' 
#' @param data a data frame containing sequences of a multinomial variable without missing
#' data, on which missing data are simulated
#' @param states_high list of states that will have a larger probability to trigger a 
#' @param pstart_high probability to start a missing data for the specified states 
#' @param pstart_low probability to start a missing data for the other states
#' @param propdata proportion of the observations on which missing data will be simulated
#' subsequent gap of missing data
#' 
#' @export
make_missing_MAR<- function(data,pstart_high=0.2,pstart_low=0.03,propdata=0.6,states_high){
  sizehalf <- round(propdata*nrow(data))
  rowsmiss <- sample(1:nrow(data),size=sizehalf,replace=FALSE)
  matrix_missing <- matrix(NA,nrow(data),ncol(data))
  for(i in 1:length(rowsmiss)){
    nmis <- ncol(data)
    while(nmis>floor(0.75*ncol(data))){
      for(j in 1:ncol(data)){
        if(j==1){
          matrix_missing[rowsmiss[i],j] <- sample(x=c(0,1),size=1,p=c(pstart_low,1-pstart_low))
        }else{
          if(matrix_missing[rowsmiss[i],j-1]==1){
            if(data[rowsmiss[i],j-1]%in%states_high){
              matrix_missing[rowsmiss[i],j] <- sample(x=c(0,1),size=1,p=c(pstart_high,1-pstart_high))
            }else{
              matrix_missing[rowsmiss[i],j] <- sample(x=c(0,1),size=1,p=c(pstart_low,1-pstart_low))
              
            }
          }else{
            matrix_missing[rowsmiss[i],j] <- sample(x=c(0,1),size=1,p=c(66,34))
          }
        }
      }
      nmis <- sum(matrix_missing[rowsmiss[i],]==0)
    }
  }
  data[matrix_missing==0] <- NA
  return(data)
}