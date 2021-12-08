# C. File containing function used several time for different part of seqimpute
#
#
################################################################################
# Consider time-dependent covariates

COtselection <- function(COtselected, COt, ncot, t1, t2, nr, nc, j, shifted){
  COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
  for (d in 1:(ncot/nc)) {
    COttemp <- cbind(COttemp, COt[,(j+shifted) + (d-1)*nc])
  }
  COtselected[t1:t2,] <- COttemp
  
  return(COtselected)
}


################################################################################
# Concatenating CDi with COtselected_i (the matrix containing the current
# time-dependent covariates) Checking if COt is NOT completely empty

COtselectionSpe <- function(CDi, COt, ncot, nc, i, j, k){
  if (ncot > 0) {
    COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
    for (d in 1:(ncot/nc)) {
      COtselected_i <- cbind(COtselected_i, COt[i,(j) + (d-1)*nc])
    }
    COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
    # Concatenating CDi and COtselected_i
    # into CDi
    CDi <- cbind(CDi, COtselected_i)
    # Transformation of the names of the columns
    # of CDi (called V1, V2, ..., "VtotV")
    colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
  }
  return(CDi)
}


################################################################################
# Creation of matrices REFORD with ORDER

REFORDInit <- function(ORDER, nr, nc){
  # The purpose of this function is to accelerate part 3.3 in which we
  # initially (i.e. with the first versions of seqimpute3.R) go "order"
  # times through the matrix ORDER
  #
  # Going one single time through the matrix ORDER, we create MaxGap
  # matrices REFORD (numbered from REFORD_1 to "REFORD_MaxGap") which
  # collect the coordinates of each corresponding values greater than 0
  # It will create MaxGap lookup matrices that will be used in point 3.3
  # to directly pinpoint the NA that have to be currently imputed
  # according to the value of the variable "order"
  # Updating MaxGap
  
  MaxGap <- max(ORDER[ORDER!=0])-(min(ORDER[ORDER!=0]) - 1)
  
  # Initialization of the REFORD matrices
  REFORD_L <- lapply(1:MaxGap, matrix, data = NA, nrow=0, ncol=2)
  
  # Return the matrix of coordinates of value of ORDER bigger than 0
  non_zero <- which(ORDER > 0, arr.ind = TRUE)
  non_zero <- non_zero[(non_zero[,1] <= nr) & (non_zero[,2] <= nc),,drop=F]
  # Updating ORDER so that it becomes a matrix with positive
  # values going from 1 to MaxGap
  ORDER[non_zero] <- ORDER[non_zero] - (min(ORDER[ORDER!=0]) - 1)
  
  # Collecting the coordinates for each value of "order"
  non_zero <- non_zero[order(non_zero[,1]),,drop=F]
  ord_cord <- ORDER[non_zero]
  
  for (i in 1:MaxGap) {
    REFORD_L[[i]] <- non_zero[which(ord_cord == i),,drop=F]
  }
  
  return(list(MaxGap, REFORD_L, ORDER))
}

################################################################################
# Creation of matrices REFORD with GapSize

REFORDInit_TI <- function(GapSize, nr, ORDER, GapSizelist){
  
  # Initialization of the REFORD matrices
  REFORD_L <- lapply(1:GapSize, matrix, data = NA, nrow = 0, ncol=2)
  
  # Return the matrix of coordinates of value of ORDER bigger than 0
  non_zero <- which(ORDER > 0, arr.ind = TRUE)
  non_zero <- non_zero[(non_zero[,1] <= nr) & (non_zero[,2]  %in% GapSizelist),,drop=F]
  
  # Collecting the coordinates for each value of "order"
  ord_cord <- ORDER[non_zero,drop=F]
  
  for (i in 1:GapSize) {
    REFORD_L[[i]] <- non_zero[which(ord_cord == i),]
  }
  
  return(REFORD_L)
}



################################################################################
# Past VI computation

PastVICompute <- function(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, COtsample, COt, pastDistrib, futureDistrib, k,timing){
  
  CDp <- matrix(NA, nrow=nr*ud, ncol=np)

  if (ncot > 0) {
    # initialisation of matrix COtselected
    COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
  }
  
  # initialisation of matrix CDdb (for past
  # distribution analysis) (Distribution Before)
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
    db <- matrix(NA, nrow=nr, ncol=k)
    # submatrix of CDdb:
    # CDdb is composed of
    # ud matrix db on top
    # of each other
  }
  
  # initialisation of matrix CDda (for future
  # distribution analysis) (Distribution After)
  if (futureDistrib) {
    CDda <- matrix(NA, nrow=nr*ud, ncol=k)
    # CDda has same dimensions as CDdb
    da <- matrix(NA, nrow=nr, ncol=k)
    # da has same dimensions as db
  }
  for (j in frameSize:nc){
    # /!\ j is initialised at
    # the very end (utmost right) of the
    # frame
    t1 <- (nr*(iter-1)+1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr*iter
    # VD
    CD[t1:t2,1] <- OD[,j-frameSize+np+1]
    if(timing==T){
      CD[t1:t2,2]<-j-frameSize+np+1
    }
    
    # past VIs
    CDp[t1:t2,] <- OD[,(j-frameSize+1):(j-frameSize+np)] 
    
    
    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize
    
    # Eventually considering time-dependent
    # covariates
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, nc, j, shifted = -frameSize+np+1)
    }
    # Past distribution (i.e. Before)
    if (pastDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[(1:(j-frameSize+np)),], factor, levels=c(1:k,NA), exclude=NULL)
      # because:
      # j-frameSize+np+1 - 1 = j-frameSize+np
      
      db_list <- lapply(tempOD,summary)
      db_matrix <- do.call(rbind,db_list)
      CDdb[t1:t2,] <- db_matrix[,1:k]/length(1:(j-frameSize+np))
    }
    # Future distribution (i.e. After)
    if (futureDistrib) {
      if ( (j-frameSize+np+2) <= nc) {  
        ODt <- t(OD)
        ODt <- as.data.frame(ODt)
        tempOD <- lapply(ODt[((j-frameSize+np+2):nc),], factor, levels=c(1:k,NA), exclude=NULL)
        # because:
        # j-frameSize+np+1 + 1
        # = j-frameSize+np+2
        
        da_list <- lapply(tempOD,summary)
        da_matrix <- do.call(rbind,da_list)
        CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np+2):nc)
      } else {
        # if index in OD exceeds OD number of
        # columns, the future distribution of
        # the k categorical variables is simply
        # null for everyone of them
        CDda[t1:t2,] <- matrix(nrow=nr,ncol=k,0)
      }
    }
    
    iter <- iter+1
  }
  
  
  # past VIs
  CD <- cbind(CD, CDp) 
  
  
  if (pastDistrib) {
    CD <- cbind(CD, CDdb)
  }
  
  if (futureDistrib) {
    CD <- cbind(CD, CDda)
  }
  
  # Conversion of CD into a data frame
  CD <- as.data.frame(CD)
  
  # Eventually concatenating CD with COs (the matrix
  # containing the covariates)
  if (all(is.na(CO))==FALSE) { # Checking if CO is NOT
    # completely empty
    # Creation of the stacked covariates
    # matrix for 3.1
    COs <- do.call("rbind", rep(list(CO), ud))
    # Concatenating CD and COs into CD
    CD <- cbind(CD, COs)
  }
  # Else, in case CO is empty (i.e. we don't consider
  # any covariate) simply continue with the current CD
  
  # Eventually concatenating CD with COtselected (the
  # matrix containing the current time-dependent
  # covariates)
  # Checking if COt is NOT completely empty
  if (ncot > 0) {
    # Concatenating CD and COtselected into CD
    CD <- cbind(CD, as.data.frame(COtselected))
  }
  return(CD)
}


################################################################################
# Future VI computation

FutureVICompute <- function(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, COtsample, COt, pastDistrib, futureDistrib, k, nf,timing){
  
  CDf <- matrix(NA, nrow=nr*ud, ncol=nf)

  if (ncot > 0) {
    # initialisation of matrix COtselected
    COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
  }
  
  # initialisation of matrix CDdb (for past
  # distribution analysis) (Distribution Before)
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
    db <- matrix(NA, nrow=nr, ncol=k)
    # submatrix of CDdb:
    # CDdb is composed of
    # ud matrix db on top
    # of each other
  }
  
  # initialisation of matrix CDda (for future
  # distribution analysis) (Distribution After)
  if (futureDistrib) {
    CDda <- matrix(NA, nrow=nr*ud, ncol=k)
    # CDda has same dimensions as CDdb
    da <- matrix(NA, nrow=nr, ncol=k)
    # da has same dimensions as db
  }
  for (j in frameSize:nc){
    # /!\ j is initialised at
    # the very end (utmost right) of the
    # frame
    t1 <- (nr*(iter-1)+1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr*iter
    # VD
    CD[t1:t2,1] <- OD[,j-frameSize+np+1]
    
    if(timing==T){
      CD[t1:t2,2]<-j-frameSize+np+1
    }
    
    # future VIs
    CDf[t1:t2,] <- OD[,(j-nf+1):j]
    
    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize
    
    # Eventually considering time-dependent
    # covariates
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, nc, j, shifted = -frameSize+np+1)
    }
    # Past distribution (i.e. Before)
    if (pastDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[(1:(j-frameSize+np)),], factor, levels=c(1:k,NA), exclude=NULL)
      # because:
      # j-frameSize+np+1 - 1 = j-frameSize+np
      
      db_list <- lapply(tempOD,summary)
      db_matrix <- do.call(rbind,db_list)
      CDdb[t1:t2,] <- db_matrix[,1:k]/length(1:(j-frameSize+np))
    }
    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j-frameSize+np+2):nc),], factor, levels=c(1:k,NA), exclude=NULL)
      # because:
      # j-frameSize+np+1 + 1
      # = j-frameSize+np+2
      
      da_list <- lapply(tempOD,summary)
      da_matrix <- do.call(rbind,da_list)
      CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np+2):nc)
    }
    iter <- iter+1
  }
  
  
  # future VIs
  CD <- cbind(CD, CDf) 
  
  if (pastDistrib) {
    CD <- cbind(CD, CDdb)
  }
  
  if (futureDistrib) {
    CD <- cbind(CD, CDda)
  }
  
  # Conversion of CD into a data frame
  CD <- as.data.frame(CD)
  
  # Eventually concatenating CD with COs (the matrix
  # containing the covariates)
  if (all(is.na(CO))==FALSE) { # Checking if CO is NOT
    # completely empty
    # Creation of the stacked covariates
    # matrix for 3.1
    COs <- do.call("rbind", rep(list(CO), ud))
    # Concatenating CD and COs into CD
    CD <- cbind(CD, COs)
  }
  # Else, in case CO is empty (i.e. we don't consider
  # any covariate) simply continue with the current CD
  
  # Eventually concatenating CD with COtselected (the
  # matrix containing the current time-dependent
  # covariates)
  # Checking if COt is NOT completely empty
  if (ncot > 0) {
    # Concatenating CD and COtselected into CD
    CD <- cbind(CD, as.data.frame(COtselected))
  }
  return(CD)
}


################################################################################
# Past or future VI computation

PastFutureVICompute <- function(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, COtsample, COt, pastDistrib, futureDistrib, k, nf, shift,timing){
                              
  CDp <- matrix(NA, nrow=nr*ud, ncol=np)
  CDf <- matrix(NA, nrow=nr*ud, ncol=nf)
  

  if (ncot > 0) {
    # initialisation of matrix COtselected
    COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
  }
  
  # initialisation of matrix CDdb (for past
  # distribution analysis) (Distribution Before)
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
    db <- matrix(NA, nrow=nr, ncol=k)
    # submatrix of CDdb:
    # CDdb is composed of
    # ud matrix db on top
    # of each other
  }
  # initialisation of matrix CDda (for future
  # distribution analysis) (Distribution After)
  if (futureDistrib) {
    CDda <- matrix(NA, nrow=nr*ud, ncol=k)
    # CDda has same dimensions as CDdb
    da <- matrix(NA, nrow=nr, ncol=k)
    # da has same dimensions as db
  }
  for (j in frameSize:nc){
    # /!\ j is initialised at
    # the very end (utmost right) of the
    # frame
    t1 <- (nr*(iter-1)+1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr*iter
    # VD
    CD[t1:t2,1] <- OD[,j-frameSize+np+1+shift]
    if(timing==T){
      CD[t1:t2,2] <- j-frameSize+np+1+shift
    }
    
    # past VIs
    CDp[t1:t2,] <- OD[,(j-frameSize+1):(j-frameSize+np)] 
    
    # future VIs
    CDf[t1:t2,] <- OD[,(j-nf+1):j]
    
    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize
    
    # Eventually considering time-dependent
    # covariates
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, nc, j, shifted = -frameSize+np+1+shift)
    }
    # Past distribution (i.e. Before)
    if (pastDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[(1:(j-frameSize+np+shift)),], factor, levels=c(1:k,NA), exclude=NULL)
      # because:
      # j-frameSize+np+1 - 1 = j-frameSize+np
      
      db_list <- lapply(tempOD,summary)
      db_matrix <- do.call(rbind,db_list)
      CDdb[t1:t2,] <- db_matrix[,1:k]/length(1:(j-frameSize+np+shift))
    }
    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j-frameSize+np+shift+2):nc),], factor, levels=c(1:k,NA), exclude=NULL)
      # because:
      # j-frameSize+np+1 + 1
      # = j-frameSize+np+2
      da_list <- lapply(tempOD,summary)
      da_matrix <- do.call(rbind,da_list)
      CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np+shift+2):nc)
       
    }
    
    iter <- iter+1
  }
  
  # past and future VIs
  CD <- cbind(CD, CDp, CDf)
  
  
  
  if (pastDistrib) {
    CD <- cbind(CD, CDdb)
  }
  
  if (futureDistrib) {
    CD <- cbind(CD, CDda)
  }
  
  # Conversion of CD into a data frame
  CD <- as.data.frame(CD)
  # Eventually concatenating CD with COs (the matrix
  # containing the covariates)
  if (all(is.na(CO))==FALSE) { # Checking if CO is NOT
    # completely empty
    # Creation of the stacked covariates
    # matrix for 3.1
    COs <- do.call("rbind", rep(list(CO), ud))
    # Concatenating CD and COs into CD
    CD <- cbind(CD, COs)
  }
  # Else, in case CO is empty (i.e. we don't consider
  # any covariate) simply continue with the current CD
  # Eventually concatenating CD with COtselected (the
  # matrix containing the current time-dependent
  # covariates)
  # Checking if COt is NOT completely empty
  if (ncot > 0) {
    # Concatenating CD and COtselected into CD
    CD <- cbind(CD, as.data.frame(COtselected))
  }
  return(CD)
}



################################################################################
# Imputation where only PAST VIs  exist

ODiImputePAST <- function(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise,timing){
  for (u in 1:nr_REFORD) {
    i <- REFORD[u,1]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- REFORD[u,2]
    # taking out the second coordinate
    # (column number in ORDER) from REFORD
    
    if(timing==F){
      CDi <- matrix(NA, nrow=k, ncol=1)
    }else{
      CDi <- matrix(NA,nrow=k,ncol=2)
      CDi[,2]<-j
    }
    
    # Matrix for past values
    vect <- LOOKUP[i,(j-np):(j-1)]
    # /!\ current pointer
    # on olumn is thus: "j"
    CDpi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)
    
    # Matrix for past distribution
    if (pastDistrib) {
      dbi <- summary(factor(LOOKUP[i,1:(j-1)], levels=c(1:k), exclude=NULL))/length(1:(j-1))
      CDdbi <- matrix(dbi[1:k], nrow=k, ncol=k, byrow=TRUE)
    }
    
    # Matrix for future distribution
    if (futureDistrib) {
      dai <- summary(factor(LOOKUP[i,(j+1):nc], levels=c(1:k), exclude=NULL))/length((j+1):nc)
      CDdai <- matrix(dai[1:k], nrow=k, ncol=k, byrow=TRUE)
    }
    
    # Concatenating CDi
    CDi <- cbind(CDi, CDpi)
    
    if (pastDistrib) {
      CDi <- cbind(CDi, CDdbi)
    }
    
    if (futureDistrib) {
      CDi <- cbind(CDi, CDdai)
    }
    
    # Conversion of CDi into a data frame
    CDi <- as.data.frame(CDi)
    # Type transformation of the columns of CDi
    # The first values of CDi must be of type factor
    # (categorical values)

    # We also account for the fact that levels that do not appear at
    # all in a given variable of CD were discarded with droplevels before
    # the fit of the mlogit model
    
    if(timing==FALSE){
      if(regr!="rf"){
        for(v in 1:(1+np+nf)){
          CDi[,v]<-factor(CDi[,v],levels=levels(CD[,v]),exclude=NULL)
        }
      }else{
        for(v in 2:(1+np+nf)){
          CDi[,v]<-factor(CDi[,v],levels=c(1:(k+1)))
          CDi[,v][is.na(CDi[,v])]<-k+1
        }
        CDi[,1]<-factor(CDi[,1],levels=c(1:k))
      }
      # The last values of CDi must be of type numeric
      # (distributions)
      if (pastDistrib | futureDistrib) {
        CDi[,(1+np+nf+1):ncol(CDi)] <- lapply(CDi[,(1+np+nf+1):ncol(CDi)],as.numeric)
      }
    }else{
      if(regr!="rf"){
        for(v in c(1,3:(2+np+nf))){
          CDi[,v]<-factor(CDi[,v],levels=levels(CD[,v]),exclude=NULL)
        }
      }else{
        for(v in c(3:(2+np+nf))){
          CDi[,v]<-factor(CDi[,v],levels=c(1:(k+1)))
          CDi[,v][is.na(CDi[,v])]<-k+1
        }
        CDi[,1]<-factor(CDi[,1],levels=c(1:k))
      }
      # The last values of CDi must be of type numeric
      # (distributions)
      if (pastDistrib | futureDistrib) {
        CDi[,(2+np+nf+1):ncol(CDi)] <- lapply(CDi[,(2+np+nf+1):ncol(CDi)],as.numeric)
      }
    }
    # Eventually concatenating CDi with COi
    # (the matrix containing the covariates)
    if (all(is.na(CO))==FALSE) {
      # Checking if CO is NOT
      # completely empty
      # Creation of the matrix COi used in 3.3
      COi <- do.call(rbind, replicate(k, as.data.frame(CO[i,]), simplify=FALSE))
      # Concatenating CDi and COi into CDi
      CDi <- cbind(CDi, COi)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }
    # Else, in case CO is empty (i.e. we don't
    # consider any covariate)
    # simply continue with the current CDi
    
    # Eventually concatenating CDi with
    # COtselected_i (the matrix containing the
    # current time-dependent covariates)
    # Checking if COt is NOT completely empty
    if (ncot > 0) {
      COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
      for (d in 1:(ncot/nc)) {
        COtselected_i <- cbind(COtselected_i, COt[i,(j) + (d-1)*nc])
      }
      COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
      # Concatenating CDi and COtselected_i
      # into CDi
      CDi <- cbind(CDi, COtselected_i)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }
    
    
    # Check for missing-values among predictors
    if (max(is.na(CDi[1,2:totV]))==0){
      ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
      
    }
  }
  return(ODi)
}



################################################################################
# Imputation where past and future VIs exist

ODiImputePF <- function(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise, shift, MaxGap, order,timing){
  for (u in 1:nr_REFORD) {
    i <- REFORD[u,1]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- REFORD[u,2]
    # taking out the second coordinate
    # (column number in ORDER) from REFORD
    if(timing==F){
      CDi <- matrix(NA,nrow=k, ncol=1)
    }else{
      CDi <- matrix(NA,nrow=k,ncol=2)
      CDi[,2]<-j
    }    
    # Matrix for past values
    shift <- as.numeric(shift)
    vect <- LOOKUP[i,(j-shift-np):(j-shift-1)]

    CDpi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)
    
    
    # Matrix for future values
    vect <- LOOKUP[i,(j-shift+MaxGap-order+1):(j-shift+MaxGap-order+nf)]
    CDfi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)
    
    # Matrix for past distribution
    if (pastDistrib) {
      dbi <- summary(factor(LOOKUP[i,1:(j-1)], levels=c(1:k), exclude=NULL))/length(1:(j-1))
      CDdbi <- matrix(dbi[1:k], nrow=k, ncol=k, byrow=TRUE)
    }
    # Matrix for future distribution
    if (futureDistrib) {
      dai <- summary(factor(LOOKUP[i,(j+1):nc], levels=c(1:k), exclude=NULL))/length((j+1):nc)
      CDdai <- matrix(dai[1:k], nrow=k, ncol=k, byrow=TRUE)
    }
    # Concatenating CDi
    CDi <- cbind(CDi, CDpi, CDfi)
    
    if (pastDistrib) {
      CDi <- cbind(CDi, CDdbi)
    }
    
    if (futureDistrib) {
      CDi <- cbind(CDi, CDdai)
    }
    
    # Conversion of CDi into a data frame
    CDi <- as.data.frame(CDi)
    # Type transformation of the columns of CDi
    # The first values of CDi must be of type factor
    # (categorical values)
    # We also account for the fact that levels that do not appear at
    # all in a given variable of CD were discarded with droplevels before
    # the fit of the mlogit model
    if(timing==FALSE){
      if(regr=="lm"|regr=="lrm"|regr=="mlogit"){
        for(v in 1:(1+np+nf)){
          CDi[,v]<-factor(CDi[,v],levels=levels(CD[,v]),exclude=NULL)
        }
      }else if(regr=="rf"){
        for(v in 2:(1+np+nf)){
          CDi[,v]<-factor(CDi[,v],levels=c(1:(k+1)))
          CDi[,v][is.na(CDi[,v])]<-k+1
        }
        CDi[,1]<-factor(CDi[,1],levels=c(1:k))
      }else{
        #glmnet
        CDi[,1] <- factor(CDi[,1],levels=c(1:k))
        for(v in 2:(1+np+nf)){
          CDi[,v]<-factor(CDi[,v],levels=levels(CD[,v]),exclude=NULL)
        }
      }
      # The last values of CDi must be of type numeric
      # (distributions)
      if (pastDistrib | futureDistrib) {
        CDi[,(1+np+nf+1):ncol(CDi)] <- lapply(CDi[,(1+np+nf+1):ncol(CDi)],as.numeric)
      }
    }else{
      if(regr=="lm"|regr=="lrm"|regr=="mlogit"){
        for(v in c(1,3:(2+np+nf))){
          CDi[,v]<-factor(CDi[,v],levels=levels(CD[,v]),exclude=NULL)
        }
      }else if(regr=="rf"){
        for(v in c(3:(2+np+nf))){
          CDi[,v]<-factor(CDi[,v],levels=c(1:(k+1)))
          CDi[,v][is.na(CDi[,v])]<-k+1
        }
        CDi[,1]<-factor(CDi[,1],levels=c(1:k))
      }else{
        #glmnet
        CDi[,1] <- factor(CDi[,1],levels=c(1:k))
        for(v in 3:(2+np+nf)){
          CDi[,v]<-factor(CDi[,v],levels=levels(CD[,v]),exclude=NULL)
        }
      }
      # The last values of CDi must be of type numeric
      # (distributions)
      if (pastDistrib | futureDistrib) {
        CDi[,(2+np+nf+1):ncol(CDi)] <- lapply(CDi[,(2+np+nf+1):ncol(CDi)],as.numeric)
      }
    }
    # Eventually concatenating CDi with COi
    # (the matrix containing the covariates)
    if (all(is.na(CO))==FALSE) {
      # Checking if CO is NOT
      # completely empty
      # Creation of the matrix COi used in 3.3
      COi <- do.call(rbind, replicate(k, as.data.frame(CO[i,]), simplify=FALSE))
      # Concatenating CDi and COi into CDi
      CDi <- cbind(CDi,COi)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }
    # Else, in case CO is empty (i.e. we don't
    # consider any covariate)
    # simply continue with the current CDi
    
    # Eventually concatenating CDi with
    # COtselected_i (the matrix containing the
    # current time-dependent covariates)
    # Checking if COt is NOT completely empty
    if (ncot > 0) {
      COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
      for (d in 1:(ncot/nc)) {
        COtselected_i <- cbind(COtselected_i, COt[i,(j) + (d-1)*nc])
      }
      COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
      # Concatenating CDi and COtselected_i into
      # CDi
      CDi <- cbind(CDi, COtselected_i)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }
    
    
    # Check for missing-values among predictors
    # (i.e. we won't impute any value on the current
    # MD if there is any NA among the VIs)
    if (max(is.na(CDi[1,2:totV]))==0){
      # checking that
      # there is no NA
      # among the current
      # VI (otherwise no
      # data will be
      # imputed for the
      # current NA)
      ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
      
    }
  }
  return(ODi)
}



################################################################################
# Impute value with the chosen regression model

RegrImpute <- function(ODi, CDi, regr, reglog, noise, i, j, k){
  if(regr=="glmnet"){
    pred <- predict(reglog,newx=CDi,type="response")
  }else if (regr == "mlogit") {
    ## Case of MULTINOMIAL REGRESSION MODEL
    
    if(nrow(CDi)!=length(reglog$freq)){
      CDi <- CDi[1:length(reglog$freq),]
    }
    
    pred <- predict(reglog,newdata=CDi)
    # Example of value returned by pred:
    # (Sytematically, the upper line
    # represents the possible categories of
    # the variable (here, k=2, so the
    # possible categories are
    # either 1" or "2"))
    #            1            2
    # 1.000000e+00 2.111739e-22
    #
    # Cumulative pred
    pred <- cumsum(pred)
    
    # We save the potential values of the variable because it
    # could happen that they have the form 1 2 4 5. This way,
    # we are sure to have the real values
    names_saved <- names(pred)
    # Corresponding example value returned
    # by the "cumulative pred":
    # 1 2
    # 1 1
    #
    # Introducing a variance "noise"
    pred <- rnorm(length(pred),pred,noise)
    # Checking that the components of vector
    # "pred" are still included in the
    # interval [0,1]
    pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
    pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
    # Imputation
    alea <- runif(1)
    # Example value returned in "alea":
    # [1] 0.2610005
    #
    sel <- as.numeric(names_saved[which(pred>=alea)])
    # Corresponding example value returned
    # in sel:
    # 1 2
    # 1 2
    #
    
    
  } else if (regr == "lm") {
    
    # Since we are performing a linear
    # regression, each element of CDi are
    # numbers and have to be considered as
    # class "numeric"
    CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))
    
    pred <- predict(reglog, CDi)
    # Introducing a variance "noise"
    pred <- rnorm(length(pred),pred,noise)
    # Rounding pred to the nearest integer
    pred <- round(pred)
    # Restricting pred to its lowest
    # value: 1
    pred <- ifelse(pred<1,1,pred)
    # Restricting pred to its highest
    # value: k
    pred <- ifelse(pred>k,k,pred)
    sel <- pred
    
    
  } else if(regr=="lrm") { # meaning (regr == "lrm")
    ## Case of ORDINAL REGRESSION MODEL
    
    # Since we are performing an ordinal
    # regression, each element of CDi are
    # numbers and have to be considered as
    # class "numeric"
    CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))
    
    pred <- predict(reglog, CDi, type="fitted.ind")
    # Testing if we are in case where k=2
    # (if this is the case, we need to
    # create the second complementary
    # probility by hand since lrm returns
    # only the first probability)
    if (k == 2) {
      pred <- c(pred,(1-pred))
    }
    # Cumulative pred
    pred <- cumsum(pred)
    # Introducing a variance "noise"
    pred <- rnorm(length(pred),pred,noise)
    # Checking that the components of vector
    # "pred" are still included in the
    # interval [0,1]
    pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
    pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
    # Imputation
    # Since we have introduce a noise on the
    # variance, it might occur that "alea"
    # is greater than the greatest value of
    # "pred". We have then to restrict
    # "alea" to the last value of "pred"
    alea <- runif(1)
    if (alea > pred[length(pred)]) {
      alea <- pred[length(pred)]
    }
    sel <- which(pred>=alea)
    
    
    
  }else if(regr=="rf"){ # regr=="rf" randomForest
    #pred <- predict(reglog,newdata=CDi,type="prob")[1,]
    pred<-predict(reglog,data=CDi,predict.all=T)$predictions[1,]
    pred<-factor(pred,levels=c(1:k))
    tab <- table(pred)
    tab <- tab/sum(tab)
    # 
    # # Example of value returned by pred:
    # # (Sytematically, the upper line
    # # represents the possible categories of
    # # the variable (here, k=2, so the
    # # possible categories are
    # # either 1" or "2"))
    # #            1            2
    # # 1.000000e+00 2.111739e-22
    # #
    # # Cumulative pred
    pred <- cumsum(tab)
  
    # Corresponding example value returned
    # by the "cumulative pred":
    # 1 2
    # 1 1
    #
    # Introducing a variance "noise"
    #pred <- rnorm(length(pred),pred,noise)
    # Checking that the components of vector
    # "pred" are still included in the
    # interval [0,1]
    #pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
    #pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
    # Imputation
    alea <- runif(1)
    # Example value returned in "alea":
    # [1] 0.2610005
    #
    sel <- which(pred>=alea)
    # Corresponding example value returned
    # in sel:
    # 1 2
    # 1 2
    #
    
  }else{ # random forest proba fed to a multinomial model
    
    
  }
  
  
  ODi[i,j] <- sel[1]
  return(ODi)
}



################################################################################
# Compute model with the chosen regression model

ComputeModel <- function(CD, regr, tot_VI, np, nf, k,num.trees,min.node.size,max.depth,timing){
  npfi <- np+nf
  # ==>> Manipulation of parameters
  
  # Conversion of CD in a data frame
  CD <- as.data.frame(CD)
  
  # Transformation of the names of the columns of CD
  # (called V1, V2, ..., "Vtot_VI")
  colnames(CD) <- paste("V", 1:ncol(CD), sep = "")
  
  
  
  if (regr == "mlogit") {
    ## Case of MULTINOMIAL REGRESSION MODEL
    
    # Linking the package mlogit
    
    # By default, every column of CD are of class "numeric".
    # Thus, there is no need to convert the columns
    # containing the distribution data to class "numeric".
    # Moreover the class of the covariates columns at the
    # very end are ALREADY set correctly and we don't need
    # to update them.
    # On the other hand, the first columns of CD (1 up to
    # 1+np+nf) have to be of class "factor" (because they
    # are the columns containing our categorical
    # data coming from OD).
    
    # Transformation of the first columns
    # (i.e. the categorical values) of CD (column 1 up to
    # column 1+np+nf) into factor
    if(timing==FALSE){
      CD[,(1:(1+npfi))] <- lapply(CD[,(1:(1+npfi))],factor, levels=c(1:k,NA), exclude=NULL)
      
      
      # We drop the levels of the variables that do not appear. Not doing
      # so create an error later on when the mlogit model is computed
      # (due to the inversion of a matrix having only zeros in a given row)
      for(n in 1:(1+npfi)){
        CD[,n]<-droplevels(CD[,n])
      }
    }else{
      CD[,c(1,3:(2+npfi))] <- lapply(CD[,c(1,3:(2+npfi))],factor, levels=c(1:k,NA), exclude=NULL)
      
      
      # We drop the levels of the variables that do not appear. Not doing
      # so create an error later on when the mlogit model is computed
      # (due to the inversion of a matrix having only zeros in a given row)
      for(n in c(1,3:(2+npfi))){
        CD[,n]<-droplevels(CD[,n])
      }
    }
    
    # Dataframe for mlogit
    NCD <- dfidx(CD, varying=NULL, choice="V1", shape="wide")
    
    
    # Computation of the multinomial model
    reflevel <- as.character(min(levels(CD[,1]),na.rm=T))
    if(tot_VI==1){
      # First case is evaluated aside
      if(reflevel=="1"){
        reglog <- mlogit(V1~0, data=NCD, reflevel="1")
      }else{
        reglog <- mlogit(V1~0, data=NCD, reflevel="2")
      }
    }
    
    if(tot_VI>1){
      # creation of "V2" ... "Vtot_VI" (to use them in the
      # formula)
      factors_character <- paste("V", 2:tot_VI, sep = "")
      # Transformation of this object from character to
      # vector (in order to be able to access its
      # components)
      factors <- as.vector(factors_character)
      # Creation of a specific formula according to the
      # value of tot_VI
      fmla <- as.formula(paste("V1~0|", paste(factors, collapse="+")))
      if(reflevel=="1"){
        reglog <- try(mlogit(fmla, data=NCD, reflevel="1"))
      }else{
        reglog <- try(mlogit(fmla, data=NCD, reflevel="2"))
      }
      if (class(reglog) == "try-error"){
        warning(paste("/!\\ A simpler model was used at some point."))
        if(timing==F){
          if(np>0){
            fmla <- as.formula(paste0("V1~0|V",1+np))
          }else{
            fmla <- as.formula("V1~0|V2")
                               
          }
        }else{
          if(np>0){
            fmla <- as.formula(paste0("V1~0|V",2+np))
          }else{
            fmla <- as.formula("V1~0|V3")
          }
        }
        print(fmla)
        if(reflevel=="1"){
          reglog <- try(mlogit(fmla, data=NCD, reflevel="1"))
        }else{
          reglog <- try(mlogit(fmla, data=NCD, reflevel="2"))
        }        
        
        if(class(reglog)=="try-error"){
          fmla <- as.formula("V1~1")
          if(reflevel=="1"){
            reglog <- mlogit(fmla, data=NCD, reflevel="1")
          }else{
            reglog <- mlogit(fmla, data=NCD, reflevel="2")
          }   
        }
      }
    }
    
    
  } else if (regr == "lm") {
    ## Case of LINEAR REGRESSION MODEL
    
    # Since we are performing a linear regression, each
    # element of CD are numbers and have then to remain as
    # class "numeric" (we don't have to perform some class
    # adjustment as for the case of the creation of a
    # multinomial model).
    
    # Computation of the linear regression model
    if(tot_VI==1){
      reglog <- lm(V1~0, data=CD)
      # first case is evaluated aside
    }
    
    if(tot_VI>1){
      # creation of "V2" ... "Vtot_VI" (to use them in the
      # formula)
      factors_character <- paste("V", 2:tot_VI, sep = "")
      # Transformation of this object from character to
      # vector (in order to be able
      # to access its components)
      factors <- as.vector(factors_character)
      # Creation of a specific formula according to the
      # value of tot_VI
      fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
      reglog <- lm(fmla, data=CD)
    }
    
    
    
    
  } else if(regr=="lrm") { # meaning (regr == "lrm")
    ## Case of ORDINAL REGRESSION MODEL
    
    # Linking to the package rms to use the function "lrm"
    
    # Since we are performing an ordinal regression, each
    # element of CD are numbers and have then to remain as
    # class "numeric" (we don't have to perform some
    # class adjustment as for the case of the creation of a
    # multinomial model).
    
    # Computation of the ordinal model
    if(tot_VI==1){
      reglog <- lrm(V1~0, data=CD)
      # first case is evaluated aside
    }
    
    if(tot_VI>1){
      # creation of "V2" ... "Vtot_VI" (to use them in the
      # formula)
      factors_character <- paste("V", 2:tot_VI, sep = "")
      # Transformation of this object from character to
      # vector (in order to be able
      # to access its components)
      factors <- as.vector(factors_character)
      # Creation of a specific formula according to the
      # value of tot_VI
      fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
      reglog <- lrm(fmla, data=CD)
    }
    
  }else if(regr=="rf"){ # Case regr=="rf" random forest
    if(timing==FALSE){
      CD[,(1:(1+npfi))] <- lapply(CD[,(1:(1+npfi))],factor, levels=c(1:k))
      for(v in 2:(1+npfi)){
        CD[,v]<-factor(CD[,v],levels=c(1:(k+1)))
        CD[,v][is.na(CD[,v])]<-k+1
      }
      CD[,1]<-factor(CD[,1],levels=c(1:k))
      CD[,1]<-droplevels(CD[,1])
      CD<-CD[!is.na(CD[,1]),]
    }else{
      CD[,c(1,3:(2+npfi))] <- lapply(CD[,c(1,3:(2+npfi))],factor, levels=c(1:k))
      for(v in 3:(2+npfi)){
        CD[,v]<-factor(CD[,v],levels=c(1:(k+1)))
        CD[,v][is.na(CD[,v])]<-k+1
      }
      CD[,1]<-factor(CD[,1],levels=c(1:k))
      CD[,1]<-droplevels(CD[,1])
      CD<-CD[!is.na(CD[,1]),]
    }
    
    factors_character <- paste("V", 2:tot_VI, sep = "")
    # Transformation of this object from character to
    # vector (in order to be able
    # to access its components)
    factors <- as.vector(factors_character)
    # Creation of a specific formula according to the
    # value of tot_VI
    fmla <- as.formula(paste("V1~", paste(factors, collapse="+")))
    #reglog <- randomForest(fmla, data=CD,ntree=100)
    reglog<-ranger(fmla,data=CD,num.trees=num.trees,min.node.size=min.node.size,
                   max.depth=max.depth)
    
  }else{
    }
  # else if(regr=="glmnet"){# random forest proba fed to a multinomial model
  #   
  #   # Linking the package mlogit
  #   
  #   # By default, every column of CD are of class "numeric".
  #   # Thus, there is no need to convert the columns
  #   # containing the distribution data to class "numeric".
  #   # Moreover the class of the covariates columns at the
  #   # very end are ALREADY set correctly and we don't need
  #   # to update them.
  #   # On the other hand, the first columns of CD (1 up to
  #   # 1+np+nf) have to be of class "factor" (because they
  #   # are the columns containing our categorical
  #   # data coming from OD).
  #   
  #   # Transformation of the first columns
  #   # (i.e. the categorical values) of CD (column 1 up to
  #   # column 1+np+nf) into factor
  #   
  #   CD <- CD[!is.na(CD$V1),]
  #   if(timing==FALSE){
  #     CD[,1] <- factor(CD[,1],levels=c(1:k))
  #     CD[,(2:(1+npfi))] <- lapply(CD[,(2:(1+npfi))],factor, levels=c(1:k,NA), exclude=NULL)
  #     
  #     
  #     # We drop the levels of the variables that do not appear. Not doing
  #     # so create an error later on when the mlogit model is computed
  #     # (due to the inversion of a matrix having only zeros in a given row)
  #     for(n in 2:(1+npfi)){
  #       CD[,n]<-droplevels(CD[,n])
  #     }
  #   }else{
  #     CD[,1] <- factor(CD[,1],levels=c(1:k))
  #     CD[,c(3:(2+npfi))] <- lapply(CD[,c(3:(2+npfi))],factor, levels=c(1:k,NA), exclude=NULL)
  #     
  #     
  #     # We drop the levels of the variables that do not appear. Not doing
  #     # so create an error later on when the mlogit model is computed
  #     # (due to the inversion of a matrix having only zeros in a given row)
  #     for(n in c(3:(2+npfi))){
  #       CD[,n]<-droplevels(CD[,n])
  #     }
  #   }
  #   
  #   # Computation of the multinomial model
  #   if(tot_VI==1){
  #     # First case is evaluated aside
  #     reglog <- mlogit(V1~0, data=NCD, reflevel="1")
  #   }
  #   
  #   if(tot_VI>1){
  #     # creation of "V2" ... "Vtot_VI" (to use them in the
  #     # formula)
  #     y <- CD[,1]
  #     x <- model.matrix(V1~.,CD)[,-1]
  #     #reglog <- mlogit(fmla, data=NCD, reflevel="1")
  #     reglog <- try(glmnet(x=x,y=y,lambda=0,family="multinomial"))
  #     if (class(reglog) == "try-error"){
  #       warning(paste("/!\\ A simpler model was used at some point."))
  #       if(timing==F){
  #         if(np>0){
  #           x<- CD[,1+np]
  #         }else{
  #           x <- CD[,2]
  #           
  #         }
  #       }else{
  #         if(np>0){
  #           x <- CD[,2+np]
  #         }else{
  #           x <- CD[,3]
  #         }
  #       }
  #       reglog <- try(glmnet(x=x,y=y,lambda=0))
  #       if(class(reglog)=="try-error"){
  #         fmla <- as.formula("V1~1")
  #         reglog <- mlogit(fmla,data=NCD,reflevel="1")
  #       }
  #     }
  #   }
  # }
  
  return(list(reglog, CD))
}


