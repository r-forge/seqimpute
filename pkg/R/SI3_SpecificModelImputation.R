# 3. Imputation using a specific model -----------------------------------------


################################################################################
#' Impute data using a specific model
#'
#' @export
ModelImputation <- function(OD, CO, COt, ODi, MaxGap, totV, totVi, regr, nc, np, nf, nr, ncot, COtsample, pastDistrib, futureDistrib, k, available, REFORD_L, noise,num.trees,min.node.size,max.depth){
  
  for (order in 1:MaxGap){ 
    # /!\ "order" corresponds to the
    # values of the components of ORDER (i.e. the number
    # of the iteration, the order in which the
    # values are going to be imputed)
    # 3.1. Building of the data matrix CD for the computation of the model ----------------------------
    CD_shift <- CDCompute(CO, OD, COt, MaxGap, order, np, nc, nr, nf, COtsample, pastDistrib, futureDistrib, ncot, nco, k)
    # 3.2. Computation of the model (Dealing with the LOCATIONS of imputation)-------------------------
    CD_shift[c("reglog","CD")] <- ComputeModel(CD_shift$CD, regr, totV, np+nf, k,num.trees,min.node.size,max.depth)
    

    # 3.3. Imputation using the just created model (Dealing with the actual VALUES to impute) ---------
    ODi <- CreatedModelImputation(order, CO, CD_shift$CD, COt, OD, ODi, pastDistrib, futureDistrib, available, REFORD_L, ncot, nc, np, nf, k, totV, regr, CD_shift$reglog, noise, CD_shift$shift, MaxGap)

  
  } 
  return(ODi)
}



################################################################################
#' Building of the data matrix CD for the computation of the model
#'
#' @export
CDCompute <- function(CO, OD, COt, MaxGap, order, np, nc, nr, nf, COtsample, pastDistrib, futureDistrib, ncot, nco, k){
  # Building of a data matrix for the computation of the
  # model
  ud <- nc-(MaxGap-order+np+nf)    # number of usable data

  
  # for each row of OD
  frameSize <- MaxGap-order+np+nf+1 # size of the current
  # mobile caracteristic frame (that
  # changes according to
  # "order") which is equal to
  # nc-ud+1
  # Structure and building of the data matrix CD
  # The first column of CD is the dependent variable (VD,
  # response variable)
  # The following columns are the independent variables
  # (VIs, explanatory variables) coming from the past
  # (np>0) or the future (nf>0) ordered by time and the
  # distribution of the possible values (i.e. all possible
  # categorical variables numbered from 1 to k and of
  # course the value NA) respectively Before and After the
  # NA to impute.
  #
  #           VD   Past VIs   Future VIS    Past distribution    Future distribution
  #
  #        /                                                                           \
  #       |                                                                            |
  # CD =  |  CD      CDp        CPf               CDdb                  CDda          |
  #       \                                                                          /
  #        \                                                                       /
  #
  # We are then going to create a specific CD according to
  # the value of np and nf

  CD <- matrix(NA,nrow=nr*ud,ncol=1) # initialization of
  # the current very left part of
  # the predictive model matrix
  # ("matrice de modele de
  # prediction") with NA
  # everywhere (/!\ for each
  # "order", we are going to
  # build such a CD)

  # Dealing with the change of shape of the prediction
  # frame (according to whether the imputed data is
  # located at the beginning (left) of a gap or at the end
  # (right)).
  # The purpose of this if statement is to detect if we
  # have to insert a shift (to jump at the end of the gap)
  # or not
  if ( (np > 0 & nf > 0) & ( (MaxGap%%2==0 & order%%2==0) | (MaxGap%%2!=0 & order%%2!=0) )){
    shift <- MaxGap - order      # jumping at the end of
    # the gap
  } else {
    shift <- 0           # no shift is needed (note that
    # no shift is needed for the
    # case of model 2 (only past)
    # and model 3 (only future))
  }
  
  iter <- 1                # initialisation of the number
  # of iterations of the
  # following for loops
  # Only PAST
  if (np>0 & nf==0) { 
    CD <- PastVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, COtsample, COt, pastDistrib, futureDistrib, k)
    # Only FUTURE
  } else if (np==0 & nf>0) {
    # only FUTURE VIs do exist
    CD <- FutureVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, COtsample, COt, pastDistrib, futureDistrib, k, nf)
    
    # PAST and FUTURE
  } else {
    # meaning np>0 and nf>0 and that, thus,
    # PAST as well as FUTURE VIs do exist
    CD <- PastFutureVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, COtsample, COt, pastDistrib, futureDistrib, k, nf, shift)
  }
  
  CD_shift <- list()
  CD_shift[c("CD","shift")] <- list(CD, shift)
  return(CD_shift)
}



################################################################################
#' Imputation using the just created model (Dealing with the actual VALUES to impute) 
#'
#' @export
CreatedModelImputation <- function(order, CO, CD, COt, OD, ODi, pastDistrib, futureDistrib, available, REFORD_L, ncot, nc, np, nf, k, totV, regr, reglog, noise, shift, MaxGap){
  # Structure and building of the data matrix CDi
  # The first column of CDi is the dependent variable (VD,
  # response variable) that we have to implement during
  # the current iteration (i.e. automatically a NA)
  # The following columns are the corresponding
  # independent variables (VIs, explanatory variables)
  # coming from the past (np>0) or the future (nf>0)
  # (ordered by time and corresponding to the current
  # predictive pattern) and the distribution of the
  # possible values (i.e. all possible categorical
  # variables numbered from 1 to k and of course
  # the value NA) respectively Before and After the NA to
  # impute.
  # (The fact that every lines of CDi are identical is
  # related to the working of the function mlogit that has
  # to have as much lines in CDi as there are categories
  # of the variable (see the parameter "k") --> so, CDi is
  # composed of k identical lines)
  #
  #            VD   Past VIs   Future VIS    Past distribution     Future distribution
  #
  #         /                                                                             \
  #        |                                                                              |
  # CDi =  |  CDi    CDpi       CPfi              CDdbi                   CDdai           |
  #        \                                                                             /
  #         \                                                                          /
  #
  # We are then going to create a specific CDi according
  # to the value of np and nf
  
  
  # Analysing the value of parameter available
  if (available==TRUE){   # we take the previously imputed
    # data into account
    LOOKUP <- ODi
  } else { # that is available == FALSE and thus we
    # don't take the previously imputed
    # data into account
    LOOKUP <- OD
  }
  
  
  # Assigning the current "REFORD_order" matrix to the
  # variable matrix REFORD
  # (according to the current value of "order")
  REFORD <- as.matrix(REFORD_L[[order]])
  if (ncol(REFORD) == 1) {
    REFORD <- t(REFORD)
  }
  nr_REFORD <- nrow(REFORD)
  
  
  if (np>0 & nf==0) { # only PAST VIs do exist
    ODi <- ODiImputePAST(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise)
    
        #---------------------------------------------------------------------------------------------
  } else if (np==0 & nf>0) {  # only FUTURE VIs do exist
    ODi <- ODiImputeFUTURE(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise)
  } else { # meaning np>0 and nf>0 and that,
    # thus, PAST as well as FUTURE VIs
    # do exist
    ODi <- ODiImputePF(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise, shift, MaxGap, order)
    
  }
  return(ODi)
}



################################################################################
#' Imputation where only FUTURE VIs exist
#'
#' @export
ODiImputeFUTURE <- function(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise){
  for (u in 1:nr_REFORD){
    i <- REFORD[u,1]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- REFORD[u,2]
    # taking out the second coordinate
    # (column number in ORDER) from REFORD
    
    CDi <- matrix(NA, nrow=k, ncol=1)
    
    # Matrix for future values
    vect <- LOOKUP[i,(j+1):(j+nf)]
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
    CDi <- cbind(CDi, CDfi)
    
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
      
      COtselected_i <- COt[i, j + ((1:(ncot/nc))-1)*nc]

      COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
      # Concatenating CDi and COtselected_i into
      # CDi
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

