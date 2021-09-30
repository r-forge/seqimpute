# 0. File containing preliminary functions used for checking the correctness and usability of the data.
#
#
################################################################################
#' In case of a factor dataset OD:
#' RECODING of OD with numbers "1", "2", etc. instead of its "words"
#'
#' @export
preliminaryChecks <- function(OD, CO, COt, np, nf, nfi, npt, k, pastDistrib, futureDistrib){
  
  dataOD <- list()
  # Updating the number of columns of CO
  dataOD["nco"] <- emptyColUpdate(CO)
  #
  # Updating the number of columns of COt
  dataOD["ncot"] <- emptyColUpdate(COt)
  
  # Deleting entire rows of OD filled only with NAs
  # (and deleting corresponding lines in CO and COt)
  dataOD[c("OD", "CO", "COt", "rowsNA")]  <- deleteNaRows(OD,CO,COt)    
  dataOD[c("nr", "nc", "COtsample", "totV", "totVt", "totVi")] <- defineNbVariables(dataOD$OD, dataOD$COt, dataOD$ncot, np, dataOD$nco, nf, nfi, npt, k, pastDistrib, futureDistrib)
  dataOD[c("OD", "ODi", "ODClass", "ODlevels")] <- factorAndMatrixChecks(dataOD$OD, dataOD$nc, k)
  
  return(dataOD)
}


################################################################################
#' Initial tests and manipulations on parameters
#'
#' @export
InitCorectControl <- function(regr, ODClass, OD, nr, nc, k, np, nf, nco, ncot, nfi, npt, pastDistrib, futureDistrib, totV, totVi, totVt, noise){

  # 0.1 Testing "regr" effectively either "mlogit", "lm" or "lrm" -----------------------------------------------------------------------------------------------
  
  if ( (regr != "mlogit") & (regr != "lm") & (regr != "lrm") & (regr!="rf") &(regr!="glmnet")) {
    stop("/!\\ regr defines the type of regression model you want to use.
               It has to be either assigned to character 'mlogit' (for multinomial
               regression),'lm' (for linear regression) or 'lrm' (for ordinal
               regression)")
  }
  
  
  
  
  
  # 0.2 Testing the class of the variables of the original dataset OD -------------------------------------------------------------------------------------------
  if ( (ODClass != "factor") & (ODClass != "numeric") ) {
    stop("/!\\ The class of the variables contained in your original dataset
           should be either 'factor' or 'numeric'")
  }
  
  
  
  
  # 0.3 Testing effectively exactly k possible categories of the variable (in case we consider categorical variables) -------------------------------------------
  if (ODClass == "factor") {
    for (i in 1:nr) {
      for (j in 1:nc) {
        if ( is.na(OD[i,j])==FALSE & (OD[i,j]<=0 | OD[i,j]>k) ) {
          stop("/!\\ Your Original dataset doesn't contain the right
                       number of k categories of the variable")
        } else {
          next
        }
      }
    }
  }
  
  
  
  
  
  
  
  # 0.4 Eventually discarding the consideration of pastDistrib and futureDistrib --------------------------------------------------------------------------------
  # Making sure that pastDistrib and futureDistrib are set to FALSE by default
  # in case OD is made of "numeric" variables
  if (ODClass=="numeric") {
    pastDistrib <- FALSE
    futureDistrib <- FALSE
    # Update of the totV variables
    # (which become then smaller since
    # we don't consider pastDistrib and futureDistrib)
    totV <- 1+np+nf+nco+(ncot/nc)
    totVi <- 1+nfi+nco+(ncot/nc)
    totVt <- 1+npt+nco+(ncot/nc)
  }
  
  
  
  
  
  
  
  
  
  # 0.5 Testing not np<0, nor nf<0 nor both ==0 -----------------------------------------------------------------------------------------------------------------
  # No VIs is not possible (we should have at least np>0 or nf>0)
  if(np==0 & nf==0)
    stop("/!\\ We can't have np as well as nf equal to '0' at the same
               time")
  # Negative value for np as well as nf raises an error
  if(np<0 | nf<0)
    stop("/!\\ np and nf can't be negative numbers")
  
  
  
  
  
  
  # 0.6 Testing np and nf are not chosen too large --------------------------------------------------------------------------------------------------------------
  # if (np+1+nf>nc)
  #   stop("/!\\ You have to choose lower value for np and nf. Your selected
  #              np and nf are too large to fit the dimensions of your data matrix")
  # 
  # 
  
  
  
  # 0.7 Testing not nfi<0, nor npt<0 ----------------------------------------------------------------------------------------------------------------------------
  # Negative value for nfi as well as npt raises an error
  if(nfi<0 | npt<0)
    stop("/!\\ nfi and npt can't be negative numbers")
  # In case nfi = 0 and/or npt = 0, the imputation of the initial gaps and/or
  # the terminal gaps respectively is simply omitted
  
  
  
  
  
  
  
  
  # 0.8 Testing the right construction of COt -------------------------------------------------------------------------------------------------------------------
  # Since COt contains the time-dependent covariates, each of them is
  # represented by a submatrix that has the same number of columns as OD.
  # The total number of columns of COt is necessarily a multiple of the number
  # of column of OD
  if (ncot%%nc != 0) {
    stop("/!\\ Be sure to have understood the notion of time-dependent covariates. Each time-dependent covariates contained in COt has to have the same number of columns as OD.")
  }
  
  
  # 0.9 Taking the absolute value of the parameter "noise" -----------------------------------------------------------------------------------------------------
  # Since "noise" is the variance of the elements of the final vector pred, it
  # can't be negative
  noise <- abs(noise)

  return(list(pastDistrib, futureDistrib, totV, totVi, totVt, noise))
}


################################################################################
#' Update the number of columns of a matrix
#'
#' @export
emptyColUpdate <- function(x) {
  if (all(is.na(x))==FALSE) { # if CO is NOT completely empty
    return(ncol(x))
  }
  # else, in case CO is completely empty
  # then, nco has to be set to "0"
  return(0)
}
  

################################################################################
#' Deleting entire rows of OD filled only with NAs
#' (and deleting corresponding lines in CO and COt)
#'
#' @export
deleteNaRows <- function(OD, CO, COt) {
  rowsNA <- c()
  for(i in 1:nrow(OD)) {
    if (all(is.na(OD[i,]))) {
      rowsNA <- c(rowsNA,i)
    }
  }
  
  if(length(rowsNA)>0){
    OD <- OD[-rowsNA,]
    if (all(is.na(CO))==FALSE) { # Checking if CO is NOT completely
      # empty and updating the covariate matrix CO as well!
      CO <- CO[-rowsNA,]
    }
    if (all(is.na(COt))==FALSE) { # Checking if COt is NOT completely
      # empty and updating the covariate matrix COt as well!
      COt <- COt[-rowsNA,]
    }
    warning(paste("/!\\ Row number",paste(rowsNA,collapse=", "),"was not imputed
                            because it only consisted of NAs."),sep='')
  } 
  return(list(OD, CO, COt, rowsNA))
}


################################################################################
#' Definition of the number of rows and columns in OD
#' (And eventually update of nr)
#'
#' @export
defineNbVariables <- function(OD, COt, ncot, np, nco, nf, nfi, npt, k, pastDistrib, futureDistrib) {
  
  nr <- nrow(OD)
  nc <- ncol(OD)
  
  COtsample <- vector()
  if (ncot > 0) {
    # Creaion of a sample of COt
    # (taking a unique sample of each time-dependent covariates and coercing
    # them into the dataframe COtsample in order to create COtselected (with
    # some correct corresponding classes on each column!) in the
    # parts 3.1)
    COtsample <- as.data.frame(matrix(nrow=nr,ncol=0))
    
    for (d in 1:(ncot/nc)) {
      COtsample <- cbind(COtsample,COt[,1 + (d-1)*nc])
    }
  }
  
  # Total number of variables in the imputation model
  totV <- 1+np+nf+nco+(ncot/nc)
  totVi <- 1+nfi+nco+(ncot/nc)
  totVt <- 1+npt+nco+(ncot/nc)
  if (pastDistrib) {
    totV <- totV + k
    totVt <- totVt + k
  }
  if (futureDistrib) {
    totV <- totV + k
    totVi <- totVi + k
  }
  return(list(nr, nc, COtsample, totV, totVt, totVi))
}


################################################################################
#' In case of a factor dataset OD:
#' RECODING of OD with numbers "1", "2", etc. instead of its "words"
#'
#' @export
factorAndMatrixChecks <- function(OD, nc, k) { 
  ODClass <- class(OD[1,1])
  #*************************************
  ODlevels <- vector()
  if (ODClass == "factor") {
    ODlevels <- levels(OD[1,1])
    OD <- as.data.frame( sapply(OD, mapvalues, from = ODlevels,
                                to = as.character(as.vector(1:length(ODlevels)))) )
  }
  #*************************************
  
  
  # Making sure that OD is a matrix and not a data frame
  # /!\ Using simply "OD <- data.matrix(OD)"
  # may convert the components of OD to a different number
  # The data.matrix() function converts factors to numbers by using their
  # internal codes.
  # That's why they're listed as factors in the data frame and have different
  # values after using data.matrix(). To create a numeric matrix in this
  # situation, we use rather the syntax:
  #********************************************
  OD <- apply(as.matrix(OD), 2, as.numeric)
  #********************************************
  # When using as.matrix(), factors become strings. Using apply() will convert
  # everything to numeric without losing the matrix structure.
  
  
  
  # Converting the columns of OD to factor
  #
  # /!\ With:
  # for (j in 1:nc) {
  #     OD[,j] <- as.factor(OD[,j])
  # }
  # Copying OD in ODi (the matrix in which we will proceed to the
  # imputations ODi <- OD
  # The values in certain components of OD change because as we convert into
  # factors, the column containing only "1" or only "2" are then substituted
  # into only "1" (because only one unique state (either "1" or "2") has been
  # identified in the column)
  # Here is the code to use to preserve original values in a variable turned
  # into a factor:
  ODi <- OD
  #
  # In case OD is constituted of factor variables, we make sure that the
  # variables of ODi are considered as factor ranging from "1" to "k"
  if (ODClass == "factor") {
    for (j in 1:nc){
      ODi[,j] <- factor(x = OD[,j], levels = c(1:k))
    }
  }
  return(list(OD, ODi, ODClass, ODlevels))
  
}

