#' Numbering NAs and types of gaps among a dataset
#'
#' \code{seqQuickLook.R} is a function aimed at providing a quick overview of the
#' frequency of the NAs and the number and size of the different types of gaps
#' spread in the original dataset \code{OD}. The user should run this function before
#' running the main function \code{seqimpute.R} in order to identify a judicious choice
#' for the values of \code{np} and \code{nf}.
#'
#' @param OD \code{matrix} object containing sequences of a variable with missing data (coded as \code{NA}).
#' @param np \code{numeric} object corresponding to the number of previous observations in the imputation model of the internal gaps (default \code{1}).
#' @param nf \code{numeric} object corresponding to the number of future observations in the imputation model of the internal gaps (default \code{0}).
#'
#' @author Andre Berchtold <andre.berchtold@@unil.ch>
#'
#' @return It returns a  \code{data.frame} object that summarizes for each type of gaps 
#' (Internal Gaps, Initial Gaps, Terminal Gaps, LEFT-
#' hand side SLG, RIGHT-hand side SLG, Both-hand side SLG),
#' the minimum length, the maximum length, the total number of gaps and 
#' the total number of NAs induced.
#'
#' @examples
#' data(OD)
#'
#' seqQuickLook(OD=OD, np=1, nf=0)
#'
#' @export


seqQuickLook <- function(OD, np=1, nf=0) {







    # 1. Initial tests on parameters --------------------------------------------------------------------------------------------------------------------


    # 1.1 Testing the class of the variables of the original dataset OD ---------------------------------------------------------------------------------
    ODClass <- class(OD[1,1])
    if ( (ODClass != "factor") & (ODClass != "numeric") ) {
        stop("/!\\ The class of the variables contained in your original dataset
         should be either 'factor' or 'numeric'")
    }


    if (ODClass == "factor") {
      k <- length(sort(unique(as.vector(as.matrix(OD)))))
    }else{
      k <- max(OD)
    }
    
    dataOD <- list()
    dataOD[c("OD", "CO", "COt", "rowsNA")]  <- deleteNaRows(OD,CO,COt)   
    
    
    
    # Naming the number of rows and columns of the original dataset OD
    nr <- nrow(dataOD$OD)
    nc <- ncol(dataOD$OD)
    
    
    
    
    
    # 1. Analysis of OD and creation of matrices ORDER, ORDER2 and ORDER3 -----------------------------------------------------------------------------------------
    dataOD[c("MaxInitGapSize", "InitGapSize",  "MaxTermGapSize", "TermGapSize", "MaxGap", "ORDER", "ORDER2", "ORDER3")] <- OrderCreation(dataOD$OD, nr, nc)
    
    
    
    # 2. Computation of the order of imputation of each MD (i.e. updating of matrix ORDER) --------------------------------------------------------------------
    if (max(dataOD$ORDER)!=0) {
      dataOD[c("ORDERSLGLeft", "ORDERSLGRight", "ORDERSLGBoth", "LongGap", "MaxGap", "REFORD_L", "ORDER")] <- ImputeOrderComputation(dataOD$ORDER, dataOD$ORDER3, dataOD$MaxGap, np, nf, nr, nc)
    }
    

    MDGapsChart <- matrix(0,6,4)
    MDGapsChart <- as.data.frame(MDGapsChart)
    rownames(MDGapsChart) <- c("Internal Gaps","Initial Gaps","Terminal Gaps","LEFT-hand side SLG","RIGHT-hand side SLG","BOTH-hand side SLG")
    colnames(MDGapsChart) <- c("MinGapSize","MaxGapSize","numbOfGaps","sumNAGaps")
    
    if(max(dataOD$ORDER)>0){
      MDGapsChart[1,1] <- min(dataOD$ORDER[dataOD$ORDER>0])
      MDGapsChart[1,2] <- max(dataOD$ORDER)
    }
    # Number of Internal Gaps
    # Transforming dataOD$ORDER in a single row vector
    dataOD$ORDERVect <- as.vector(t(dataOD$ORDER))
    # Transforming this single row vector into class character
    dataOD$ORDERVectChar <- paste(dataOD$ORDERVect,collapse="")
    # Identifying the patterns "0 1" (this is the signature look of an internal gap
    # (it always indicates the beginning of an internal gap!))
    MDGapsChart[1,3] <- str_count(dataOD$ORDERVectChar,pattern="01")
    MDGapsChart[1,4] <- sum(dataOD$ORDER>0)
    
    
    if(max(dataOD$InitGapSize)>0){
      MDGapsChart[2,1] <- min(dataOD$InitGapSize[dataOD$InitGapSize>0])
      MDGapsChart[2,2] <- max(dataOD$InitGapSize)
      MDGapsChart[2,3] <- sum(table(dataOD$InitGapSize[dataOD$InitGapSize>0]))
      MDGapsChart[2,4] <- sum(dataOD$InitGapSize)
      
    }
    
    if(max(dataOD$TermGapSize)>0){
      MDGapsChart[3,1] <- min(dataOD$TermGapSize[dataOD$TermGapSize>0])
      MDGapsChart[3,2] <- max(dataOD$TermGapSize)
      MDGapsChart[3,3] <- sum(table(dataOD$TermGapSize[dataOD$TermGapSize>0]))
      MDGapsChart[3,4] <- sum(dataOD$TermGapSize)
    }
    
    LeftGap <- rowSums(dataOD$ORDERSLGLeft>0)
    if(max(LeftGap)>0){
      MDGapsChart[4,1] <- min(LeftGap[LeftGap>0])
      MDGapsChart[4,2] <- max(LeftGap)
      MDGapsChart[4,3] <- sum(table(LeftGap[LeftGap>0]))
      MDGapsChart[4,4] <- sum(LeftGap)
    }
    
    RightGap <- rowSums(dataOD$ORDERSLGRight>0)
    if(max(RightGap)>0){
      MDGapsChart[5,1] <- min(RightGap[RightGap>0])
      MDGapsChart[5,2] <- max(RightGap)
      MDGapsChart[5,3] <- sum(table(RightGap[RightGap>0]))
      MDGapsChart[5,4] <- sum(RightGap)
    }
    
    BothGap <- rowSums(dataOD$ORDERSLGBoth>0)
    if(max(BothGap)>0){
      MDGapsChart[6,1] <- min(BothGap[BothGap>0])
      MDGapsChart[6,2] <- max(BothGap)
      MDGapsChart[6,3] <- sum(table(BothGap[BothGap>0]))
      MDGapsChart[6,4] <- sum(BothGap)
    }
    
    return(MDGapsChart)
    

}
