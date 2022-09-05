#' Imputation of missing data in sequence analysis
#'
#' Multiple imputation of missing data present in a dataset through the prediction based
#' on either a multinomial or a random forest regression model.
#' Covariates and time-dependant covariates can be included in the model.
#' The prediction of the missing values is based on the theory of Prof. Brendan
#' Halpin. It considers a various amount of surrounding available information to
#' perform the prediction process.
#' In fact, we can among others specify \code{np} (the number of past variables
#' taken into account) and \code{nf} (the number of future information taken
#' into account).
#'
#' @details The imputation process is divided into several steps. According to the location of the gaps of NA among the original dataset, we have defined 5 types of gaps:
#'
#' - Internal Gaps (simple usual gaps)
#'
#' - Initial Gaps (gaps situated at the very beginning of a sequence)
#'
#' - Terminal Gaps (gaps situaed at the very end of a sequence)
#'
#' - Left-hand side SLG (Specially Located Gaps) (gaps of which the beginning location is included in the interval \code{[0,np]}
#' but the ending location is not included in the interval \code{[ncol(OD)-nf,ncol(OD)]})
#'
#' - Right-hand side SLG (Specially Located Gaps) (gaps of which the ending location is included in the interval \code{[ncol(OD)-nf,ncol(OD)]} 
#' but the beginning location is not included in the interval \code{[0,np]})
#' 
#' - Both-hand side SLG (Specially Located Gaps) (gaps of which the beginning location is included in the interval \code{[0,np]}
#' and the ending location is included in the interval \code{[ncol(OD)-nf,ncol(OD)]} )
#'
#' Order of imputation of the gaps types:
#'     1. Internal Gaps
#'     2. Initial Gaps
#'     3. Terminal Gaps
#'     4. Left-hand side SLG
#'     5. Right-hand side SLG
#'     6. Both-hand side SLG
#'
#' @param OD either a data frame containing sequences of a multinomial variable with missing data (coded as \code{NA}) or
#' a state sequence object built with the TraMineR package
#' @param regr a character specifying the imputation method. If \code{regr="multinom"}, multinomial models are used,
#' while if \code{regr="rf"}, random forest models are used.
#' @param np number of previous observations in the imputation model of the internal gaps.
#' @param nf number of future observations in the imputation model of the internal gaps.
#' @param nfi number of future observations in the imputation model of the initial gaps.
#' @param npt number of previous observations in the imputation model of the terminal gaps.
#' @param available a logical value allowing the user to choose whether to consider the already imputed data in the predictive model (\code{available = TRUE}) or not (\code{available = FALSE}).
#' @param CO a data frame containing some covariates among which the user can choose in order to specify his model more accurately.
#' @param COt a data frame object containing some time-dependent covariates that help specifying the predictive model more accurately.
#' @param pastDistrib a logical indicating if the past distribution should be used as predictor in the imputation model.
#' @param futureDistrib a logical indicating if the futur distribution should be used as predictor in the imputation model.
#' @param mi number of multiple imputations  (default: \code{1}).
#' @param mice.return a logical indicating whether an object of class \code{mids}, that can be directly used by the \code{mice} package, should be returned
#' by the algorithm. By default, a data frame with the imputed datasets stacked vertically is returned.
#' @param include logical. If a dataframe is returned (\code{mice.return = FALSE}), indicates if the original
#' dataset should be included or not. This parameter does not apply if \code{mice.return=TRUE}.
#' @param noise \code{numeric} object adding a noise on the predicted variable \code{pred} determined by the multinomial model 
#' (by introducing a variance \code{noise} for each components of the vector \code{pred}) (the user can choose any value for \code{noise}, but we recommend to choose a rather relatively small value situated in the interval \code{[0.005-0.03]}).
#' @param ParExec logical. If \code{TRUE}, the multiple imputations are run in parallell. This allows faster run time depending of how many core the processor has. 
#' @param ncores integer. Number of cores to be used for the parallel computation. If no value is set for this parameter, the number of cores will be set
#' to the maximum number of CPU cores minus 1. 
#' @param SetRNGSeed an integer that is used to set the seed in the case of parallel computation. Note that setting \code{set.seed()} alone before the seqimpute function won't work in case
#' of parallel computation.
#' @param num.trees random forest parameter setting the number of trees of each random forest model.
#' @param min.node.size random forest parameter setting the minimum node size for each random forest model.
#' @param max.depth random forest parameter setting the maximal depth tree for each random forest model.
#' @param verbose logical. If \code{TRUE}, seqimpute will print history and warnings on console. Use \code{verbose=FALSE} for silent computation.
#' 
#' @author Andre Berchtold <andre.berchtold@@unil.ch> Kevin Emery Anthony Guinchard Kamyar Taher
#'
#' @return Returns either an S3 object of class \code{mids} if \code{mice.return = TRUE}
#' or a dataframe, where the imputed dataset are stacked vertically. In the second case,
#' two columns are added: \code{.imp} integer that refers to the imputation number
#' (0 corresponding to the original dataset if \code{include=TRUE}) and \code{.id} character corresponding to
#' the rownames of the dataset to impute. 
#'
#' @examples
#' # Default single imputation
#' RESULT <- seqimpute(OD=OD, np=1, nf=0, nfi=1, npt=1, mi=2)
#' 
#' # Seqimpute used with parallelisation
#' RESULT <- seqimpute(OD=OD, np=1, nf=1, nfi=1, npt=1, mi=2, ParExec=TRUE, SetRNGSeed=17,ncores=2)
#' 
#' @references HALPIN, Brendan (2012). Multiple imputation for life-course sequence data. Working Paper WP2012-01, Department of Sociology, University of Limerick. http://hdl.handle.net/10344/3639.
#' @references HALPIN, Brendan (2013). Imputing sequence data: Extensions to initial and terminal gaps, Stata's. Working Paper WP2013-01, Department of Sociology, University of Limerick. http://hdl.handle.net/10344/3620
#'
#' 
#' @export
seqimpute <- function(OD, regr="multinom", np=1, nf=0, nfi=1, npt=1,
                      available=TRUE, CO=matrix(NA,nrow=1,ncol=1),
                      COt=matrix(NA,nrow=1,ncol=1), pastDistrib=FALSE,
                      futureDistrib=FALSE, mi=1, mice.return=FALSE, include=FALSE, noise=0, SetRNGSeed=FALSE, ParExec=FALSE,ncores=NULL
                      ,num.trees=10,min.node.size=NULL,max.depth=NULL,verbose=TRUE) {
  
 
  #***************************************************************************
  if(sum(is.na(OD))==0){
    if(verbose==T){
      message("This dataset has no missing values!")
    }
    return(OD)
    
  }

  #k <- n_distinct(data.frame(newcol = c(t(OD)), stringsAsFactors=FALSE),na.rm = TRUE)
  
  rownamesDataset <- rownames(OD)
  nrowsDataset <- nrow(OD)
  
  if(mice.return==TRUE){
    include <- TRUE
  }
  
  
  # 0. Initial tests and manipulations on parameters ------------------------------------------------------------------------------------------------------------
  dataOD <- preliminaryChecks(OD, CO, COt, np, nf, nfi, npt, pastDistrib, futureDistrib)
  dataOD[c("pastDistrib", "futureDistrib", "totV", "totVi", "totVt", "noise")] <- InitCorectControl(regr, dataOD$ODClass, dataOD$OD, dataOD$nr, dataOD$nc, dataOD$k, np, nf, dataOD$nco, dataOD$ncot, nfi, npt,  pastDistrib, futureDistrib, dataOD$totV, dataOD$totVi, dataOD$totVt, noise)
  
  
  # 1. Analysis of OD and creation of matrices ORDER, ORDER2 and ORDER3 -----------------------------------------------------------------------------------------
  dataOD[c("MaxInitGapSize", "InitGapSize",  "MaxTermGapSize", "TermGapSize", "MaxGap", "ORDER", "ORDER2", "ORDER3")] <- OrderCreation(dataOD$OD, dataOD$nr, dataOD$nc)
  
  
  
  # 2. Computation of the order of imputation of each MD (i.e. updating of matrix ORDER) --------------------------------------------------------------------
  if (max(dataOD$ORDER)!=0) {
    dataOD[c("ORDERSLGLeft", "ORDERSLGRight", "ORDERSLGBoth", "LongGap", "MaxGap", "REFORD_L", "ORDER")] <- ImputeOrderComputation(dataOD$ORDER, dataOD$ORDER3, dataOD$MaxGap, np, nf, dataOD$nr, dataOD$nc)
  }
  
  #Setting parallel or sequential backend and  random seed
  if (ParExec & (parallel::detectCores() > 2 & mi>1)){
    if(is.null(ncores)){
      Ncpus <- parallel::detectCores() - 1
    }else{
      Ncpus <- min(ncores,parallel::detectCores() - 1)
    }
    cl <- parallel::makeCluster(Ncpus)
    doSNOW::registerDoSNOW(cl) #registerDoParallel doesn't have compatibility with ProgressBar
    if(SetRNGSeed){
      doRNG::registerDoRNG(SetRNGSeed)
    }
    # set progress bar for parallel processing
    pb <- txtProgressBar(max = mi, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # condition used to run code part needed for parallel processing
    ParParams = TRUE
  }else{ 
    if (ParExec & mi==1){
      if(verbose==T){
        message(paste("/!\\ The number of mi iteration is 1, parallel processing is only available for mi > 1."))
      }
    } else if (ParExec){
      if(verbose==T){
        message(paste("/!\\ The number of cores of your processor does not allow paralell processing, at least 3 cores are needed."))
      }
    }
    if(SetRNGSeed){
      set.seed(SetRNGSeed)
    }
    
    foreach::registerDoSEQ()
    opts = NULL
    
    # condition used to run code part needed for sequential processing
    ParParams  = FALSE 
  }
  
  
  #Beginning of the multiple imputation (imputing "mi" times)
  o <- NULL
  RESULT <- foreach(o=1:mi, .inorder = TRUE, .combine = "rbind", .options.snow = opts) %dopar% {
    if (!ParParams){
      # Parallel and sequential execution of foreach don't use the same casting mechanism, this one is used for sequential execution.
      if(verbose==T){
        cat("iteration :",o,"/",mi,"\n")
      }
    }
    # Trying to catch the potential singularity error (part 1/2)
    # (comment this part of code (as well as the one at the very end of
    # seqimpute.R) to debug more easily and see the full message of any
    # occuring error)
    #********************************************************************************************************************************
        #************************************************************************************************************************
        # 3. Imputation using a specific model --------------------------------------------------------------------------------------------------------
    if (max(dataOD$ORDER)!=0) {
      # Otherwise if there is only 0 in ORDER,
      # there is no need to impute internal gaps
      # and we directly jump to the imputation of
      # external gaps (i.e. points 4. and 5.)
      if(verbose==T){
        print("Imputation of the internal gaps...")
      }
      dataOD[["ODi"]]  <- ModelImputation(dataOD$OD, dataOD$CO, dataOD$COt, dataOD$ODi, dataOD$MaxGap, dataOD$totV, dataOD$totVi, regr, dataOD$nc, np, nf, dataOD$nr, dataOD$ncot, dataOD$COtsample, dataOD$pastDistrib, dataOD$futureDistrib, dataOD$k, available, dataOD$REFORD_L, dataOD$noise,num.trees,min.node.size,max.depth)
      
    }
    # 4. Imputing initial NAs -------------------------------------------------------------------------------------------------------------------------
    if ((nfi != 0) & (dataOD$MaxInitGapSize != 0)) {
      if(verbose==T){
        print("Imputation of the initial gaps...")
      }
      # # we only impute the initial gaps if nfi > 0
      dataOD[["ODi"]] <- ImputingInitialNAs(dataOD$CO, dataOD$COt, dataOD$OD, dataOD$ODi, dataOD$totVi, dataOD$COtsample, dataOD$futureDistrib, dataOD$InitGapSize, dataOD$MaxInitGapSize, dataOD$nr, dataOD$nc, dataOD$ud, dataOD$nco, dataOD$ncot, nfi, regr, dataOD$k, available, dataOD$noise,num.trees,min.node.size,max.depth)
    }
    # 5. Imputing terminal NAs ------------------------------------------------------------------------------------------------------------------------
    if ((npt != 0) & (dataOD$MaxTermGapSize != 0)){
      # we only impute the terminal
      # gaps if npt > 0
      if(verbose==T){
        print("Imputation of the terminal gaps...")
      }
      dataOD[["ODi"]]  <- ImputingTerminalNAs(dataOD$ODi, dataOD$CO, dataOD$OD, dataOD$COt, dataOD$COtsample, dataOD$MaxTermGapSize, dataOD$TermGapSize, dataOD$pastDistrib, regr, npt, dataOD$nco, dataOD$ncot, dataOD$totVt, dataOD$nr, dataOD$nc, dataOD$ud, available, dataOD$k, dataOD$noise,num.trees,min.node.size,max.depth)
    }
    #if (max(dataOD$ORDER)!=0) {
    # 6. Imputing SLG NAs -----------------------------------------------------------------------------------------------------------------------------
    # left-hand side SLG
    if (max(dataOD$ORDERSLGLeft)!=0) {
      # Checking if we have to impute
      # left-hand side SLG
      if(verbose==T){
        print("Imputation of the left-hand side SLG...")
      }
      dataOD[["ODi"]] <- LSLGNAsImpute(dataOD$OD, dataOD$ODi, dataOD$CO, dataOD$COt, dataOD$COtsample, dataOD$ORDERSLGLeft, dataOD$pastDistrib, dataOD$futureDistrib, regr, np, dataOD$nr, nf, dataOD$nc, dataOD$ud, dataOD$ncot, dataOD$nco, dataOD$k, dataOD$noise, available,num.trees,min.node.size,max.depth)
      
    }
    # right-hand side SLG
    if (max(dataOD$ORDERSLGRight)!=0) {
      # Checking if we have to impute right-hand
      # side SLG
      if(verbose==T){
        print("Imputation of the right-hand side SLG...")
      }
      
      dataOD[["ODi"]] <- RSLGNAsImpute(dataOD$OD, dataOD$ODi, dataOD$CO, dataOD$COt, dataOD$COtsample, dataOD$ORDERSLGRight, dataOD$pastDistrib, dataOD$futureDistrib, regr, np, dataOD$nr, nf, dataOD$nc, dataOD$ud, dataOD$ncot, dataOD$nco, dataOD$k, dataOD$noise, available,num.trees,min.node.size,max.depth)
      
    }
    # Checking if we have to impute
    # Both-hand side SLG
    if (dataOD$LongGap) {
      if(verbose==T){
        print("Imputation of the both-hand side SLG...")
      }
      for(h in 2:np){
        print(h)
        if(sum(dataOD$ORDERSLGBoth[,h-1]==0&dataOD$ORDERSLGBoth[,h]!=0)>0){
          tt <- which(dataOD$ORDERSLGBoth[,h-1]==0&dataOD$ORDERSLGBoth[,h]!=0)
          tmpORDER <- matrix(0,nrow(dataOD$ORDERSLGBoth),ncol(dataOD$ORDERSLGBoth))
          tmpORDER[tt,h:ncol(dataOD$ORDERSLGBoth)] <- dataOD$ORDERSLGBoth[tt,h:ncol(dataOD$ORDERSLGBoth)]
          dataOD[["ODi"]] <- RSLGNAsImpute(dataOD$OD, dataOD$ODi, dataOD$CO, dataOD$COt,dataOD$COtsample, tmpORDER, dataOD$pastDistrib, dataOD$futureDistrib, regr, h-1, dataOD$nr, nf, dataOD$nc, dataOD$ud, dataOD$ncot, dataOD$nco, dataOD$k, dataOD$noise, available,num.trees,min.node.size,max.depth)
          
        }
      }
      
    }
    #}
    
    #****************************************************************************************************************************************
    
    # Updating the matrix RESULT used to store the multiple imputations
    return(cbind(replicate(dataOD$nr,o), dataOD$ODi))
    
  }
  if (ParParams){
    parallel::stopCluster(cl)
  }
  RESULT <- rbind(cbind(replicate(dataOD$nr,0),dataOD$OD), RESULT)

  # X. Final conversions ----------------------------------------------------------------------------------------------------------------------------------------
  RESULT <- FinalResultConvert(RESULT, dataOD$ODClass, dataOD$ODlevels, rownamesDataset, nrowsDataset, dataOD$nr, dataOD$nc, dataOD$rowsNA, include, mi, mice.return)
  
  
  # Rearrange dataframe by order of the mi imputation as parallel computation may not return the values in sequential order.
  # if (ParParams){
  #   RESULT <- RESULT[order(RESULT$.imp),]
  # }
  
  return(RESULT)
}
