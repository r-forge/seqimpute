#' Spotting impossible transitions among a dataset
#'
#' The purpose of \code{seqImpTrans.R} is to spot the impossible transitions
#' in a dataset. These could have been introduced by mistake during data
#' collection or not intelligently imputed during a first run of
#' \code{seqimpute.R}.
#'
#' @param OD \code{matrix} object containing sequences of a variable with missing data (coded as \code{NA}).
#' @param k \code{numeric} object corresponding to the number of categories of the variable numbered from \code{1} to \code{k}.
#' @param impTrans \code{character} vector gathering the impossible transitions. For example: impTrans <- c("1->3","1->4","2->1","4->1","4->3")
#'
#' @author Andre Berchtold <andre.berchtold@@unil.ch>, Kevin Emery
#'
#' @return It returns a message containing the impTransOverview \code{data.frame} object that gathers the occurences of each type of impossible transitions.
#' @return rowMat \code{matrix} object containing the row coordinates of the impossible transitions (stored in \code{seqImpTransList[1]}).
#' @return colMat \code{matrix} object containing the column coordinates of the impossible transitions (stored in \code{seqImpTransList[2]}).
#'
#' @examples
#' data(OD)
#'
#' seqImpTransList <- seqImpTrans(OD=OD, k=2, impTrans=c("yes->no"))
#'
#' @keywords multinomial logistic regression, linear regression, ordinal regression, missing data, impossible transitions
#' 
#' @importFrom stringr str_count
#' @importFrom stringr str_detect
#' @importFrom stringr str_locate
#' @importFrom stringr str_locate_all
#' 
#' @importFrom graphics plot
#' 
#' @importFrom stats as.formula
#' @importFrom stats cutree
#' @importFrom stats lm
#' @importFrom stats predict
#' @importFrom stats rnorm
#' @importFrom stats runif
#' 
#' @importFrom utils capture.output
#' 
#' @importFrom swfscMisc na.count
#' 
#' @importFrom Amelia missmap
#' 
#' @importFrom TraMineR seqdef
#' @importFrom TraMineR seqfplot
#' @importFrom TraMineR seqdplot
#' @importFrom TraMineR seqsubm
#' @importFrom TraMineR seqdist
#' 
#' @importFrom cluster agnes
#' 
#' @importFrom plyr mapvalues
#' 
#' @importFrom dfidx dfidx
#' @importFrom mlogit mlogit
#' 
#' @importFrom rms lrm
#'
#' @export


seqImpTrans <- function(OD, k, impTrans){
  
  
  # Naming the number of rows and columns of OD
  nr <- nrow(OD)
  nc <- ncol(OD)
  
  
  # test
  
  # 1. Initial tests on parameters --------------------------------------------------------------------------------------------------------------------
  # 1.1 Testing the class of the variables of the original dataset OD ---------------------------------------------------------------------------------
  ODClass <- class(OD[1,1])
  if ( (ODClass != "factor") & (ODClass != "numeric") ) {
    stop("/!\\ The class of the variables contained in your original dataset
         should be either 'factor' or 'numeric'")
  }
  
  # 1.2 Testing effectively exactly k possible categories of the multinomial variable -----------------------------------------------------------------
  if (class(OD[1,1]) == "numeric") {
    for (i in 1:nr) {
      for (j in 1:nc) {
        if ( is.na(OD[i,j])==FALSE & (OD[i,j]<=0 | OD[i,j]>k) ) {
          stop("/!\\ Your dataset doesn't contain the right number of k categories of the multinomial variable")
        } else {
          next
        }
      }
    }
  } else { # Meaning that our values are of type "factor"
    for (i in 1:nr) {
      for (j in 1:nc) {
        if (length(levels(OD[i,j])) != k) {
          stop("/!\\ Your dataset doesn't contain the right number of k categories of the multinomial variable")
        } else {
          next
        }
      }
    }
  }
  
  
  # 1.3 Indicating the existence of entire rows of OD filled only with NAs to the user ----------------------------------------------------------------
  i <- 1
  numbOfNAFilledLines <- 0
  while (i <= nrow(OD)) {
    if (all(is.na(OD[i,]))) {
      OD <- OD[-i,]
      numbOfNAFilledLines <- numbOfNAFilledLines + 1
      warning(paste("/!\\ Row number",i,"of OD consists only of NAs."),sep='')
    }
    i <- i+1
  }
  if (numbOfNAFilledLines == 1) {
    warning(paste("Your data matrix contains 1 row filled solely with NAs. This row won't be included into the next statistics about NAs."),sep='')
  }
  if (numbOfNAFilledLines > 1) {
    warning(paste("Your data matrix contains",numbOfNAFilledLines,"rows filled solely with NAs."),sep='',"\n",
            "  These rows have been removed and thus haven't been included in","\n",
            "  the analysis of the impossible transitions among your dataset.")
  }
  # Updating the number of rows in OD
  nr <- nrow(OD)
  
  
  # 1.4 Test on input parameter impTrans -----------------------------------------
  for (i in 1:length(impTrans)) {
    if (!str_detect(impTrans[i],"->")) {
      stop("/!\\ Warning, you should construct your impossible transition(s) vector impTrans with little arrows as follow: impTrans <- c('...->...', '...->...', etc.).")
    }
    # Testing if we are effectively analyzing a transition or not
    locDash <- str_locate(impTrans[i],"-")
    firstState <- substr(impTrans[i],1,locDash-1)
    locSpike <- str_locate(impTrans[i],">")
    secondState <- substr(impTrans[i],locSpike+1,nchar(impTrans[1]))
    # if (firstState==secondState) {
    # stop("/!\\ You have typed in two same states (two times '",firstState,"') on both sides of the arrow. This doesn't correspond to a transition.")
    # }
  }
  
  
  
  # 2. Finding number of gaps -------------------------------------------------------------------------------------------------------------------------
  
  # 2.1 Number of initial gaps ------------------------------------------------------------------------------------------------------------------------
  #
  # Creation of matrix ORDER
  ORDER <- matrix(0,nr,nc) # initialization of matrix ORDER with 0 everywhere
  SEL <- is.na(OD)==TRUE   # creation of matrix SEL, constituted of TRUE where there is MD in OD and of FALSE everywhere else
  ORDER[SEL] <- 1          # setting some 1 in ORDER at the location where in SEL we have some TRUE
  
  # Creation of vector InitGapSize (i.e. a vector containing the size of the initial gaps of each line)
  InitGapSize <- vector()
  for (i in 1:nr) {
    if (ORDER[i,1]==0) {
      InitGapSize[i] <- 0
    } else {
      InitGapSize[i] <- 1
      for (j in 2:nc) {
        if (ORDER[i,j]==1) {
          InitGapSize[i] <- InitGapSize[i] + 1
        } else {
          break
        }
      }
    }
  }
  MaxInitGapSize <- max(InitGapSize)
  # Number of initial gaps
  numbOfInitGaps <- length(InitGapSize[InitGapSize!=0])
  
  # Creation of matrix ORDERI
  ORDERI <- matrix(0,nr,nc)
  for (i in 1:nr) {
    if (InitGapSize[i]!=0) {
      ORDERI[i,1:InitGapSize[i]] <- c(MaxInitGapSize:(MaxInitGapSize+1-InitGapSize[i]))
    } else {
      next
    }
  }
  # Replacing each value of ORDERI greater than '0' by '1'
  ORDERI[ORDERI > 0] <- 1
  
  
  
  
  
  # 2.2 Number of terminal gaps -----------------------------------------------------------------------------------------------------------------------
  # Creation of vector TermGapSize (i.e. a vector containing the size of the terminal gaps of each line)
  TermGapSize <- vector()
  for (i in 1:nr) {
    if (ORDER[i,nc]==0) {
      TermGapSize[i] <- 0
    } else {
      TermGapSize[i] <- 1
      for (j in (nc-1):1) {
        if (ORDER[i,j]==1) {
          TermGapSize[i] <- TermGapSize[i] + 1
        } else {
          break
        }
      }
    }
  }
  MaxTermGapSize <- max(TermGapSize)
  # Number of terminal gaps
  numbOfTermGaps <- length(TermGapSize[TermGapSize!=0])
  
  # Creation of matrix ORDERT
  ORDERT <- matrix(0,nr,nc)
  for (i in 1:nr) {
    if (TermGapSize[i]!=0) {
      ORDERT[i,(nc-TermGapSize[i]+1):nc] <- c((MaxTermGapSize+1-TermGapSize[i]):MaxTermGapSize)
    } else {
      next
    }
  }
  # Replacing each value of ORDERT greater than '0' by '1'
  ORDERT[ORDERT > 0] <- 1
  
  
  
  
  
  # 2.3 Number of internal gaps -----------------------------------------------------------------------------------------------------------------------
  
  # /!\ Final version of the matrix ORDER that we use through point 3.1 to 3.3 of the program
  ORDER <- ORDER - ORDERI - ORDERT
  
  # Number of Internal Gaps
  # Transforming ORDER in a single row vector
  ORDERVect <- as.vector(t(ORDER))
  # Transforming this single row vector into class character
  ORDERVectChar <- paste(ORDERVect,collapse="")
  # Identifying the patterns "0 1" (this is the signature look of an internal gap
  # (it always indicates the beginning of an internal gap!))
  numbOfInternGaps <- str_count(ORDERVectChar,pattern="01")
  
  
  
  
  # Warning message returning the number of gaps of each type included in OD
  if ( (numbOfInitGaps > 0) | (numbOfTermGaps > 0) | (numbOfInternGaps > 0) ) {
    warning("/!\\ We have detected ",numbOfInitGaps," initial gap(s), ",numbOfTermGaps," terminal gap(s)","\n",
            "    and ",numbOfInternGaps," internal gap(s) in your dataset.","\n",
            "    Be aware that these gaps may hide other impossible transitions!")
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  # 3. Spotting impossible transitions ----------------------------------------------------------------------------------------------------------------
  
  
  
  ## Setup
  #
  # Transforming every line of OD into class character and adding a dash inbetween
  # every state
  dashes <- replicate(nc,"->")
  CharAndDashes <- function(OD) {
    mytestChar <- paste(as.vector(rbind(OD,dashes)),collapse="")
  }
  ODCharAndDashes <- apply(OD,1,CharAndDashes)
  
  
  
  
  
  
  # 3.1 Number of listed impossible transitions -------------------------------------------------------------------------------------------------------
  # Identifying the patterns of the impossible transitions in each line of OD
  countImpTrans <- function(ODCharAndDashes) {
    str_count(ODCharAndDashes,pattern=paste0("(?=",impTrans,")"))
  }
  numbOfImpTransByRow <- sapply(ODCharAndDashes,countImpTrans)
  if(length(impTrans)>1) {
    numbOfImpTrans <- rowSums(numbOfImpTransByRow)
  } else {
    numbOfImpTrans <- sum(numbOfImpTransByRow)
  }
  # Interrupting the program in case no impossible transitions among impTrans have
  # been found
  if (sum(numbOfImpTrans) == 0) {
    stop(" /!\\ Warning, no impossible transitions have been found. Your input vector impTrans doesn't contain any transitions present in your dataset.")
  }
  
  if(length(impTrans)==1){
    numbOfImpTransByRow <- t(as.matrix(numbOfImpTransByRow))
  }
  
  
  
  
  
  startLocMat <- matrix(NA,nrow=length(impTrans),ncol=nrow(OD))
  for(i in 1:length(impTrans)){
    Tmplist <- str_locate_all(ODCharAndDashes, paste0("(?=",impTrans[i],")"))
    for(j in 1:length(ODCharAndDashes)){
      tempMat <- Tmplist[[j]]
      if (nrow(tempMat) > 0) {
        if (nrow(tempMat) > 1) {
          tempMat <- paste(tempMat[,1],collapse = " ")
        }
        startLocMat[i,j] <- as.matrix(tempMat[1])
      }
    }
  }
  
  
  
  
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }
  GetRealColPosition <- function(ODCharAndDashes,RowPos,ColPos) {
    ODCharAndDashesSubstring <- substr(ODCharAndDashes[RowPos],1,ColPos)
    realColPosition <- str_count(ODCharAndDashesSubstring,pattern=">")+1
    return(realColPosition)
  }
  list_rows <- list()
  list_cols <- list()

  for(i in 1:length(impTrans)){
    if(numbOfImpTrans[i]>0){
      TmpRows<-c()
      TmpCols<-c()
      for(j in 1:ncol(startLocMat)){
        if(numbOfImpTransByRow[i,j]>0){
          detectTemp <- str_detect(startLocMat[i,j]," ")
          if(detectTemp==FALSE){
            TmpRows <- c(TmpRows,j)
            TmpCols <- c(TmpCols, GetRealColPosition(ODCharAndDashes = ODCharAndDashes, RowPos = j, ColPos =  as.numeric(startLocMat[i,j])))
          }else{
            hiddenCol <- str_count(startLocMat[i,j]," ")
            TmpRows <- c(TmpRows,rep(j,hiddenCol+1))
            extractedNumbers <- Numextract(startLocMat[i,j])
            TmpVec <- c()
            for(k in 1:length(extractedNumbers)){
              TmpVec[k] <- GetRealColPosition(ODCharAndDashes = ODCharAndDashes, RowPos = j, ColPos =  as.numeric(extractedNumbers[k]))
            }
            TmpCols <- c(TmpCols,TmpVec)
          }
        }
      }
      list_rows[[i]]<-TmpRows
      list_cols[[i]]<-TmpCols
    }
  }
  
  
  MaxNumTransitions <- max(sapply(list_cols,function(x) length(x)))

  rowMat <- matrix(NA,length(impTrans),MaxNumTransitions)
  colMat <- matrix(NA,length(impTrans),MaxNumTransitions)
  
  # Converting into a dataframe
  rowMat <- as.data.frame(rowMat)
  # Renaming the columns of rowMat
  colnames_rowMat <- paste(1:ncol(rowMat),")",sep='')
  rownames(rowMat) <- impTrans
  colnames(rowMat) <- colnames_rowMat
  
  # Converting into a dataframe
  colMat <- as.data.frame(colMat)
  # Renaming the columns of colMat
  colnames_colMat <- paste(1:ncol(colMat),")",sep='')
  rownames(colMat) <- impTrans
  colnames(colMat) <- colnames_colMat
  
  for(i in 1:length(impTrans)){
    if(numbOfImpTrans[i]>0){
      rowMat[i,1:length(list_rows[[i]])] <- list_rows[[i]]
      colMat[i,1:length(list_cols[[i]])] <- list_cols[[i]]
    }
  }
  
  if(length(impTrans)>1){
    rowMat <- rowMat[rowSums(!is.na(rowMat))>0,]
    colMat <- colMat[rowSums(!is.na(colMat))>0,]
  }
  
  
  
  # 3.3 Summary data frame ----------------------------------------------------------------------------------------------------------------------------
  #
  # Update of impTrans with an arrow
  impTransOverview <- data.frame(c(impTrans,"","Total:"),c(numbOfImpTrans,"",sum(numbOfImpTrans)))
  colnames(impTransOverview) <- c("Impossible transitions", "Occurence")
  
  
  message("Impossible transitions summarizing table:","\n",
          "-----------------------------------------","\n",
          "Note:","\n",
          "The location of these spotted impossible transitions can","\n",
          "be found by looking in parallel at the matrices 'rowMat' and","\n",
          "'colMat': the row in which the corresponding in OD row (resp. the column)","\n",
          "index appears indicates the case of impossible transition you are","\n",
          "looking at and the column informs about the rank of this specific","\n",
          "impossible transition (i.e. if it is the first (1)) impossible transition","\n",
          "of this kind that is met or the second (2)), etc.)","\n","\n",
          
          paste0(capture.output(impTransOverview), collapse = "\n"),"\n")
  
  
  
  
  
  
  
  
  
  seqImpTransList <- list("rowMat" = rowMat, "colMat" = colMat)
  return(seqImpTransList)
  
  
  
  }
