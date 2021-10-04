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
#' @return It returns a message containing the MDGapsChart \code{data.frame} object that summarizes the number and the
#' size of each type of gaps (Internal Gaps, Initial Gaps, Terminal Gaps, LEFT-
#' hand side SLG, RIGHT-hand side SLG) present among the original dataset.
#'
#' @examples
#' data(OD)
#'
#' seqQuickLook(OD=OD, np=1, nf=0)
#'
#' @export


seqQuickLook <- function(OD, np=1, nf=0) {


    k <- n_distinct(data.frame(newcol = c(t(OD)), stringsAsFactors=FALSE),na.rm = TRUE)



    # Naming the number of rows and columns of the original dataset OD
    nr <- nrow(OD)
    nc <- ncol(OD)









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
                "  These rows haven't been included in the missing data overview of OD.",
                "\n")
    }
    # Updating the number of rows in OD
    nr <- nrow(OD)




















    # 2. Indexing NAs among the different types of gaps -------------------------------------------------------------------------------------------------

    # 2.1 NAs among Initial Gaps ------------------------------------------------------------------------------------------------------------------------
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
    MinInitGapSize <- min(InitGapSize[InitGapSize!=0])
    MaxInitGapSize <- max(InitGapSize)
    # Number of initial gaps
    numbOfInitGaps <- length(InitGapSize[InitGapSize!=0])
    sumNAInitGaps <- sum(InitGapSize)










    # 2.2 NAs among Terminal Gaps -----------------------------------------------------------------------------------------------------------------------
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
    MinTermGapSize <- min(TermGapSize[TermGapSize!=0])
    MaxTermGapSize <- max(TermGapSize)
    # Number of terminal gaps
    numbOfTermGaps <- length(TermGapSize[TermGapSize!=0])
    sumNATermGaps <- sum(TermGapSize)












    # 2.3 NAs among SLG (Specially Located Gaps) --------------------------------------------------------------------------------------------------------

    # Updating of ORDER with "0" on every external NAs
    # (The purpose of this modification of ORDER is that we don't take into account external NAs at this moment of the program.
    # We will first impute internal gaps and consider external gaps further (as far as nfi and npt are greater than 0))
    for (i in 1:nr) {
        if (InitGapSize[i]!=0) {
            ORDER[i,1:InitGapSize[i]] <- vector("numeric",InitGapSize[i])
        }

        if (TermGapSize[i]!=0) {
            ORDER[i,(nc-TermGapSize[i]+1):nc] <- vector("numeric",TermGapSize[i])
        }

        else {
            next
        }
    }



    # Creation of matrices ORDER2 and ORDER3
    ORDER2 <- ORDER                                     # initially both ORDER2 and
    ORDER3 <- ORDER                                     # ORDER3 are equal to ORDER

    for (i in 1:nr){                                      # in matrix ORDER, we go line by line...
        for (j in 2:nc){                                    # ... from column to column (actually exactly as we read a book) (and /!\ beginning from column 2!)
            if (ORDER[i,j-1]==1 & ORDER[i,j]==1){             # if the previous value of ORDER is equal to 1 and the its actual value is also equal to 1, then...
                ORDER2[i,j] <- ORDER2[i,j-1]+1                  # ... construction of ORDER2 for this iteration: the current (j) corresponding value of ORDER2 is assigned by the its previous (j-1) value incremented by 1
                ORDER3[i,(j-(ORDER2[i,j]-1)):j] <- ORDER2[i,j]  # ... construction of ORDER3 for this iteration: all the values in ORDER3 from the beginning of the gap (j-(ORDER2[i,j]-1)) up to the current (j) location are assigned by the current (j) corresponding value in ORDER2
            }
        }
    }


    # Updating ORDER with "0" on every NAs belonging to a Specially Located Gap (SLG)
    # (The purpose of this modification of ORDER is that we don't take into account SLG NAs at this moment of the program.
    # We will first impute internal gaps, external gaps and consider SLG at the very end
    # (as far as some SLG have been detected)

    # Creation of ORDERSLG (ORDERSLGLeft and ORDERSLGRight)

    # Initialization of matrix in which we will store the SLG
    ORDERSLG <- matrix(0,nrow=nr,ncol=nc)

    # Initialization of the range in which SLG could be found
    tempMinGapLeft <- matrix(0,nrow=nr,ncol=nc)
    tempMaxGapLeft <- matrix(0,nrow=nr,ncol=nc)
    tempMinGapRight <- matrix(0,nrow=nr,ncol=nc)
    tempMaxGapRight <- matrix(0,nrow=nr,ncol=nc)

    for (i in 1:nr) {                           # we will go through each line of ORDER

        if (np > 1) {                           # if np > 1, it may be possible that SLG on the left-hand side of OD exist
            j <- 2

            while (j <= np) {
                jump <- 1

                if (ORDER[i,j]>0) {
                    tempMinGapLeft[i,j] <- j

                    while (ORDER[i,j]>0) {
                        ORDERSLG[i,j] <- ORDER[i,j]
                        j <- j+1
                    }

                    tempMaxGapLeft[i,j] <- j-1

                    jump <- max(tempMaxGapLeft[i,]) - max(tempMinGapLeft[i,])
                }

                j <- j+jump
            }
        }

        if (nf > 1) {                           # if nf > 1, it may be possible that SLG on the right-hand side of OD exist
            j <- nc-1

            while ((nc-j+1) <= nf) {
                jump <- 1

                if (ORDER[i,j]>0) {
                    tempMinGapRight[i,j] <- j

                    while (ORDER[i,j]>0) {
                        ORDERSLG[i,j] <- ORDER[i,j]
                        j <- j-1
                    }

                    tempMaxGapRight[i,j] <- j+1

                    jump <- max(tempMinGapRight[i,]) - max(tempMaxGapRight[i,])
                }

                j <- j-jump
            }
        }

    }

    # Extracting extrema from tempMinGapLeft, tempMaxGapLeft, tempMinGapRight and tempMaxGapRight
    # And creation of ORDERSLGLeft and ORDERSLGRight (matrices for both groups of SLG (i.e. one on
    # the left- and the other one on the right-hand side of the matrix ORDERSLG)
    ORDERSLGLeft <- matrix(nrow=nr,ncol=nc,0)
    if (max(tempMinGapLeft)!=0 & max(tempMaxGapLeft)!=0) {
        minGapLeft <- min(tempMinGapLeft[tempMinGapLeft!=0])
        maxGapLeft <- max(tempMaxGapLeft[tempMaxGapLeft!=0])
        ORDERSLGLeft[,minGapLeft:maxGapLeft] <- ORDERSLG[,minGapLeft:maxGapLeft]
    }

    ORDERSLGRight <- matrix(nrow=nr,ncol=nc,0)
    if (max(tempMinGapRight!=0) & max(tempMaxGapRight)!=0) {
        minGapRight <- max(tempMinGapRight[tempMinGapRight!=0])
        maxGapRight <- min(tempMaxGapRight[tempMaxGapRight!=0])
        ORDERSLGRight[,maxGapRight:minGapRight] <- ORDERSLG[,maxGapRight:minGapRight]
    }





    # 2.3.1 NAs among LEFT-hand side SLG ----------------------------------------------------------------------------------------------------------------

    # First LEFT-hand side SLG location
    MaxElColORDERSLGLeft <- apply(ORDERSLGLeft,2,max)
    for (j in 1:nc) {
        if (MaxElColORDERSLGLeft[j]>0) {
            LeftSLGBeginningLocation <- j
            break
        }
    }

    # Creation of matrices ORDERSLGLeft2 and ORDERSLGLeft3
    ORDERSLGLeft2 <- ORDERSLGLeft
    ORDERSLGLeft3 <- ORDERSLGLeft
    for (i in 1:nr){
        for (j in 2:nc){
            if (ORDERSLGLeft[i,j-1]==1 & ORDERSLGLeft[i,j]==1){
                ORDERSLGLeft2[i,j] <- ORDERSLGLeft2[i,j-1]+1
                ORDERSLGLeft3[i,(j-(ORDERSLGLeft2[i,j]-1)):j] <- ORDERSLGLeft2[i,j]
            }
        }
    }


    if (max(ORDERSLGLeft3) == 0) {
        MinSLGLeftGapSize <- 0
    } else {
        MinSLGLeftGapSize <- min(ORDERSLGLeft3[ORDERSLGLeft3!=0])
    }
    MaxSLGLeftGapSize <- max(max(ORDERSLGLeft2))

    # Number of RIGHT-hand side SLG Gaps
    # Transforming ORDERSLGRight in a single row vector
    ORDERSLGLeftVect <- as.vector(t(ORDERSLGLeft))
    # Transforming this single row vector into class character
    ORDERSLGLeftVectChar <- paste(ORDERSLGLeftVect,collapse="")
    # Identifying the patterns "0 1" (this is the signature look of an internal gap
    # (it always indicates the beginning of an internal gap!))
    numbOfSLGLeftGaps <- str_count(ORDERSLGLeftVectChar,pattern="01")

    sumNASLGLeftGaps <- sum(ORDERSLGLeft)








    # 2.3.2 NAs among RIGHT-hand side SLG ---------------------------------------------------------------------------------------------------------------

    # First RIGHT-hand side SLG location
    MaxElColORDERSLGRight <- apply(ORDERSLGRight,2,max)
    for (j in nc:1) {
        if (MaxElColORDERSLGRight[j]>0) {
            RightSLGBeginningLocation <- nc-j+1
            break
        }
    }

    # Creation of matrices ORDERSLGRight2 and ORDERSLGRight3
    ORDERSLGRight2 <- ORDERSLGRight
    ORDERSLGRight3 <- ORDERSLGRight
    for (i in 1:nr){
        for (j in 2:nc){
            if (ORDERSLGRight[i,j-1]==1 & ORDERSLGRight[i,j]==1){
                ORDERSLGRight2[i,j] <- ORDERSLGRight2[i,j-1]+1
                ORDERSLGRight3[i,(j-(ORDERSLGRight2[i,j]-1)):j] <- ORDERSLGRight2[i,j]
            }
        }
    }


    if (max(ORDERSLGRight3) == 0) {
        MinSLGRightGapSize <- 0
    } else {
        MinSLGRightGapSize <- min(ORDERSLGRight3[ORDERSLGRight3!=0])
    }
    MaxSLGRightGapSize <- max(max(ORDERSLGRight2))

    # Number of right-hand side SLG Gaps
    # Transforming ORDERSLGRight in a single row vector
    ORDERSLGRightVect <- as.vector(t(ORDERSLGRight))
    # Transforming this single row vector into class character
    ORDERSLGRightVectChar <- paste(ORDERSLGRightVect,collapse="")
    # Identifying the patterns "0 1" (this is the signature look of an internal gap
    # (it always indicates the beginning of an internal gap!))
    numbOfSLGRightGaps <- str_count(ORDERSLGRightVectChar,pattern="01")

    sumNASLGRightGaps <- sum(ORDERSLGRight)












    # 2.4 NAs among Internal Gaps -----------------------------------------------------------------------------------------------------------------------

    # /!\ Final version of the matrix ORDER that we use through point 3.1 to 3.3 of the program
    ORDER <- ORDER - ORDERSLGLeft - ORDERSLGRight


    MinInternGapSize <- min(ORDER3[ORDER3!=0])
    MaxInternGapSize <- max(max(ORDER2)) # renders us the size of the greatest gap in OD

    # Number of Internal Gaps
    # Transforming ORDER in a single row vector
    ORDERVect <- as.vector(t(ORDER))
    # Transforming this single row vector into class character
    ORDERVectChar <- paste(ORDERVect,collapse="")
    # Identifying the patterns "0 1" (this is the signature look of an internal gap
    # (it always indicates the beginning of an internal gap!))
    numbOfInternGaps <- str_count(ORDERVectChar,pattern="01")

    sumNAInternGaps <- sum(ORDER)











    # 3. Overall NA indexing ----------------------------------------------------------------------------------------------------------------------------

    # Replacing OD as a single vector
    # We have then like one single column of nrxnc rows
    # (every initial column of OD are then stacked on each other)
    # i.e. every nr+1 (i.e. 1, nr+1, (2*nr)+1, (3*nr)+1, etc.) element we are situated like
    # at the top beginning of a new column of OD
    ODvector <- as.vector(OD)

    NaReplacedBy1s_col <- ifelse(is.na(ODvector),1,0)







    # 3.1 Total number of NAs and occurrence percentage -------------------------------------------------------------------------------------------------
    totSumNA <- sum(NaReplacedBy1s_col)

    percNA <- paste(round(totSumNA/(nr*nc)*100,3),"%")

    totNumbOfGaps <- numbOfInternGaps + numbOfInitGaps + numbOfTermGaps + numbOfSLGLeftGaps + numbOfSLGRightGaps








    # 3.2 Number of NAs by column and occurrence percentage ---------------------------------------------------------------------------------------------
    # sumNA vectors (one for every column of OD)
    shift <- 1
    for (i in 1:nc) {
        assign(paste("sumNA_col_",i,sep=''),sum(NaReplacedBy1s_col[shift:(i*nr)]))
        tempObject = get(paste0("sumNA_col_",i,sep=''))
        assign(paste("percNA_col_",i,sep=''),round(tempObject/nr,3))
        shift <- shift+nr
    }









    # 3.3 Number of NAs by row and occurence percentage -------------------------------------------------------------------------------------------------
    # sumNA vectors (one for every row of OD)
    ODtvector <- as.vector(t(OD)) # transposing OD and taking its corresponding vector composed of its rows
    NaReplacedBy1s_row <- ifelse(is.na(ODtvector),1,0)
    # Every nc+1 (i.e. 1, nc+1, (2*nc)+1, (3*nc)+1, etc.) element we are situated like
    # at the beginning of a new row of OD
    shift <- 1
    for (i in 1:nr) {
        assign(paste("sumNA_row_",i,sep=''),sum(NaReplacedBy1s_row[shift:(i*nc)]))
        tempObject = get(paste0("sumNA_row_",i,sep=''))
        assign(paste("percNA_row_",i,sep=''),round(tempObject/nc,3))
        shift <- shift+nc
    }








    # 3.4 Number of sequences (rows) wiht and without NAs
    numbOfSeqWithNA <- na.count(apply(OD,1,max))
    numbOfSeqWithoutNA <- nr - numbOfSeqWithNA

























    # 4. Summarizing NAs analysis results and giving tips on parameters setting -------------------------------------------------------------------------

    # 4.1 Overall Missing data information --------------------------------------------------------------------------------------------------------------
    #print(paste("Total number of missing data among OD: ",totSumNA))
    print(paste("Missing data occurrence percentage:     ",percNA))
    #print(paste("Total number of gaps among OD:         ",totNumbOfGaps))

    # 4.2 Missing data among gaps summarizing chart -----------------------------------------------------------------------------------------------------
    TypeOfGaps <- c("Internal Gaps","Initial Gaps","Terminal Gaps","LEFT-hand side SLG","RIGHT-hand side SLG","Total")
    MinGapSize <- c(MinInternGapSize,MinInitGapSize,MinTermGapSize,MinSLGLeftGapSize,MinSLGRightGapSize," ")
    MaxGapSize <- c(MaxInternGapSize,MaxInitGapSize,MaxTermGapSize,MaxSLGLeftGapSize,MaxSLGRightGapSize," ")
    numbOfGaps <- c(numbOfInternGaps,numbOfInitGaps,numbOfTermGaps,numbOfSLGLeftGaps,numbOfSLGRightGaps,sum(numbOfInternGaps,numbOfInitGaps,numbOfTermGaps,numbOfSLGLeftGaps,numbOfSLGRightGaps))
    sumNAGaps <- c(sumNAInternGaps,sumNAInitGaps,sumNATermGaps,sumNASLGLeftGaps,sumNASLGRightGaps,sum(sumNAInternGaps,sumNAInitGaps,sumNATermGaps,sumNASLGLeftGaps,sumNASLGRightGaps))


    MDGapsChart <- data.frame(TypeOfGaps,MinGapSize,MaxGapSize,numbOfGaps,sumNAGaps)
    #print(MDGapsChart)




    # 4.3 Warnings concerning SLG -----------------------------------------------------------------------------------------------------------------------
    if (max(ORDERSLGLeft) > 0) {
        warning("Suggestion on np:","\n",
                "-----------------------------------------","\n",
                "Since the first SLG on the left begins in column ",LeftSLGBeginningLocation,",","\n",
                "it would be recommanded to choose np strictly inferior to this value...","\n",
                "\n")
    }

    if (max(ORDERSLGRight) > 0) {
        warning("Suggestion on nf:","\n",
                "-----------------------------------------","\n",
                "Since the first SLG on the right-hand side begins in column ",RightSLGBeginningLocation," (regarding from the right),","\n",
                "it would be recommanded to choose nf strictly inferior to this value...")
    }



















    # 4.4 Displaying interesting results to the user ----------------------------------------------------------------------------------------------------

    message(paste("Your OD is made up of ",numbOfSeqWithNA," sequences containing at least one NA,","\n",
                  "and ",numbOfSeqWithoutNA," sequences containing absolutely no NA.", sep=''),"\n","\n",
            "Summarizing Gap Chart:","\n",
            "----------------------","\n",
            paste0(capture.output(MDGapsChart), collapse = "\n"),"\n")









}
