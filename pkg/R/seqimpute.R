#' Imputation of missing data in sequence analysis
#'
#' Imputation of missing data present in a dataset through the prediction based
#' on either a multinomial, a linear or an ordinal regression model.
#' In order to specify even more the prediction, fixed as well as time-dependant
#' covariates be included in the model.
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
#' - Left-hand side SLG (Specially Located Gaps) (gaps of which the beginning location is included in the interval \code{[0,np]})
#'
#' - Right-hand side SLG (Specially Located Gaps) (gaps of which the ending location is included in the interval \code{[ncol(OD)-nf,ncol(OD)]})
#'
#' Order of imputation of the gaps types:
#'     1. Internal Gaps
#'     2. Initial Gaps
#'     3. Terminal Gaps
#'     4. Left-hand side SLG
#'     5. Right-hand side SLG
#'
#' @param OD \code{matrix} object containing sequences of a multinomial variable with missing data (coded as \code{NA}).
#' @param regr \code{character} object corresponding to the type of regression model the user want to use to compute. The prediction (either multinomial with "\code{mlogit}", linear with "\code{lm}" or ordinal with "\code{lrm}") (default \code{mlogit}).
#' @param k \code{numeric} object corresponding to the number of categories of the multinomial variable numbered from \code{1} to \code{k}.
#' @param np \code{numeric} object corresponding to the number of previous observations in the imputation model of the internal gaps (default \code{1}).
#' @param nf \code{numeric} object corresponding to the number of future observations in the imputation model of the internal gaps (default \code{0}).
#' @param nfi \code{numeric} object corresponding to the number of future observations in the imputation model of the initial gaps (default \code{1}).
#' @param npt \code{numeric} object corresponding to the number of previous observations in the imputation model of the terminal gaps (default \code{1}).
#' @param available \code{logical} object allowing the user to choose whether to consider the already imputed data in our predictive model (\code{available = TRUE}) or not (\code{available = FALSE}) (default \code{TRUE}).
#' @param CO \code{data.frame} object containing some covariates among which the user can choose in order to specify his model more accurately (default empty matrix 1x1 (\code{matrix(NA,nrow=1,ncol=1)})).
#' @param COt \code{data.frame} object containing some time-dependent covariates that help specifying the predictive model more accurately (default empty matrix 1x1 (\code{matrix(NA,nrow=1,ncol=1)})).
#' @param pastDistrib \code{logical} object allowing to take account of the past distribution in the multinomial logistic regression model or not (default \code{FALSE}).
#' @param futureDistrib \code{logical} object allowing to take account of the future distribution in the multinomial logistic regression model or not (default \code{FALSE}).
#' @param mi \code{numeric} object corresponding to the number of imputations the program is going to perform (default: \code{1}).
#' @param mi.return \code{numeric} object being either "\code{1}" (omitting \code{OD}) or "\code{2}" (including \code{OD}) and corresponding to the two possible returned formats of the final matrix \code{RESULT} (default: \code{1}).
#' @param noise \code{numeric} object adding a noise on the predicted variable \code{pred} determined by the multinomial model (by introducing a variance \code{noise} for each components of the vector \code{pred}) (the user can choose any value for \code{noise}, but we recommend to choose a rather relatively small value situated in the interval \code{[0.005-0.03]}) (default \code{0}).
#'
#' @author Andre Berchtold <andre.berchtold@@unil.ch>
#'
#' @return RESULT \code{data.frame} object containing the imputed original dataset.
#'
#' @examples
#' library(seqimpute)
#' data(OD, CO, COt, package="seqimpute")
#'
#' RESULT <- seqimpute(OD=OD, k=2, np=1, nf=0, nfi=1, npt=1, available=TRUE, CO=CO, COt=COt, pastDistrib=FALSE, futureDistrib=FALSE, mi=1, mi.return=1, noise=0)
#'
#' @references HALPIN, Brendan, March 2013. Imputing Sequence Data : Extensions to initial and terminal gaps, Stata's mi. Unviversity of Limerick Department of Sociology Working Paper Series. Working Paper WP2013-01, p.3. Available at : http://www.ul.ie/sociology/pubs/wp2013-01.pdf
#'
#' @keywords multinomial logistic regression, missing data, missing patterns
#'
#' @import mlogit
#' @import rms
#' @import stringr
#' @import Amelia
#' @import mlbench
#' @import TramineR 
#' @import cluster
#' @import swfscMisc 
#' @import plyr
#' @export


seqimpute <- function(OD, regr="mlogit", k, np=1, nf=0, nfi=1, npt=1,
                      available=TRUE, CO=matrix(NA,nrow=1,ncol=1),
                      COt=matrix(NA,nrow=1,ncol=1), pastDistrib=FALSE,
                      futureDistrib=FALSE, mi=1, mi.return=1, noise=0) {



# test




    # Selecting the columns of CO the user finally wants to use in his model
    #*******************************************************************************
    # if (ncol(CO) > 1) {          # if nco is greater than "1", it means that
    #     # the covariate matrix CO is composed of several columns,
    #     # that is, the user has to choose which column of CO (i.e. which
    #     # covariates) he wants to take in his modelS
    #     colNamesCO <- colnames(CO)
    #     takenCO <- matrix(nrow=nrow(CO),ncol=0)
    #
    #     message(paste("We have identified",ncol(CO),"covariates in your
    #                       covariate matrix CO."))
    #     message("Please, select your covariates among your covariate
    #                 matrix CO:")
    #     for (i in 1:ncol(CO)) {
    #         covChoice <- readline(prompt=paste(i,". Type in [y] followed by
    #                 [Enter] if you desire to consider the covariate ",
    #                                            colNamesCO[i]," in your model: ","\n","(otherwise simply press
    #                 [Enter] or tap anything and then press [Enter]) ",sep=''))
    #         if (covChoice == 'y') {
    #             takenCO <- cbind(takenCO,CO[i])
    #         } else {
    #             next
    #         }
    #     }
    #
    #     message("Here is a preview of your selected covariate(s):")
    #     if (all(is.na(CO))==FALSE) { # if CO is NOT completely empty
    #         print(head(CO))
    #     } else {
    #         print("NA (no covariates selected)")
    #     }
    #     invisible(readline(prompt="Press [Enter] to continue or relaunch the
    #                            program to change this selection of covariates..."))
    # }
    # # Otherwise, if the user sets CO containing only one single column
    # # (or that he has set nothing for CO and that CO has still its default
    # # value (i.e. an empty matrix of dimension 1x1)), it obviously means
    # # that he wants to use this specific covariate...
    #
    # CO has to remain as a data frame!
    #*******************************************************************************
    #
    # Updating the number of columns of CO
    if (all(is.na(CO))==FALSE) { # if CO is NOT completely empty
        nco <- ncol(CO)
    } else {                     # else, in case CO is completely empty
        # then, nco has to be set to "0"
        nco <- 0
    }
    #
    # Updating the number of columns of COt
    if (all(is.na(COt))==FALSE) { # if COt is NOT completely empty
        # (i.e. "if (ncot > 0)")
        ncot <- ncol(COt)
    } else {                     # else, in case COt is completely empty
        # then, ncot has to be set to "0"
        ncot <- 0
    }












    # Deleting entire rows of OD filled only with NAs
    # (and deleting corresponding lines in CO and COt)
    i <- 1
    while (i <= nrow(OD)) {
        if (all(is.na(OD[i,]))) {
            OD <- OD[-i,]
            if (all(is.na(CO))==FALSE) { # Checking if CO is NOT completely
                # empty and updating the covariate matrix CO as well!
                CO <- CO[-i,]
            }
            if (all(is.na(COt))==FALSE) { # Checking if COt is NOT completely
                # empty and updating the covariate matrix COt as well!
                COt <- COt[-i,]
            }
            warning(paste("/!\\ Row number",i,"of OD has been completely erased
                          because it only consisted of NAs."),sep='')
        }
        i <- i+1
    }










    # Definition of the number of rows and columns in OD
    # (And eventually update of nr)
    nr <- nrow(OD)
    nc <- ncol(OD)






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







    # In case of a factor dataset OD:
    # RECODING of OD with numbers "1", "2", etc. instead of its "words"
    ODClass <- class(OD[1,1])
    #*************************************
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




















    # 0. Initial tests and manipulations on parameters ------------------------------------------------------------------------------------------------------------

    # 0.1 Testing "regr" effectively either "mlogit", "lm" or "lrm" -----------------------------------------------------------------------------------------------
    if ( (regr != "mlogit") & (regr != "lm") & (regr != "lrm")) {
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
    if (np+1+nf>nc)
        stop("/!\\ You have to choose lower value for np and nf. Your selected
             np and nf are too large to fit the dimensions of your data matrix")





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









    # 0.9 Testing effectively mi.return == 1 or mi.retunr ==2 -----------------------------------------------------------------------------------------------------
    if (mi.return!=1 & mi.return!=2)
        stop("/!\\ mi.return can only take the values of '1' or '2'")





    # 0.10 Taking the absolute value of the parameter "noise" -----------------------------------------------------------------------------------------------------
    # Since "noise" is the variance of the elements of the final vector pred, it
    # can't be negative
    noise <- abs(noise)
















    # 1. Analysis of OD and creation of matrices ORDER, ORDER2 and ORDER3 -----------------------------------------------------------------------------------------

    # Coding of missing data in function of the length of the gap
    # ORDER: matrix of the same size of OD giving the imputation order of each
    # MD (0 for observed data and 1 for MD) --> there will be 1 everywhere there
    # is MD and 0 everywhere else
    # ORDER2: matrix of the same size of OD numbering MD into each gap (for
    # example from 1 to 3 for a gap of length 3)
    # ORDER3: matrix of the same size of OD replacing each MD by the length of
    # the gap it belongs to

    # Creation of matrix ORDER
    ORDER <- matrix(0,nr,nc) # initialization of matrix ORDER with 0 everywhere
    SEL <- is.na(OD)==TRUE   # creation of matrix SEL, constituted of TRUE where
    # there is MD in OD and of FALSE everywhere else
    ORDER[SEL] <- 1          # setting some 1 in ORDER at the location where in
    # SEL we have some TRUE





    # Creation of vector InitGapSize (i.e. a vector containing the size of the
    # initial gaps of each line)
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


    # Creation of vector TermGapSize (i.e. a vector containing the size of the
    # terminal gaps of each line)
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


    # Updating of ORDER with "0" on every external NAs
    # (The purpose of this modification of ORDER is that we don't take into
    # account external NAs at this moment of the program.
    # We will first impute internal gaps and consider external gaps further
    # (as far as nfi and npt are greater than 0))
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
    ORDER2 <- ORDER                         # initially both ORDER2 and
    ORDER3 <- ORDER                         # ORDER3 are equal to ORDER

    for (i in 1:nr){                            # in matrix ORDER, we go line
        # by line...
        for (j in 2:nc){                        # ... from column to column
            # (actually exactly as we
            # read a book) (and /!\
            # beginning from column 2!)
            if (ORDER[i,j-1]==1 & ORDER[i,j]==1){   # if the previous value
                # of ORDER is equal to 1
                # and its current value
                # is also equal to 1,
                # then...
                ORDER2[i,j] <- ORDER2[i,j-1]+1      # ... construction of
                # ORDER2 for this
                # iteration: the current
                # (j) corresponding
                # value of ORDER2 is
                # assigned by its
                # previous (j-1) value
                # incremented by 1
                # ... construction of ORDER3 for this iteration: all the values
                # in ORDER3 from the beginning of the gap (j-(ORDER2[i,j]-1)) up
                # to the current (j) location are assigned by the current (j)
                # corresponding value in ORDER2
                ORDER3[i,(j-(ORDER2[i,j]-1)):j] <- ORDER2[i,j]
            }
        }
    }

    MaxGap <- max(max(ORDER2)) # renders us the size of the greatest gap in OD














    if (max(ORDER)!=0) {

        # 2. Computation of the order of imputation of each MD (i.e. updating of matrix ORDER) --------------------------------------------------------------------

        # 2.1. Model 1: use of previous and future observations ---------------------------------------------------------------------------------------------------
        if (np>0 & nf>0){

            ord <- integer(MaxGap)          # initialization of the vector "ord"

            # Creation of the longest vector of the matrix
            ord[1] <- 1
            iter_even <-0
            iter_uneven <- 0
            for (i in 2:MaxGap) {
                if (i%%2==0) {
                    shift <- MaxGap -2 -3*iter_even
                    iter_even <- iter_even + 1
                } else {
                    shift <- -1 -iter_uneven
                    iter_uneven <- iter_uneven + 1
                }
                index <- i + shift
                ord[index] <- i
            }

            ifelse(MaxGap%%2==0, ord<-ord, ord<-rev(ord)) # reverse the order of
            # ord in case we are
            # in an even-first
            # case (that is in the
            # case of an uneven
            # MaxGap)


            # Creation of every shorter vector based on "ord" (taking some parts
            # of "ord")
            for (i in 1:nr){
                j <- 1
                while (j<=nc){
                    if (ORDER3[i,j] != 0){           # meeting a value in ORDER3
                        if (ORDER3[i,j]%%2 == 0){
                            ORDER[i,j:(j+ORDER3[i,j]-1)] <- ord[ (floor(MaxGap/2)-ORDER3[i,j]/2+1) : (floor(MaxGap/2)+ORDER3[i,j]/2) ]
                        }
                        else{
                            ORDER[i,j:(j+ORDER3[i,j]-1)] <- ord[ (floor(MaxGap/2)-floor(ORDER3[i,j]/2)+1) : (floor(MaxGap/2)+ceiling(ORDER3[i,j]/2)) ]
                        }
                        j <- j+ORDER3[i,j]+1
                    }
                    else{
                        j <- j+1
                    }
                }
            }
        }




        # 2.2. Model 2: use of previous observations only ---------------------------------------------------------------------------------------------------------
        if (np>0 & nf==0){            # Verifying that we are in case of model 2
            for (i in 1:nr){          # Beginning from row 1, we will go row by
                # row in the matrix
                j <- 1                # Setting the counter of the columns to 1
                while (j<=nc){              # We will cover the columns from
                    # left to right (until we reach the
                    # end of the row)
                    if (ORDER3[i,j]>0){                 # If we meet a component
                        # of ORDER3, it means
                        # that we are entering a
                        # gap of length equal to
                        # this component
                        numb <- ORDER3[i,j]              # We then store this
                        # size of gap in "numb"
                        ord <- c((MaxGap-numb+1):MaxGap) # We create a vector
                        # "ord" of length numb
                        # but going from where
                        # we are in the
                        # "residual part to
                        # fill" to MaxGap
                        ORDER[i,j:(j+numb-1)] <- ord     # We then insert the
                        # vector "ord" at the
                        # location where it has
                        # to go on the current
                        # row of the final
                        # matrix ORDER
                        j <- j+numb+1       # We increment the counter of the
                        # columns in order to skip the
                        # entire gap ("+numb") and to skip
                        # the directly next position as well
                        # ("+1") (because we know that the
                        # following location just after a
                        # gap won't be another gap
                        # (inevitably!)), before continuing
                        # to analyze the rest of the current
                        # row
                    }
                    else {                  # Otherwise...
                        j <- j+1            # ... we just simply increment the
                        # counter of the columns to go to
                        # the next column during the next
                        # iteration
                    }
                }
            }
        }




        # 2.3. Model 3: use of future observations only -----------------------------------------------------------------------------------------------------------
        if (np==0 & nf>0){                  # Verifying that we are effectively
            # in case of model 3
            for (i in 1:nr){                # Beginning from row 1, we will go
                # row by row in the matrix
                j <- nc                       # Setting the counter of
                # the columns to the end (i.e. the
                # total number of columns "nc")
                while (j>=1){                 # We will cover the columns from
                    # right to left (until we reach
                    # the beginning of the row)
                    if (ORDER3[i,j]>0){                  # If we meet a
                        # component of ORDER3,
                        # it means that we are
                        # entering a gap of
                        # length equal to this
                        # component
                        numb <- ORDER3[i,j]               # We then store this
                        # size of gap in
                        # "numb"
                        ord <- c(MaxGap:(MaxGap-numb+1))  # We create a vector
                        # "ord" of length numb
                        # but going from
                        # MaxGap to where we
                        # are in the "residual
                        # part to fill"
                        ORDER[i,(j-numb+1):j] <- ord      # We then insert the
                        # vector "ord" at the
                        # location where it
                        # has to go (that is
                        # "j-numb+1" elements
                        # more to the left) on
                        # the current row of
                        # the final matrix
                        # ORDER
                        j <- j-numb-1   # We decrement the counter of the
                        # columns in order to skip the entire
                        # gap going to the left ("-numb") and to
                        # skip the directly next position as
                        # well (still going to the left) ("-1")
                        # (because we know that the previous
                        # location just before a gap won't be
                        # another gap (inevitably!)), before
                        # continuing to analyze the "rest" of
                        # the current row
                    }
                    else {              # Otherwise...
                        j <- j-1        # ... we just simply decrement the
                        # counter of the columns to go to the
                        # previous column (i.e. one more column
                        # on the left) during the next iteration
                    }
                }
            }
        }











        # Updating ORDER with "0" on every NAs belonging to a Specially Located
        # Gap (SLG)
        # (The purpose of this modification of ORDER is that we don't take into
        # account SLG NAs at this moment of the program.
        # We will first impute internal gaps, external gaps and consider SLG at
        # the very end
        # (as far as some SLG have been detected)

        # 6.1 Creation of ORDERSLG (ORDERSLGLeft and ORDERSLGRight)

        # Initialization of matrix in which we will store the SLG
        ORDERSLG <- matrix(0,nrow=nr,ncol=nc)

        # Initialization of the range in which SLG could be found
        tempMinGapLeft <- matrix(0,nrow=nr,ncol=nc)
        tempMaxGapLeft <- matrix(0,nrow=nr,ncol=nc)
        tempMinGapRight <- matrix(0,nrow=nr,ncol=nc)
        tempMaxGapRight <- matrix(0,nrow=nr,ncol=nc)

        for (i in 1:nr) {                # we will go through each line of ORDER

            if (np > 1) {                # if np > 1, it may be possible that
                # SLG on the left-hand side of OD exist
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

            if (nf > 1) {               # if nf > 1, it may be possible that SLG
                # on the right-hand side of OD exist
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

        # Extracting extrema from tempMinGapLeft, tempMaxGapLeft,
        # tempMinGapRight and tempMaxGapRight
        # And creation of ORDERSLGLeft and ORDERSLGRight (matrices for both
        # groups of SLG (i.e. one on the left- and the other one on the right-
        # hand side of the matrix ORDERSLG)
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

        # /!\ Final version of the matrix ORDER that we use through point 3.1 to
        # 3.3 of the program
        ORDER <- ORDER - ORDERSLGLeft - ORDERSLGRight













        # 2.4. Creation of matrices REFORD ------------------------------------------------------------------------------------------------------------------------
        # The purpose of this section is to accelerate part 3.3 in which we
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
        for(order in 1:MaxGap) {
            assign(paste("REFORD_",order,sep=''),matrix(nrow=0,ncol=2))
        }

        ORDERinit <- ORDER

        for (i in 1:nr) {
            for (j in 1:nc) {
                if (ORDER[i,j] > 0) {

                    # Updating ORDER so that it becomes a matrix with positive
                    # values going from 1 to MaxGap
                    ORDER[i,j] <- ORDER[i,j] - (min(ORDERinit[ORDERinit!=0]) - 1)

                    # Collecting the coordinates for each value of "order"
                    coord <- t(matrix(c(i,j)))

                    tempObject = get(paste0("REFORD_",ORDER[coord]))
                    update <- rbind(tempObject,coord)

                    assign( (paste("REFORD_",ORDER[coord],sep='')),update )
                }
            }
        }








    }
















    # Initialization of the matrix in which we are going to store the results of
    # the multiple imputations
    RESULT <- cbind(replicate(nr,0),OD)


    # Beginning of the multiple imputation (imputing "mi" times)
    for (o in 1:mi) {




        # Trying to catch the potential singularity error (part 1/2)
        # (comment this part of code (as well as the one at the very end of
        # seqimpute.R) to debug more easily and see the full message of any
        # occuring error)
        #********************************************************************************************************************************
        tryCatch( # Trying to catch the potential singularity error
            # in order to display a more accurate comment on it

            {
                #************************************************************************************************************************




                if (max(ORDER)!=0) {
                    # Otherwise if there is only 0 in ORDER,
                    # there is no need to impute internal gaps
                    # and we directly jump to the imputation of
                    # external gaps (i.e. points 4. and 5.)








                    # 3. Imputation using a specific model --------------------------------------------------------------------------------------------------------

                    # 3.1 Building of the data matrix CD for the computation of the model -------------------------------------------------------------------------

                    for (order in 1:MaxGap){ # /!\ "order" corresponds to the
                        # values of the components of ORDER (i.e. the number
                        # of the iteration, the order in which the
                        # values are going to be imputed)

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
                        if (np>0 & nf==0) {              # only PAST VIs do
                                                         # exist
                            # initialisation of matrix CDp
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
                                # /!\ current
                                # pointer on column
                                # is thus:
                                # "j-frameSize
                                # +np+1"
                                # Past VIs
                                CDp[t1:t2,] <- OD[,(j-frameSize+1):(j-frameSize+np)]

                                # Eventually considering time-dependent
                                # covariates
                                if (ncot > 0) {
                                    COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                                    for (d in 1:(ncot/nc)) {
                                        COttemp <- cbind(COttemp, COt[,(j-frameSize+np+1) + (d-1)*nc])
                                    }
                                    COtselected[t1:t2,] <- COttemp
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
                                    if ( (j-frameSize+np+2) <= nc ) {
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


                            # Concatenating CD
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


                            # Only FUTURE
                        } else if (np==0 & nf>0) {
                            # only FUTURE VIs do exist
                            # initialisation of matrix CDf
                            CDf <- matrix(NA, nrow=nr*ud, ncol=nf)

                            if (ncot > 0) {
                                # initialisation of matrix COtselected
                                COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
                            }

                            # initialisation of matrix CDdb
                            # (for past distribution
                            # analysis) (Distribution Before)
                            if (pastDistrib) {
                                CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
                                db <- matrix(NA, nrow=nr, ncol=k)
                                # submatrix of CDdb:
                                # CDdb is composed of
                                # ud matrix db on top
                                # of each other
                            }

                            # initialisation of matrix CDda
                            # (for future distribution
                            # analysis) (Distribution After)
                            if (futureDistrib) {
                                CDda <- matrix(NA, nrow=nr*ud, ncol=k)
                                # CDda has same
                                # dimensions as
                                # CDdb
                                da <- matrix(NA, nrow=nr, ncol=k)
                                # da has same
                                # dimensions as
                                # db
                            }


                            for (j in frameSize:nc){
                                t1 <- (nr*(iter-1)+1)
                                t2 <- nr*iter
                                # VD
                                CD[t1:t2,1] <- OD[,j-frameSize+np+1]
                                # /!\ current
                                # pointer on column
                                # is thus: "j-
                                # frameSize+np+1"
                                # Future VIs
                                CDf[t1:t2,] <- OD[,(j-nf+1):j]

                                # Eventually considering time-dependent
                                # covariates
                                if (ncot > 0) {
                                    COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                                    for (d in 1:(ncot/nc)) {
                                        COttemp <- cbind(COttemp, COt[,(j-frameSize+np+1) + (d-1)*nc])
                                    }
                                    COtselected[t1:t2,] <- COttemp
                                }

                                # Past distribution (i.e. Before)
                                if (pastDistrib) {
                                    ODt <- t(OD)
                                    ODt <- as.data.frame(ODt)
                                    tempOD <- lapply(ODt[(1:(j-frameSize+np)),], factor, levels=c(1:k,NA), exclude=NULL)
                                    # because: j-frameSize+np+1 - 1
                                    # = j-frameSize+np

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
                                    # j-frameSize+np+1 + 1 = j-frameSize+np+2

                                    da_list <- lapply(tempOD,summary)
                                    da_matrix <- do.call(rbind,da_list)
                                    CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np+2):nc)
                                }

                                iter <- iter+1
                            }

                            # Concatenating CD
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
                            # any covariate)
                            # simply continue with the current CD

                            # Eventually concatenating CD with COtselected (the
                            # matrix containing the current time-dependent
                            # covariates)
                            # Checking if COt is NOT completely empty
                            if (ncot > 0) {
                                # Concatenating CD and COtselected into CD
                                CD <- cbind(CD, as.data.frame(COtselected))
                            }


                            # PAST and FUTURE
                        } else {
                            # meaning np>0 and nf>0 and that, thus,
                            # PAST as well as FUTURE VIs do exist
                            # initialisation of matrices CDp and CDf
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
                                # CDda has same
                                # dimensions as
                                # CDdb
                                da <- matrix(NA, nrow=nr, ncol=k)
                                # da has same
                                # dimensions as
                                # db
                            }


                            for (j in frameSize:nc){
                                t1 <- (nr*(iter-1)+1)
                                t2 <- nr*iter
                                # VD
                                CD[t1:t2,1] <- OD[,j-frameSize+np+1+shift]
                                # /!\ current
                                # pointer on
                                # column is
                                # thus: "j-
                                # frameSize+
                                # np+1+shift"
                                # Past VIs
                                CDp[t1:t2,] <- OD[,(j-frameSize+1):(j-frameSize+np)]
                                # Future VIs
                                CDf[t1:t2,] <- OD[,(j-nf+1):j]

                                # Eventually considering time-dependent
                                # covariates
                                if (ncot > 0) {
                                    COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                                    for (d in 1:(ncot/nc)) {
                                        COttemp <- cbind(COttemp, COt[,(j-frameSize+np+1+shift) + (d-1)*nc])
                                    }
                                    COtselected[t1:t2,] <- COttemp
                                }

                                # Past distribution (i.e. Before)
                                if (pastDistrib) {
                                    ODt <- t(OD)
                                    ODt <- as.data.frame(ODt)
                                    tempOD <- lapply(ODt[(1:(j-frameSize+np+shift)),], factor, levels=c(1:k,NA), exclude=NULL)
                                    # because:
                                    # j-frameSize+np+1+shift - 1 = j-frameSize
                                    # +np+shift

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
                                    # j-frameSize+np+1+shift + 1 = j-frameSize+
                                    # np+shift+2

                                    da_list <- lapply(tempOD,summary)
                                    da_matrix <- do.call(rbind,da_list)
                                    CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np+shift+2):nc)
                                }

                                iter <- iter+1
                            }

                            # Concatenating CD
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
                            # any covariate)
                            # simply continue with the current CD

                            # Eventually concatenating CD with COtselected (the
                            # matrix containing the current time-dependent
                            # covariates)
                            # Checking if COt is NOT completely empty
                            if (ncot > 0) {
                                # Concatenating CD and COtselected into CD
                                CD <- cbind(CD, as.data.frame(COtselected))
                            }


                        }











                        # 3.2 Computation of the model (Dealing with the LOCATIONS of imputation) -----------------------------------------------------------------




                        #*********************************************************************************
                        # Initially "computeModel.R"
                        #*********************************************************************************


                        # ==>> Manipulation of parameters

                        # Conversion of CD in a data frame
                        CD <- as.data.frame(CD)

                        # Transformation of the names of the columns of CD
                        # (called V1, V2, ..., "VtotV")
                        colnames(CD) <- paste("V", 1:ncol(CD), sep = "")










                        if (regr == "mlogit") {
                            ## Case of MULTINOMIAL REGRESSION MODEL

                            # Linking the package mlogit

                            # By default, every column of CD are of class
                            # "numeric".
                            # Thus, there is no need to convert the columns
                            # containing the distribution data to
                            # class "numeric".
                            # Moreover the class of the covariates columns at
                            # the very end are ALREADY set correctly and we
                            # don't need to update them.
                            # On the other hand, the first columns of CD
                            # (1 up to 1+np+nf) have to be of class "factor"
                            # (because they are the columns containing our
                            # categorical data coming from OD).

                            # Transformation of the first columns (i.e. the
                            # categorical values) of CD (column 1 up to column
                            # 1+np+nf) into factor
                            CD[,(1:(1+np+nf))] <- lapply(CD[,(1:(1+np+nf))],factor, levels=c(1:k,NA), exclude=NULL)


                            # Dataframe for mlogit
                            NCD <- mlogit.data(CD, varying=NULL, choice="V1", shape="wide")


                            # Computation of the multinomial model
                            if(totV==1){
                                # First case is evaluated aside
                                reglog_3 <- mlogit(V1~0, data=NCD, reflevel="1")
                            }

                            if(totV>1){
                                # creation of "V2" ... "VtotV" (to use them in
                                # the formula)
                                factors_character <- paste("V", 2:totV, sep = "")
                                # Transformation of this object from character
                                # to vector (in order to be able to access its
                                # components)
                                factors <- as.vector(factors_character)
                                # Creation of a specific formula according to
                                # the value of totV
                                fmla <- as.formula(paste("V1~0|", paste(factors, collapse="+")))
                                reglog_3 <- mlogit(fmla, data=NCD, reflevel="1")
                            }









                        } else if (regr == "lm") {
                            ## Case of LINEAR REGRESSION MODEL

                            # Since we are performing a linear regression, each
                            # element of CD are numbers and have then to remain
                            # as class "numeric" (we don't have to perform some
                            # class adjustment as for the case of the creation
                            # of a multinomial model).

                            # Computation of the linear regression model
                            if(totV==1){
                                reglog_3 <- lm(V1~0, data=CD)  # first case is
                                # evaluated aside
                            }

                            if(totV>1){
                                # creation of "V2" ... "VtotV" (to use them in
                                # the formula)
                                factors_character <- paste("V", 2:totV, sep = "")
                                # Transformation of this object from character
                                # to vector (in order to be able
                                # to access its components)
                                factors <- as.vector(factors_character)
                                # Creation of a specific formula according to
                                # the value of totV
                                fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                                reglog_3 <- lm(fmla, data=CD)
                            }











                        } else { # meaning (regr == "lrm")
                            ## Case of ORDINAL REGRESSION MODEL

                            # Linking to the package rms to use the function
                            # "lrm"

                            # Since we are performing an ordinal regression,
                            # each element of CD are numbers and have then to
                            # remain as class "numeric" (we don't have to
                            # perform some class adjustment as for the case of
                            # the creation of a multinomial model).

                            # Computation of the ordinal model
                            if(totV==1){
                                reglog_3 <- lrm(V1~0, data=CD) # first case is
                                # evaluated aside
                            }

                            if(totV>1){
                                # creation of "V2" ... "VtotV" (to use them in
                                # the formula)
                                factors_character <- paste("V", 2:totV, sep = "")
                                # Transformation of this object from character
                                # to vector (in order to be able
                                # to access its components)
                                factors <- as.vector(factors_character)
                                # Creation of a specific formula according to
                                # the value of totV
                                fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                                reglog_3 <- lrm(fmla, data=CD)
                            }

                        }




                        #*********************************************************************************
                        #*********************************************************************************













                        # 3.3 Imputation using the just created model (Dealing with the actual VALUES to impute) --------------------------------------------------

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
                        tempObject = get(paste0("REFORD_",order))
                        REFORD <- as.matrix(tempObject)
                        if (ncol(REFORD) == 1) {
                            REFORD <- t(REFORD)
                        }
                        nr_REFORD <- nrow(REFORD)



                        if (np>0 & nf==0) {             # only PAST VIs do exist
                            for (u in 1:nr_REFORD) {
                                i <- REFORD[u,1]
                                # taking out the first coordinate
                                # (row number in ORDER) from REFORD
                                j <- REFORD[u,2]
                                # taking out the second coordinate
                                # (column number in ORDER) from REFORD

                                CDi <- matrix(NA, nrow=k, ncol=1)

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
                                CDi[,(1:(1+np+nf))] <- lapply(CDi[,(1:(1+np+nf))],factor, levels=c(1:k,NA), exclude=NULL)

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


                                    #*******************************************************************
                                    # Initially "imputeValue.R"
                                    #*******************************************************************

                                    if (regr == "mlogit") {
                                        ## Case of MULTINOMIAL REGRESSION MODEL


                                        pred <- predict(reglog_3,newdata=CDi)
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
                                        sel <- which(pred>=alea)
                                        # Corresponding example value returned
                                        # in sel:
                                        # 1 2
                                        # 1 2
                                        #















                                    } else if (regr == "lm") {
                                        ## Case of LINEAR REGRESSION MODEL

                                        # Since we are performing a linear
                                        # regression, each element of CDi are
                                        # numbers and have to be considered as
                                        # class "numeric"
                                        CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                        pred <- predict(reglog_3, CDi)
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













                                    } else { # meaning (regr == "lrm")
                                        ## Case of ORDINAL REGRESSION MODEL

                                        # Since we are performing an ordinal
                                        # regression, each element of CDi are
                                        # numbers and have to be considered as
                                        # class "numeric"
                                        CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                        pred <- predict(reglog_3, CDi, type="fitted.ind")
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



                                    }

                                    #*******************************************************************
                                    #*******************************************************************



                                    ODi[i,j] <- sel[1]
                                }
                            }

                        } else if (np==0 & nf>0) {  # only FUTURE VIs do exist
                            for (u in 1:nr_REFORD) {
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
                                CDi[,(1:(1+np+nf))] <- lapply(CDi[,(1:(1+np+nf))],factor, levels=c(1:k,NA), exclude=NULL)
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
                                if (max(is.na(CDi[1,2:totV]))==0){


                                    #*******************************************************************
                                    # Initially "imputeValue.R"
                                    #*******************************************************************

                                    if (regr == "mlogit") {
                                        ## Case of MULTINOMIAL REGRESSION MODEL


                                        pred <- predict(reglog_3,newdata=CDi)
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
                                        sel <- which(pred>=alea)
                                        # Corresponding example value returned
                                        # in sel:
                                        # 1 2
                                        # 1 2
                                        #















                                    } else if (regr == "lm") {
                                        ## Case of LINEAR REGRESSION MODEL

                                        # Since we are performing a linear
                                        # regression, each element of CDi are
                                        # numbers and have to be considered as
                                        # class "numeric"
                                        CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                        pred <- predict(reglog_3, CDi)
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













                                    } else { # meaning (regr == "lrm")
                                        ## Case of ORDINAL REGRESSION MODEL

                                        # Since we are performing an ordinal
                                        # regression, each element of CDi are
                                        # numbers and have to be considered as
                                        # class "numeric"
                                        CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                        pred <- predict(reglog_3, CDi, type="fitted.ind")
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



                                    }

                                    #*******************************************************************
                                    #*******************************************************************


                                    ODi[i,j] <- sel[1]
                                }
                            }

                        } else { # meaning np>0 and nf>0 and that,
                            # thus, PAST as well as FUTURE VIs
                            # do exist
                            for (u in 1:nr_REFORD) {
                                i <- REFORD[u,1]
                                # taking out the first coordinate
                                # (row number in ORDER) from REFORD
                                j <- REFORD[u,2]
                                # taking out the second coordinate
                                # (column number in ORDER) from REFORD

                                CDi <- matrix(NA, nrow=k, ncol=1)

                                # Matrix for past values
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
                                CDi[,(1:(1+np+nf))] <- lapply(CDi[,(1:(1+np+nf))],factor, levels=c(1:k,NA), exclude=NULL)
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


                                    #*******************************************************************
                                    # Initially "imputeValue.R"
                                    #*******************************************************************

                                    if (regr == "mlogit") {
                                        ## Case of MULTINOMIAL REGRESSION MODEL


                                        pred <- predict(reglog_3,newdata=CDi)
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
                                        sel <- which(pred>=alea)
                                        # Corresponding example value returned
                                        # in sel:
                                        # 1 2
                                        # 1 2
                                        #















                                    } else if (regr == "lm") {
                                        ## Case of LINEAR REGRESSION MODEL

                                        # Since we are performing a linear
                                        # regression, each element of CDi are
                                        # numbers and have to be considered as
                                        # class "numeric"
                                        CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                        pred <- predict(reglog_3, CDi)
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













                                    } else { # meaning (regr == "lrm")
                                        ## Case of ORDINAL REGRESSION MODEL

                                        # Since we are performing an ordinal
                                        # regression, each element of CDi are
                                        # numbers and have to be considered as
                                        # class "numeric"
                                        CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                        pred <- predict(reglog_3, CDi, type="fitted.ind")
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



                                    }

                                    #*******************************************************************
                                    #*******************************************************************


                                    ODi[i,j] <- sel[1]
                                }
                            }
                        }

                    }


                }



















                # 4. Imputing initial NAs -------------------------------------------------------------------------------------------------------------------------

                if ((nfi != 0) & (MaxInitGapSize != 0)) {
                    # # we only impute the initial gaps if nfi > 0



                    # 4.1.-2. Creation of ORDERI ------------------------------------------------------------------------------------------------------------------

                    # Creation of matrix ORDERI
                    ORDERI <- matrix(0,nr,nc)
                    for (i in 1:nr) {
                        if (InitGapSize[i]!=0) {
                            ORDERI[i,1:InitGapSize[i]] <- c(MaxInitGapSize:(MaxInitGapSize+1-InitGapSize[i]))
                        } else {
                            next
                        }
                    }

                    # In a similar manner to part 2.4., we go here one single
                    # time through a reduced version of ORDERI and we create
                    # MaxInitGapSize REFORDI matrices collecting the coordinates
                    # of each corresponding values in ORDERI which are greater
                    # than 0


                    # Initialization of the REFORDI matrices
                    for(order in 1:MaxInitGapSize) {
                        assign(paste("REFORDI_",order,sep=''),matrix(nrow=0,ncol=2))
                    }
                    # Collecting the coordinates for each value of "order"
                    for (i in 1:nr) {
                        for (j in MaxInitGapSize:1) {
                            if (ORDERI[i,j] > 0) {

                                coord <- t(matrix(c(i,j)))

                                tempObject = get(paste0("REFORDI_",ORDERI[coord]))
                                update <- rbind(tempObject,coord)

                                assign( (paste("REFORDI_",ORDERI[coord],sep='')),update )

                            }
                        }
                    }






                    # 4.3. Imputation using a specific model ------------------------------------------------------------------------------------------------------

                    # 4.3.1 Building of the data matrix CD for the computation of the model -----------------------------------------------------------------------

                    # For Initial Gaps
                    # we will impute single NA after single NA going from the
                    # center of OD towards its very left border
                    # We here only take observations from the FUTURE into
                    # account --> nfi
                    # But we will create one single regression
                    # model used for the imputation of all
                    # the NAs belonging to an "initial gap"
                    frameSize <- 1+nfi
                    ud <- nc-frameSize+1

                    CD <- matrix(NA, nrow=nr*ud, ncol=1)

                    # initialisation of matrix CDf
                    CDf <- matrix(NA, nrow=nr*ud, ncol=nfi)

                    if (ncot > 0) {
                        # initialisation of matrix COtselected
                        COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
                    }

                    # initialisation of matrix CDda (for future distribution
                    # analysis)
                    # (Distribution After)
                    if (futureDistrib) {
                        CDda <- matrix(NA, nrow=nr*ud, ncol=k)
                        da <- matrix(NA, nrow=nr, ncol=k)
                        # submatrix of CDda: CDda is
                        # composed of ud matrix da on
                        # top of each other
                    }

                    iter <- 1

                    for (j in frameSize:nc) {

                        t1 <- (nr*(iter-1)+1)
                        t2 <- nr*iter

                        # VD
                        CD[t1:t2,1] <- OD[,j-frameSize+1]
                        # /!\ current pointer on
                        # column is thus:
                        # "j-frameSize+1"

                        # Future VIs
                        CDf[t1:t2,] <- OD[,(j-nfi+1):j]

                        # Eventually considering time-dependent covariates
                        if (ncot > 0) {
                            COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                            for (d in 1:(ncot/nc)) {
                                COttemp <- cbind(COttemp, COt[,(j-frameSize+1) + (d-1)*nc])
                            }
                            COtselected[t1:t2,] <- COttemp
                        }

                        # Future distribution (i.e. After)
                        if (futureDistrib) {
                            ODt <- t(OD)
                            ODt <- as.data.frame(ODt)
                            tempOD <- lapply(ODt[((j-frameSize+2):nc),], factor, levels=c(1:k,NA), exclude=NULL)
                            da_list <- lapply(tempOD,summary)
                            da_matrix <- do.call(rbind,da_list)
                            CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+2):nc)
                        }

                        iter <- iter+1
                    }

                    # Concatenating CD
                    CD <- cbind(CD, CDf)

                    if (futureDistrib) {
                        CD <- cbind(CD, CDda)
                    }

                    # Conversion of CD into a data frame
                    CD <- as.data.frame(CD)

                    # Eventually concatenating CD with COs (the matrix
                    # containing the covariates)
                    if (all(is.na(CO))==FALSE) {
                        # Checking if CO is NOT completely
                        # empty
                        # Creation of the stacked covariates matrix for 3.1
                        COs <- do.call("rbind", rep(list(CO), ud))
                        # Concatenating CD and COs into CD
                        CD <- cbind(CD, COs)
                    }
                    # Else, in case CO is empty (i.e. we don't consider any
                    # covariate)
                    # simply continue with the current CD

                    # Eventually concatenating CD with COtselected (the
                    # matrix containing the current time-dependent
                    # covariates)
                    # Checking if COt is NOT completely empty
                    if (ncot > 0) {
                        # Concatenating CD and COtselected into CD
                        CD <- cbind(CD, as.data.frame(COtselected))
                    }





                    # 4.3.2 Computation of the model (Dealing with the LOCATIONS of imputation) -------------------------------------------------------------------


                    #*********************************************************************************
                    # Initially "computeModel.R"
                    #*********************************************************************************





                    # ==>> Manipulation of parameters

                    # Conversion of CD in a data frame
                    CD <- as.data.frame(CD)

                    # Transformation of the names of the columns of CD
                    # (called V1, V2, ..., "VtotVi")
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
                        CD[,(1:(1+0+nfi))] <- lapply(CD[,(1:(1+0+nfi))],factor, levels=c(1:k,NA), exclude=NULL)


                        # Dataframe for mlogit
                        NCD <- mlogit.data(CD, varying=NULL, choice="V1", shape="wide")


                        # Computation of the multinomial model
                        if(totVi==1){
                            # First case is evaluated aside
                            reglog_4 <- mlogit(V1~0, data=NCD, reflevel="1")
                        }

                        if(totVi>1){
                            # creation of "V2" ... "VtotVi" (to use them in the
                            # formula)
                            factors_character <- paste("V", 2:totVi, sep = "")
                            # Transformation of this object from character to
                            # vector (in order to be able to access its
                            # components)
                            factors <- as.vector(factors_character)
                            # Creation of a specific formula according to the
                            # value of totVi
                            fmla <- as.formula(paste("V1~0|", paste(factors, collapse="+")))
                            reglog_4 <- mlogit(fmla, data=NCD, reflevel="1")
                        }









                    } else if (regr == "lm") {
                        ## Case of LINEAR REGRESSION MODEL

                        # Since we are performing a linear regression, each
                        # element of CD are numbers and have then to remain as
                        # class "numeric" (we don't have to perform some class
                        # adjustment as for the case of the creation of a
                        # multinomial model).

                        # Computation of the linear regression model
                        if(totVi==1){
                            reglog_4 <- lm(V1~0, data=CD)
                            # first case is evaluated aside
                        }

                        if(totVi>1){
                            # creation of "V2" ... "VtotVi" (to use them in the
                            # formula)
                            factors_character <- paste("V", 2:totVi, sep = "")
                            # Transformation of this object from character to
                            # vector (in order to be able
                            # to access its components)
                            factors <- as.vector(factors_character)
                            # Creation of a specific formula according to the
                            # value of totVi
                            fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                            reglog_4 <- lm(fmla, data=CD)
                        }











                    } else { # meaning (regr == "lrm")
                        ## Case of ORDINAL REGRESSION MODEL

                        # Linking to the package rms to use the function "lrm"

                        # Since we are performing an ordinal regression, each
                        # element of CD are numbers and have then to remain as
                        # class "numeric" (we don't have to perform some
                        # class adjustment as for the case of the creation of a
                        # multinomial model).

                        # Computation of the ordinal model
                        if(totVi==1){
                            reglog_4 <- lrm(V1~0, data=CD)
                            # first case is evaluated aside
                        }

                        if(totVi>1){
                            # creation of "V2" ... "VtotVi" (to use them in the
                            # formula)
                            factors_character <- paste("V", 2:totVi, sep = "")
                            # Transformation of this object from character to
                            # vector (in order to be able
                            # to access its components)
                            factors <- as.vector(factors_character)
                            # Creation of a specific formula according to the
                            # value of totVi
                            fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                            reglog_4 <- lrm(fmla, data=CD)
                        }

                    }


                    #*********************************************************************************
                    #*********************************************************************************







                    # 4.3.3 Imputation using the just created model (Dealing with the actual VALUES to impute) ----------------------------------------------------

                    # Conversion of ODi from data.frame to matrix
                    ODi <- as.matrix(ODi)

                    # Only FUTURE VIs are useful
                    for (order in 1:MaxInitGapSize){

                        # Analysing the value of parameter available
                        if (available==TRUE){
                            # we take the previously imputed data
                            # into account
                            LOOKUP <- ODi
                        } else {
                            # that is available == FALSE and thus we don't take
                            # the previously imputed data into account
                            LOOKUP <- OD
                        }


                        # Assigning the current "REFORDI_order" matrix to the
                        # variable matrix REFORDI
                        # (according to the current value of "order")
                        tempObject = get(paste0("REFORDI_",order))
                        REFORDI <- as.matrix(tempObject)
                        if (ncol(REFORDI) == 1) {
                            REFORDI <- t(REFORDI)
                        }
                        nr_REFORDI <- nrow(REFORDI)

                        for (u in 1:nr_REFORDI) {
                            i <- REFORDI[u,1]
                            # taking out the first coordinate (row
                            # number in ORDER) from REFORDI
                            j <- REFORDI[u,2]
                            # taking out the second coordinate
                            # (column number in ORDER) from REFORDI

                            CDi <- matrix(NA,nrow=k,ncol=1)

                            # Matrix for future values
                            vect <- LOOKUP[i,(j+1):(j+nfi)]
                            CDfi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)

                            # Matrix for future distribution
                            if (futureDistrib) {
                                dai <- summary(factor(LOOKUP[i,(j+1):nc], levels=c(1:k), exclude=NULL))/length((j+1):nc)
                                CDdai <- matrix(dai[1:k], nrow=k, ncol=k, byrow=TRUE)
                            }

                            # Concatenating CDi
                            CDi <- cbind(CDi, CDfi)

                            if (futureDistrib) {
                                CDi <- cbind(CDi, CDdai)
                            }

                            # Conversion of CDi into a data frame
                            CDi <- as.data.frame(CDi)

                            # Type transformation of the columns of CDi
                            # The first values of CDi must be of type factor
                            # (categorical values)
                            CDi[,(1:(1+nfi))] <- lapply(CDi[,(1:(1+nfi))],factor, levels=c(1:k,NA), exclude=NULL)
                            # The last values of CDi must be of type numeric
                            # (distributions)
                            if (futureDistrib) {
                                CDi[,(1+nfi+1):ncol(CDi)] <- lapply(CDi[,(1+nfi+1):ncol(CDi)],as.numeric)
                            }
                            # Eventually concatenating CDi with COi (the matrix
                            # containing the covariates)
                            if (all(is.na(CO))==FALSE) { # Checking if CO is NOT
                                # completely empty
                                # Creation of the matrix COi used in 3.3
                                COi <- do.call(rbind, replicate(k, as.data.frame(CO[i,]), simplify=FALSE))
                                # Concatenating CDi and COi into CDi
                                CDi <- cbind(CDi, COi)
                                # Transformation of the names of the columns of
                                # CDi (called V1, V2, ..., "VtotV")
                                colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                            }
                            # Else, in case CO is empty (i.e. we don't consider
                            # any covariate) simply continue with the current
                            # CDi

                            # Eventually concatenating CDi with COtselected_i
                            # (the matrix containing the current time-dependent
                            # covariates)
                            # Checking if COt is NOT completely empty
                            if (ncot > 0) {
                                COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
                                for (d in 1:(ncot/nc)) {
                                    COtselected_i <- cbind(COtselected_i, COt[i,(j) + (d-1)*nc])
                                }
                                COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
                                # Concatenating CDi and COtselected_i into CDi
                                CDi <- cbind(CDi, COtselected_i)
                                # Transformation of the names of the columns of
                                # CDi (called V1, V2, ..., "VtotV")
                                colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                            }


                            # Check for missing-values among predictors
                            if (max(is.na(CDi[1,2:totVi]))==0){


                                #*******************************************************************
                                # Initially "imputeValue.R"
                                #*******************************************************************

                                if (regr == "mlogit") {
                                    ## Case of MULTINOMIAL REGRESSION MODEL


                                    pred <- predict(reglog_4,newdata=CDi)
                                    # Example of value returned by pred:
                                    # (Sytematically, the upper line represents
                                    # the possible categories of the variable
                                    # (here, k=2, so the possible categories are
                                    # either 1" or "2"))
                                    #            1            2
                                    # 1.000000e+00 2.111739e-22
                                    #
                                    # Cumulative pred
                                    pred <- cumsum(pred)
                                    # Corresponding example value returned by in
                                    # the "cumulative pred":
                                    # 1 2
                                    # 1 1
                                    #
                                    # Introducing a variance "noise"
                                    pred <- rnorm(length(pred),pred,noise)
                                    # Checking that the components of vector
                                    # "pred" are still included in the interval
                                    # [0,1]
                                    pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                    pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                    # Imputation
                                    alea <- runif(1)
                                    # Example value returned in "alea":
                                    # [1] 0.2610005
                                    #
                                    sel <- which(pred>=alea)
                                    # Corresponding example value returned in
                                    # sel:
                                    # 1 2
                                    # 1 2
                                    #















                                } else if (regr == "lm") {
                                    ## Case of LINEAR REGRESSION MODEL

                                    # Since we are performing a linear
                                    # regression, each element of CDi are
                                    # numbers and have to
                                    # be considered as class "numeric"
                                    CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                    pred <- predict(reglog_4, CDi)
                                    # Introducing a variance "noise"
                                    pred <- rnorm(length(pred),pred,noise)
                                    # Rounding pred to the nearest integer
                                    pred <- round(pred)
                                    # Restricting pred to its lowest value: 1
                                    pred <- ifelse(pred<1,1,pred)
                                    # Restricting pred to its highest value: k
                                    pred <- ifelse(pred>k,k,pred)
                                    sel <- pred













                                } else { # meaning (regr == "lrm")
                                    ## Case of ORDINAL REGRESSION MODEL

                                    # Since we are performing an ordinal
                                    # regression, each element of CDi are
                                    # numbers and have to
                                    # be considered as class "numeric"
                                    CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                    pred <- predict(reglog_4, CDi, type="fitted.ind")
                                    # Testing if we are in case where k=2
                                    # (if this is the case, we need to create
                                    # the second complementary probility by hand
                                    # since lrm returns only the first
                                    # probability)
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
                                    # "pred". We have then to restrict "alea"
                                    # to the last value of "pred"
                                    alea <- runif(1)
                                    if (alea > pred[length(pred)]) {
                                        alea <- pred[length(pred)]
                                    }
                                    sel <- which(pred>=alea)



                                }

                                #*******************************************************************
                                #*******************************************************************


                                ODi[i,j] <- sel[1]
                            }
                        }
                    }
                }



















                # 5. Imputing terminal NAs ------------------------------------------------------------------------------------------------------------------------

                if ((npt != 0) & (MaxTermGapSize != 0)) {
                    # we only impute the terminal
                    # gaps if npt > 0



                    # 5.1.-2. Creation of ORDERT ------------------------------------------------------------------------------------------------------------------

                    # Creation of matrix ORDERT
                    ORDERT <- matrix(0,nr,nc)
                    for (i in 1:nr) {
                        if (TermGapSize[i]!=0) {
                            ORDERT[i,(nc-TermGapSize[i]+1):nc] <- c((MaxTermGapSize+1-TermGapSize[i]):MaxTermGapSize)
                        } else {
                            next
                        }
                    }

                    # In a similar manner to part 2.4., we go here one single
                    # time through a reduced version of ORDERT and we create
                    # MaxTermGapSize REFORDT matrices collecting the coordinates
                    # of each corresponding values in ORDERT which are greater
                    # than 0


                    # Initialization of the REFORDT matrices
                    for(order in 1:MaxTermGapSize) {
                        assign(paste("REFORDT_",order,sep=''),matrix(nrow=0,ncol=2))
                    }
                    # Collecting the coordinates for each value of "order"
                    for (i in 1:nr) {
                        for (j in (nc-MaxTermGapSize+1):nc) {
                            if (ORDERT[i,j] > 0) {

                                coord <- t(matrix(c(i,j)))

                                tempObject = get(paste0("REFORDT_",ORDERT[coord]))
                                update <- rbind(tempObject,coord)

                                assign( (paste("REFORDT_",ORDERT[coord],sep='')),update )

                            }
                        }
                    }






                    # 5.3. Imputation using a specific model ------------------------------------------------------------------------------------------------------

                    # 5.3.1 Building of the data matrix CD for the computation of the model -----------------------------------------------------------------------

                    # For Terminal Gaps
                    # we will impute single NA after single NA going from the
                    # center of OD towards its very right border
                    # We here only take observations from the PAST
                    # into account --> npt
                    # But we will create one single regression
                    # model used for the imputation of all the NAs belonging to
                    # a "terminal gap"
                    frameSize <- npt+1
                    ud <- nc-frameSize+1

                    CD <- matrix(NA, nrow=nr*ud, ncol=1)

                    # initialisation of matrix CDp
                    CDp <- matrix(NA, nrow=nr*ud, ncol=npt)

                    if (ncot > 0) {
                        # initialisation of matrix COtselected
                        COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
                    }

                    # initialisation of matrix CDdb
                    # (for past distribution analysis)
                    # (Distribution Before)
                    if (pastDistrib) {
                        CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
                        db <- matrix(NA, nrow=nr, ncol=k)
                        # submatrix of CDdb: CDdb is
                        # composed of ud matrix db on
                        # top of each other
                    }

                    iter <- 1

                    for (j in frameSize:nc) {

                        t1 <- (nr*(iter-1)+1)
                        t2 <- nr*iter

                        # VD
                        CD[t1:t2,1] <- OD[,j]
                        # /!\ current pointer on column is thus: "j"

                        # Past VIs
                        CDp[t1:t2,] <- OD[,(j-npt):(j-1)]

                        # Eventually considering time-dependent covariates
                        if (ncot > 0) {
                            COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                            for (d in 1:(ncot/nc)) {
                                COttemp <- cbind(COttemp, COt[,(j) + (d-1)*nc])
                            }
                            COtselected[t1:t2,] <- COttemp
                        }

                        # Past distribution (i.e. Before)
                        if (pastDistrib) {
                            ODt <- t(OD)
                            ODt <- as.data.frame(ODt)
                            tempOD <- lapply(ODt[(1:(j-1)),], factor, levels=c(1:k,NA), exclude=NULL)
                            db_list <- lapply(tempOD,summary)
                            db_matrix <- do.call(rbind,db_list)
                            CDdb[t1:t2,] <- db_matrix[,1:k]/length(1:(j-1))
                        }

                        iter <- iter+1
                    }

                    # Concatening CD
                    CD <- cbind(CD, CDp)

                    if (pastDistrib) {
                        CD <- cbind(CD, CDdb)
                    }

                    # Conversion of CD into a data frame
                    CD <- as.data.frame(CD)

                    # Eventually concatenating CD with COs
                    # (the matrix containing the covariates)
                    if (all(is.na(CO))==FALSE) {
                        # Checking if CO is NOT completely empty
                        # Creation of the stacked covariates matrix for 3.1
                        COs <- do.call("rbind", rep(list(CO), ud))
                        # Concatenating CD and COs into CD
                        CD <- cbind(CD, COs)
                    }
                    # Else, in case CO is empty (i.e. we don't consider any
                    # covariate) simply continue with the current CD

                    # Eventually concatenating CD with COtselected (the
                    # matrix containing the current time-dependent
                    # covariates)
                    # Checking if COt is NOT completely empty
                    if (ncot > 0) {
                        # Concatenating CD and COtselected into CD
                        CD <- cbind(CD, as.data.frame(COtselected))
                    }





                    # 5.3.2 Computation of the model (Dealing with the LOCATIONS of imputation) -------------------------------------------------------------------



                    #*********************************************************************************
                    # Initially "computeModel.R"
                    #*********************************************************************************





                    # ==>> Manipulation of parameters

                    # Conversion of CD in a data frame
                    CD <- as.data.frame(CD)

                    # Transformation of the names of the columns of CD
                    # (called V1, V2, ..., "VtotVt")
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
                        # On the other hand, the first columns of CD
                        # (1 up to 1+np+nf) have to be of class "factor"
                        # (because they are the columns containing our
                        # categorical data coming from OD).

                        # Transformation of the first columns (i.e. the
                        # categorical values) of CD (column 1 up to column
                        # 1+np+nf) into factor
                        CD[,(1:(1+npt+0))] <- lapply(CD[,(1:(1+npt+0))],factor, levels=c(1:k,NA), exclude=NULL)


                        # Dataframe for mlogit
                        NCD <- mlogit.data(CD, varying=NULL, choice="V1", shape="wide")


                        # Computation of the multinomial model
                        if(totVt==1){
                            # First case is evaluated aside
                            reglog_5 <- mlogit(V1~0, data=NCD, reflevel="1")
                        }

                        if(totVt>1){
                            # creation of "V2" ... "VtotVt" (to use them in the
                            # formula)
                            factors_character <- paste("V", 2:totVt, sep = "")
                            # Transformation of this object from character to
                            # vector (in order to be able to access its
                            # components)
                            factors <- as.vector(factors_character)
                            # Creation of a specific formula according to the
                            # value of totVt
                            fmla <- as.formula(paste("V1~0|", paste(factors, collapse="+")))
                            reglog_5 <- mlogit(fmla, data=NCD, reflevel="1")
                        }









                    } else if (regr == "lm") {
                        ## Case of LINEAR REGRESSION MODEL

                        # Since we are performing a linear regression, each
                        # element of CD are numbers and have then to remain as
                        # class "numeric" (we don't have to perform some class
                        # adjustment as for the case of the creation of a
                        # multinomial model).

                        # Computation of the linear regression model
                        if(totVt==1){
                            reglog_5 <- lm(V1~0, data=CD)
                            # first case is evaluated aside
                        }

                        if(totVt>1){
                            # creation of "V2" ... "VtotVt" (to use them in the
                            # formula)
                            factors_character <- paste("V", 2:totVt, sep = "")
                            # Transformation of this object from character to
                            # vector (in order to be able to access its
                            # components)
                            factors <- as.vector(factors_character)
                            # Creation of a specific formula according to the
                            # value of totVt
                            fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                            reglog_5 <- lm(fmla, data=CD)
                        }











                    } else { # meaning (regr == "lrm")
                        ## Case of ORDINAL REGRESSION MODEL

                        # Linking to the package rms to use the function "lrm"

                        # Since we are performing an ordinal regression, each
                        # element of CD are numbers and have then to remain as
                        # class "numeric" (we don't have to perform some class
                        # adjustment as for the case of the creation of a
                        # multinomial model).

                        # Computation of the ordinal model
                        if(totVt==1){
                            reglog_5 <- lrm(V1~0, data=CD)
                            # first case is evaluated aside
                        }

                        if(totVt>1){
                            # creation of "V2" ... "VtotVt" (to use them in the
                            # formula)
                            factors_character <- paste("V", 2:totVt, sep = "")
                            # Transformation of this object from character to
                            # vector (in order to be able
                            # to access its components)
                            factors <- as.vector(factors_character)
                            # Creation of a specific formula according to the
                            # value of totVt
                            fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                            reglog_5 <- lrm(fmla, data=CD)
                        }

                    }


                    #*********************************************************************************
                    #*********************************************************************************



                    # 5.3.3 Imputation using the just created model (Dealing with the actual VALUES to impute) ----------------------------------------------------

                    # Conversion of ODi from data.frame to matrix
                    ODi <- as.matrix(ODi)

                    # Only PAST VIs are useful
                    for (order in 1:MaxTermGapSize){

                        # Analysing the value of parameter available
                        if (available==TRUE){
                            # we take the previously
                            # imputed data into account
                            LOOKUP <- ODi
                        } else { # that is available == FALSE and thus we
                            # don't take the previously imputed
                            # data into account
                            LOOKUP <- OD
                        }


                        # Assigning the current "REFORDT_order" matrix to the
                        # variable matrix REFORDT
                        # (according to the current value of "order")
                        tempObject = get(paste0("REFORDT_",order))
                        REFORDT <- as.matrix(tempObject)
                        if (ncol(REFORDT) == 1) {
                            REFORDT <- t(REFORDT)
                        }
                        nr_REFORDT <- nrow(REFORDT)

                        for (u in 1:nr_REFORDT) {
                            i <- REFORDT[u,1] # taking out the first coordinate
                            # (row number in ORDER) from REFORDT
                            j <- REFORDT[u,2] # taking out the second coordinate
                            # (column number in ORDER) from REFORDT

                            CDi <- matrix(NA,nrow=k,ncol=1)

                            # Matrix for past values
                            vect <- LOOKUP[i,(j-npt):(j-1)]
                            CDpi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)

                            # Matrix for past distribution
                            if (pastDistrib) {
                                dbi <- summary(factor(LOOKUP[i,1:(j-1)], levels=c(1:k), exclude=NULL))/length(1:(j-1))
                                CDdbi <- matrix(dbi[1:k], nrow=k, ncol=k, byrow=TRUE)
                            }

                            # Concatenating CDi
                            CDi <- cbind(CDi, CDpi)

                            if (pastDistrib) {
                                CDi <- cbind(CDi, CDdbi)
                            }

                            # Conversion of CDi into a data frame
                            CDi <- as.data.frame(CDi)

                            # Type transformation of the columns of CDi
                            # The first values of CDi must be of type factor
                            # (categorical values)
                            CDi[,(1:(1+npt))] <- lapply(CDi[,(1:(1+npt))],factor, levels=c(1:k,NA), exclude=NULL)
                            # The last values of CDi must be of type numeric
                            # (distributions)
                            if (pastDistrib) {
                                CDi[,(1+npt+1):ncol(CDi)] <- lapply(CDi[,(1+npt+1):ncol(CDi)],as.numeric)
                            }
                            # Eventually concatenating CDi with COi (the matrix
                            # containing the covariates)
                            if (all(is.na(CO))==FALSE) { # Checking if CO is NOT
                                # completely empty
                                # Creation of the matrix COi used in 3.3
                                COi <- do.call(rbind, replicate(k, as.data.frame(CO[i,]), simplify=FALSE))
                                # Concatenating CDi and COi into CDi
                                CDi <- cbind(CDi, COi)
                                # Transformation of the names of the columns of
                                # CDi (called V1, V2, ..., "VtotV")
                                colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                            }
                            # Else, in case CO is empty (i.e. we don't consider
                            # any covariate)
                            # simply continue with the current CDi

                            # Eventually concatenating CDi with COtselected_i
                            # (the matrix containing the current time-dependent
                            # covariates)
                            # Checking if COt is NOT completely empty
                            if (ncot > 0) {
                                COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
                                for (d in 1:(ncot/nc)) {
                                    COtselected_i <- cbind(COtselected_i, COt[i,(j) + (d-1)*nc])
                                }
                                COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
                                # Concatenating CDi and COtselected_i into CDi
                                CDi <- cbind(CDi, COtselected_i)
                                # Transformation of the names of the columns of
                                # CDi (called V1, V2, ..., "VtotV")
                                colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                            }


                            # Check for missing-values among predictors
                            if (max(is.na(CDi[1,2:totVt]))==0){


                                #*******************************************************************
                                # Initially "imputeValue.R"
                                #*******************************************************************

                                if (regr == "mlogit") {
                                    ## Case of MULTINOMIAL REGRESSION MODEL


                                    pred <- predict(reglog_5,newdata=CDi)
                                    # Example of value returned by pred:
                                    # (Sytematically, the upper line represents
                                    # the possible categories of the variable
                                    # (here, k=2, so the possible categories are
                                    # either 1" or "2"))
                                    #            1            2
                                    # 1.000000e+00 2.111739e-22
                                    #
                                    # Cumulative pred
                                    pred <- cumsum(pred)
                                    # Corresponding example value returned by
                                    # in the "cumulative pred":
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
                                    sel <- which(pred>=alea)
                                    # Corresponding example value returned in
                                    # sel:
                                    # 1 2
                                    # 1 2
                                    #















                                } else if (regr == "lm") {
                                    ## Case of LINEAR REGRESSION MODEL

                                    # Since we are performing a linear
                                    # regression, each element of CDi are
                                    # numbers and have to be considered as class
                                    # "numeric"
                                    CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                    pred <- predict(reglog_5, CDi)
                                    # Introducing a variance "noise"
                                    pred <- rnorm(length(pred),pred,noise)
                                    # Rounding pred to the nearest integer
                                    pred <- round(pred)
                                    # Restricting pred to its lowest value: 1
                                    pred <- ifelse(pred<1,1,pred)
                                    # Restricting pred to its highest value: k
                                    pred <- ifelse(pred>k,k,pred)
                                    sel <- pred













                                } else { # meaning (regr == "lrm")
                                    ## Case of ORDINAL REGRESSION MODEL

                                    # Since we are performing an ordinal
                                    # regression, each element of CDi are
                                    # numbers and have to
                                    # be considered as class "numeric"
                                    CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                    pred <- predict(reglog_5, CDi, type="fitted.ind")
                                    # Testing if we are in case where k=2 (if
                                    # this is the case, we need to create
                                    # the second complementary probility by hand
                                    # since lrm returns only the first
                                    # probability)
                                    if (k == 2) {
                                        pred <- c(pred,(1-pred))
                                    }
                                    # Cumulative pred
                                    pred <- cumsum(pred)
                                    # Introducing a variance "noise"
                                    pred <- rnorm(length(pred),pred,noise)
                                    # Checking that the components of vector
                                    # "pred" are still included in the interval
                                    # [0,1]
                                    pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                    pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                    # Imputation
                                    # Since we have introduce a noise on the
                                    # variance, it might occur that "alea"
                                    # is greater than the greatest value of
                                    # "pred". We have then to restrict "alea"
                                    # to the last value of "pred"
                                    alea <- runif(1)
                                    if (alea > pred[length(pred)]) {
                                        alea <- pred[length(pred)]
                                    }
                                    sel <- which(pred>=alea)



                                }

                                #*******************************************************************
                                #*******************************************************************


                                ODi[i,j] <- sel[1]
                            }
                        }
                    }
                }

























                # 6. Imputing SLG NAs -----------------------------------------------------------------------------------------------------------------------------

                if (max(ORDERSLGLeft)!=0) {     # Checking if we have to impute
                    # left-hand side SLG

                    warning("/!\\ Specially Located Gaps (SLG) have been detected on the left-hand side of OD.","\n","For certain missing data groups close to the border of OD, np may have been automatically reduced.","\n","If you don't want this to happen, please choose a lower value for np.")

                    # 6.2.LEFT Computation of the order of imputation of each MD ----------------------------------------------------------------------------------

                    # Creation of the temporary SLG matrices for the left-hand
                    # side of OD
                    for (h in 2:np) {
                        if (max(ORDERSLGLeft[,h])>0) {
                            # Checking if a gap begins
                            # somewhere on the current column
                            # If that is not the case, there is
                            # no need to perform
                            # the entire following procedure
                            # and we simply can go
                            # to the next column of ORDERSLGLeft

                            # initialization of a new temporary
                            # ORDERSLGLeft_temp matrix
                            ORDERSLGLeft_temp <- matrix(0,nrow=nr,ncol=nc)

                            for (i in 1:nr) { # going from top to bottom
                                j <- h        # storing the current column
                                # we are checking (and
                                # reinitializing j)

                                if (ORDERSLGLeft[i,j]>0 & ORDERSLGLeft[i,j-1]==0) {

                                    while (ORDERSLGLeft[i,j]>0) {
                                        # going from left to right
                                        ORDERSLGLeft_temp[i,j] <- ORDERSLGLeft[i,j]
                                        j <- j+1
                                    }
                                }
                            }

                            if (max(ORDERSLGLeft_temp)==0) {
                                next
                            }


                            # Update of np and totV
                            np_temp <- h-1
                            totV_temp <- 1+np_temp+nf+nco+(ncot/nc)
                            if (pastDistrib) {
                                totV_temp <- totV_temp + k
                            }
                            if (futureDistrib) {
                                totV_temp <- totV_temp + k
                            }





                            # (i.e. updating matrix ORDERSLGLeft_temp with
                            # coherent value of "order"
                            # (going from 1 to MaxGapSLGLeft))

                            # Adjusting the matrix ORDERSLGLeft_temp and
                            # storing the coordinates of the NA to impute

                            # In a similar manner to part 2.4., we go here one
                            # single time through the reduced version
                            # ORDERSLGLeft_temp of ORDERSLG and we create
                            # MaxGapSLGLeft "REFORDSLGRLeft_" matrices
                            # collecting the coordinates of each corresponding
                            # values in ORDERSLGLeft_temp which are greater
                            # than 0


                            # REFORDSLGLeft matrices
                            # Initialization of the REFORDSLGLeft matrices
                            MaxGapSLGLeft <- max(ORDERSLGLeft_temp[ORDERSLGLeft_temp!=0])-(min(ORDERSLGLeft_temp[ORDERSLGLeft_temp!=0]) - 1)

                            for(order in 1:MaxGapSLGLeft) {
                                assign(paste("REFORDSLGLeft_",order,sep=''),matrix(nrow=0,ncol=2))
                            }

                            ORDERSLGLeft_temp_init <- ORDERSLGLeft_temp

                            for (i in 1:nr) {
                                for (j in 1:nc) {
                                    if (ORDERSLGLeft_temp[i,j] > 0) {

                                        # Updating ORDERSLGLeft_temp so that it
                                        # becomes a matrix with positive values
                                        # going from 1 to MaxGapSLGLeft
                                        ORDERSLGLeft_temp[i,j] <- ORDERSLGLeft_temp[i,j] - (min(ORDERSLGLeft_temp_init[ORDERSLGLeft_temp_init!=0]) - 1)

                                        # Collecting the coordinates for each
                                        # value of "order"
                                        coord <- t(matrix(c(i,j)))

                                        tempObject = get(paste0("REFORDSLGLeft_",ORDERSLGLeft_temp[coord]))
                                        update <- rbind(tempObject,coord)

                                        assign( (paste("REFORDSLGLeft_",ORDERSLGLeft_temp[coord],sep='')),update )
                                    }
                                }
                            }









                            # 6.3.LEFT Imputation of the missing data listed by ORDERSLGLeft_temp and ORDERSLGRight_temp using a specific model -------------------


                            # 6.3.1.LEFT Building of the different data matrices CD -------------------------------------------------------------------------------
                            # for the computation of the model for every SLG
                            # on the left-hand side of OD

                            for (order in 1:MaxGapSLGLeft) { # /!\ "order"
                                # corresponds to the values
                                # of the components of
                                # ORDERSLGLeft_temp
                                # (i.e. the number of the iteration, the order
                                # in which the values are going to be imputed)


                                # Building of a data matrix for the computation
                                # of the model
                                ud <- nc-(MaxGapSLGLeft-order+np_temp+nf)
                                # number of usable data for each row of OD
                                # size of the current mobile caracteristic frame
                                # (that changes according to "order") which is
                                # equal to nc-ud+1
                                frameSize <- MaxGapSLGLeft-order+np_temp+nf+1
                                # Structure and building of the data matrix CD
                                # The first column of CD is the dependent
                                # variable (VD, response variable)
                                # The following columns are the independent
                                # variables (VIs, explanatory variables) coming
                                # from the past (np_temp>0) or the future (nf>0)
                                # ordered by time and the distribution of
                                # the possible values (i.e. all
                                # possible categorical variables numbered from
                                # 1 to k and of course the value NA)
                                # respectively Before and After the NA to impute
                                #
                                #           VD   Past VIs   Future VIS    Past distribution    Future distribution
                                #
                                #        /                                                                           \
                                #       |                                                                            |
                                # CD =  |  CD      CDp        CPf               CDdb                  CDda          |
                                #       \                                                                          /
                                #        \                                                                       /
                                #
                                # We are then going to create a specific CD
                                # according to the value of np_temp and nf

                                # initialization of the current very left part
                                # of the predictive model matrix ("matrice de
                                # modele de prediction") with NA everywhere
                                # (/!\ for each "order", we are going to build
                                # such a CD)
                                CD <- matrix(NA,nrow=nr*ud,ncol=1)

                                # Dealing with the change of shape of the
                                # prediction frame (according to whether the
                                # imputed data is located at the beginning
                                # (left) of a gap or at the end (right))
                                # The purpose of this if statement is to detect
                                # if we have to insert a shift (to jump at the
                                # end of the gap) or not
                                if ( (np_temp > 0 & nf > 0) & ( (MaxGapSLGLeft%%2==0 & order%%2==0) | (MaxGapSLGLeft%%2!=0 & order%%2!=0) )){
                                    shift <- MaxGapSLGLeft - order
                                    # jumping at the end of the gap
                                } else {
                                    shift <- 0
                                    # no shift is needed (note that no shift
                                    # is needed for the case of model 2
                                    # (only past) and model 3 (only future))
                                }

                                iter <- 1
                                # initialisation of the number of
                                # iterations of the following for loops


                                # For left SLG, naturally only the cases "only
                                # PAST" and "PAST and FUTURE" are possible to
                                # meet (because np_temp has to be greater than
                                # 0, otherwise, it would mean that we are not
                                # in the case of a SLG and that we
                                # can impute it as a usual internal gap)

                                # Only PAST
                                if (np_temp>0 & nf==0) {
                                    # only PAST VIs do exist
                                    # initialisation of matrix CDp
                                    CDp <- matrix(NA, nrow=nr*ud, ncol=np_temp)

                                    if (ncot > 0) {
                                        # initialisation of matrix COtselected
                                        COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
                                    }

                                    # initialisation of matrix CDdb (for past
                                    # distribution analysis)
                                    # (Distribution Before)
                                    if (pastDistrib) {
                                        CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
                                        # submatrix of CDdb: CDdb is composed
                                        # of ud matrix db on top of each other
                                        db <- matrix(NA, nrow=nr, ncol=k)
                                    }

                                    # initialisation of matrix CDda (for future
                                    # distribution analysis)
                                    # (Distribution After)
                                    if (futureDistrib) {
                                        # CDda has same dimensions as CDdb
                                        CDda <- matrix(NA, nrow=nr*ud, ncol=k)
                                        # da has same dimensions as db
                                        da <- matrix(NA, nrow=nr, ncol=k)
                                    }

                                    for (j in frameSize:nc){
                                        # /!\ j is initialised at
                                        # the very end (utmost
                                        # right) of the frame
                                        t1 <- (nr*(iter-1)+1) # determining the
                                        # locations of the time
                                        # span (always nr) of the
                                        # piled up blocks of CD
                                        t2 <- nr*iter
                                        # VD
                                        CD[t1:t2,1] <- OD[,j-frameSize+np_temp+1]
                                        # /!\ current pointer on column is thus:
                                        # "j-frameSize+np_temp+1"

                                        # Past VIs
                                        CDp[t1:t2,] <- OD[,(j-frameSize+1):(j-frameSize+np_temp)]

                                        # Eventually considering time-dependent
                                        # covariates
                                        if (ncot > 0) {
                                            COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                                            for (d in 1:(ncot/nc)) {
                                                COttemp <- cbind(COttemp, COt[,(j-frameSize+np_temp+1) + (d-1)*nc])
                                            }
                                            COtselected[t1:t2,] <- COttemp
                                        }

                                        # Past distribution (i.e. Before)
                                        if (pastDistrib) {
                                            ODt <- t(OD)
                                            ODt <- as.data.frame(ODt)
                                            tempOD <- lapply(ODt[(1:(j-frameSize+np_temp)),], factor, levels=c(1:k,NA), exclude=NULL)
                                            # because:
                                            # j-frameSize+np_temp+1 - 1 =
                                            # j-frameSize+np_temp

                                            db_list <- lapply(tempOD,summary)
                                            db_matrix <- do.call(rbind,db_list)
                                            CDdb[t1:t2,] <- db_matrix[,1:k]/length(1:(j-frameSize+np_temp))
                                        }

                                        # Future distribution (i.e. After)
                                        if (futureDistrib) {
                                            if ( (j-frameSize+np_temp+2) <= nc ) {
                                                ODt <- t(OD)
                                                ODt <- as.data.frame(ODt)
                                                tempOD <- lapply(ODt[((j-frameSize+np_temp+2):nc),], factor, levels=c(1:k,NA), exclude=NULL)
                                                # because:
                                                # j-frameSize+np_temp+1 + 1 =
                                                # j-frameSize+np_temp+2

                                                da_list <- lapply(tempOD,summary)
                                                da_matrix <- do.call(rbind,da_list)
                                                CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np_temp+2):nc)
                                            } else {
                                                # if index in OD exceeds OD
                                                # number of columns, the future
                                                # distribution of the k
                                                # categorical variables is
                                                # simply null
                                                CDda[t1:t2,] <- matrix(nrow=nr,ncol=k,0)
                                            }
                                        }

                                        iter <- iter+1
                                    }

                                    # Concatenating CD
                                    CD <- cbind(CD, CDp)

                                    if (pastDistrib) {
                                        CD <- cbind(CD, CDdb)
                                    }

                                    if (futureDistrib) {
                                        CD <- cbind(CD, CDda)
                                    }

                                    # Conversion of CD into a data frame
                                    CD <- as.data.frame(CD)

                                    # Eventually concatenating CD with COs
                                    # (the matrix containing the covariates)
                                    if (all(is.na(CO))==FALSE) {
                                        # Checking if CO is NOT
                                        # completely empty
                                        # Creation of the stacked covariates
                                        # matrix for 3.1
                                        COs <- do.call("rbind", rep(list(CO), ud))
                                        # Concatenating CD and COs into CD
                                        CD <- cbind(CD, COs)
                                    }
                                    # Else, in case CO is empty (i.e. we don't
                                    # consider any covariate)
                                    # simply continue with the current CD

                                    # Eventually concatenating CD with
                                    # COtselected (the matrix containing the
                                    # current time-dependent covariates)
                                    # Checking if COt is NOT completely empty
                                    if (ncot > 0) {
                                        # Concatenating CD and COtselected into
                                        # CD
                                        CD <- cbind(CD, as.data.frame(COtselected))
                                    }


                                    # PAST and FUTURE
                                } else {
                                    # meaning np_temp>0 and nf>0 and that,
                                    # thus, PAST as well as FUTURE VIs do exist
                                    # initialisation of matrices CDp and CDf
                                    CDp <- matrix(NA, nrow=nr*ud, ncol=np_temp)
                                    CDf <- matrix(NA, nrow=nr*ud, ncol=nf)

                                    if (ncot > 0) {
                                        # initialisation of matrix COtselected
                                        COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
                                    }

                                    # initialisation of matrix CDdb (for past
                                    # distribution analysis)
                                    # (Distribution Before)
                                    if (pastDistrib) {
                                        CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
                                        db <- matrix(NA, nrow=nr, ncol=k)
                                        # submatrix of CDdb: CDdb is composed
                                        # of ud matrix db on top of each other
                                    }

                                    # initialisation of matrix CDda (for future
                                    # distribution analysis)
                                    # (Distribution After)
                                    if (futureDistrib) {
                                        CDda <- matrix(NA, nrow=nr*ud, ncol=k)
                                        # CDda has same dimensions as CDdb
                                        da <- matrix(NA, nrow=nr, ncol=k)
                                        # da has same dimensions as db
                                    }


                                    for (j in frameSize:nc){
                                        t1 <- (nr*(iter-1)+1)
                                        t2 <- nr*iter
                                        # VD
                                        CD[t1:t2,1] <- OD[,j-frameSize+np_temp+1+shift]
                                        # /!\ current pointer on column is thus:
                                        # "j-frameSize+np_temp+1+shift"

                                        # Past VIs
                                        CDp[t1:t2,] <- OD[,(j-frameSize+1):(j-frameSize+np_temp)]
                                        # Future VIs
                                        CDf[t1:t2,] <- OD[,(j-nf+1):j]

                                        # Eventually considering time-dependent
                                        # covariates
                                        if (ncot > 0) {
                                            COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                                            for (d in 1:(ncot/nc)) {
                                                COttemp <- cbind(COttemp, COt[,(j-frameSize+np_temp+1+shift) + (d-1)*nc])
                                            }
                                            COtselected[t1:t2,] <- COttemp
                                        }

                                        # Past distribution (i.e. Before)
                                        if (pastDistrib) {
                                            ODt <- t(OD)
                                            ODt <- as.data.frame(ODt)
                                            tempOD <- lapply(ODt[(1:(j-frameSize+np_temp+shift)),], factor, levels=c(1:k,NA), exclude=NULL)
                                            # because:
                                            # j-frameSize+np_temp+1+shift - 1 =
                                            # j-frameSize+np_temp+shift

                                            db_list <- lapply(tempOD,summary)
                                            db_matrix <- do.call(rbind,db_list)
                                            CDdb[t1:t2,] <- db_matrix[,1:k]/length(1:(j-frameSize+np_temp+shift))
                                        }

                                        # Future distribution (i.e. After)
                                        if (futureDistrib) {
                                            ODt <- t(OD)
                                            ODt <- as.data.frame(ODt)
                                            tempOD <- lapply(ODt[((j-frameSize+np_temp+shift+2):nc),], factor, levels=c(1:k,NA), exclude=NULL)
                                            # because:
                                            # j-frameSize+np_temp+1+shift + 1 =
                                            # j-frameSize+np_temp+shift+2
                                            da_list <- lapply(tempOD,summary)
                                            da_matrix <- do.call(rbind,da_list)
                                            CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np_temp+shift+2):nc)
                                        }

                                        iter <- iter+1
                                    }

                                    # Concatenating CD
                                    CD <- cbind(CD, CDp, CDf)

                                    if (pastDistrib) {
                                        CD <- cbind(CD, CDdb)
                                    }

                                    if (futureDistrib) {
                                        CD <- cbind(CD, CDda)
                                    }

                                    # Conversion of CD into a data frame
                                    CD <- as.data.frame(CD)

                                    # Eventually concatenating CD with COs
                                    # (the matrix containing the covariates)
                                    if (all(is.na(CO))==FALSE) {
                                        # Checking if CO is NOT
                                        # completely empty
                                        # Creation of the stacked
                                        # covariates matrix for 3.1
                                        COs <- do.call("rbind", rep(list(CO), ud))
                                        # Concatenating CD and COs into CD
                                        CD <- cbind(CD, COs)
                                    }
                                    # Else, in case CO is empty (i.e. we don't
                                    # consider any covariate)
                                    # simply continue with the current CD

                                    # Eventually concatenating CD with
                                    # COtselected (the matrix containing the
                                    # current time-dependent covariates)
                                    # Checking if COt is NOT completely empty
                                    if (ncot > 0) {
                                        # Concatenating CD and COtselected into
                                        # CD
                                        CD <- cbind(CD, as.data.frame(COtselected))
                                    }


                                }







                                # 6.3.2.LEFT Computation of the model (Dealing with the LOCATIONS of imputation) --------------------------------------------------



                                #*********************************************************************************
                                # Initially "computeModel.R"
                                #*********************************************************************************




                                # ==>> Manipulation of parameters

                                # Conversion of CD in a data frame
                                CD <- as.data.frame(CD)

                                # Transformation of the names of the columns
                                # of CD (called V1, V2, ..., "VtotV_temp")
                                colnames(CD) <- paste("V", 1:ncol(CD), sep = "")










                                if (regr == "mlogit") {
                                    ## Case of MULTINOMIAL REGRESSION MODEL

                                    # Linking the package mlogit

                                    # By default, every column of CD are of
                                    # class "numeric". Thus, there is no need to
                                    # convert the columns containing the
                                    # distribution data to class "numeric".
                                    # Moreover the class of the covariates
                                    # columns at the very end are ALREADY
                                    # set correctly and we don't need to
                                    # update them.
                                    # On the other hand, the first columns of CD
                                    # (1 up to 1+np+nf) have to be of
                                    # class "factor" (because they are
                                    # the columns containing our categorical
                                    # data coming from OD).

                                    # Transformation of the first columns (i.e.
                                    # the categorical values) of CD (column 1
                                    # up to column 1+np+nf) into factor
                                    CD[,(1:(1+np_temp+nf))] <- lapply(CD[,(1:(1+np_temp+nf))],factor, levels=c(1:k,NA), exclude=NULL)


                                    # Dataframe for mlogit
                                    NCD <- mlogit.data(CD, varying=NULL, choice="V1", shape="wide")


                                    # Computation of the multinomial model
                                    if(totV_temp==1){
                                        # First case is evaluated aside
                                        reglog_6_left <- mlogit(V1~0, data=NCD, reflevel="1")
                                    }

                                    if(totV_temp>1){
                                        # creation of "V2" ... "VtotV_temp"
                                        # (to use them in the formula)
                                        factors_character <- paste("V", 2:totV_temp, sep = "")
                                        # Transformation of this object
                                        # from character to vector (in order to
                                        # be able to access its components)
                                        factors <- as.vector(factors_character)
                                        # Creation of a specific formula
                                        # according to the value of totV_temp
                                        fmla <- as.formula(paste("V1~0|", paste(factors, collapse="+")))
                                        reglog_6_left <- mlogit(fmla, data=NCD, reflevel="1")
                                    }









                                } else if (regr == "lm") {
                                    ## Case of LINEAR REGRESSION MODEL

                                    # Since we are performing a linear
                                    # regression, each element of CD are numbers
                                    # and have then to remain as class
                                    # "numeric" (we don't have to perform
                                    # some class adjustment as for the case
                                    # of the creation of a multinomial model).

                                    # Computation of the linear regression model
                                    if(totV_temp==1){
                                        reglog_6_left <- lm(V1~0, data=CD)
                                        # first case is evaluated aside
                                    }

                                    if(totV_temp>1){
                                        # creation of "V2" ... "VtotV_temp"
                                        # (to use them in the formula)
                                        factors_character <- paste("V", 2:totV_temp, sep = "")
                                        # Transformation of this object from
                                        # character to vector (in order to be
                                        # able to access its components)
                                        factors <- as.vector(factors_character)
                                        # Creation of a specific formula
                                        # according to the value of totV_temp
                                        fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                                        reglog_6_left <- lm(fmla, data=CD)
                                    }











                                } else { # meaning (regr == "lrm")
                                    ## Case of ORDINAL REGRESSION MODEL

                                    # Linking to the package rms to use the
                                    # function "lrm"

                                    # Since we are performing an ordinal
                                    # regression, each element of CD are numbers
                                    # and have then to remain as class "numeric"
                                    # (we don't have to perform some
                                    # class adjustment as for the case
                                    # of the creation of a multinomial model).

                                    # Computation of the ordinal model
                                    if(totV_temp==1){
                                        reglog_6_left <- lrm(V1~0, data=CD)
                                        # first case is evaluated aside
                                    }

                                    if(totV_temp>1){
                                        # creation of "V2" ... "VtotV_temp"
                                        # (to use them in the formula)
                                        factors_character <- paste("V", 2:totV_temp, sep = "")
                                        # Transformation of this object from
                                        # character to vector (in order to be
                                        # able to access its components)
                                        factors <- as.vector(factors_character)
                                        # Creation of a specific formula
                                        # according to the value of totV_temp
                                        fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                                        reglog_6_left <- lrm(fmla, data=CD)
                                    }

                                }


                                #*********************************************************************************
                                #*********************************************************************************




                                # 6.3.3.LEFT Imputation using the just created model (Dealing with the actual VALUES to impute) -----------------------------------

                                # Structure and building of the data matrix CDi
                                # The first column of CDi is the dependent
                                # variable (VD, response variable) that we have
                                # to implement during the current iteration
                                # (i.e. automatically an NA).
                                # The following columns are the corresponding
                                # independent variables
                                # (VIs, explanatory variables)
                                # coming from the past (np_temp>0) or the future
                                # (nf>0) (ordered by time and corresponding to
                                # the current predictive pattern) and the
                                # distribution of the possible values (i.e. all
                                # possible categorical variables numbered from 1
                                # to k and of course the value NA) respectively
                                # Before and After the NA to impute
                                # (The fact that every lines of CDi
                                # are identical is related to the working of the
                                # function mlogit that has to have as much lines
                                # in CDi as there are categories of the variable
                                # (see the parameter "k")
                                # --> so, CDi is composed of k identical lines)
                                #
                                #            VD   Past VIs   Future VIS    Past distribution     Future distribution
                                #
                                #         /                                                                             \
                                #        |                                                                              |
                                # CDi =  |  CDi    CDpi       CPfi              CDdbi                   CDdai           |
                                #        \                                                                             /
                                #         \                                                                          /
                                #
                                # We are then going to create a specific CDi
                                # according to the value of np_temp and nf





                                # Analysing the value of parameter available
                                if (available==TRUE){
                                    # we take the previously imputed data
                                    # into account
                                    LOOKUP <- ODi
                                } else {
                                    # that is available == FALSE and thus we
                                    # don't take the previously imputed data
                                    # into account
                                    LOOKUP <- OD
                                }



                                # Assigning the current "REFORDSLGLeft_order"
                                # matrix to the variable matrix REFORDSLGLeft
                                # (according to the current value of "order")
                                tempObject = get(paste0("REFORDSLGLeft_",order))
                                REFORDSLGLeft <- as.matrix(tempObject)
                                if (ncol(REFORDSLGLeft) == 1) {
                                    REFORDSLGLeft <- t(REFORDSLGLeft)
                                }
                                nr_REFORD <- nrow(REFORDSLGLeft)



                                if (np_temp>0 & nf==0){
                                    # only PAST VIs do exist
                                    for (u in 1:nr_REFORD) {
                                        q <- REFORDSLGLeft[u,1]
                                        # taking out the first
                                        # coordinate (row number
                                        # in ORDER) from
                                        # REFORDSLGLeft
                                        w <- REFORDSLGLeft[u,2]
                                        # taking out the second
                                        # coordinate (column
                                        # number in ORDER) from
                                        # REFORDSLGLeft

                                        CDi <- matrix(NA, nrow=k, ncol=1)

                                        # Matrix for past values
                                        vect <- LOOKUP[q,(w-np_temp):(w-1)]
                                        # /!\ current pointer on column is
                                        # thus: "w"

                                        CDpi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)

                                        # Matrix for past distribution
                                        if (pastDistrib) {
                                            dbi <- summary(factor(LOOKUP[q,1:(w-1)], levels=c(1:k), exclude=NULL))/length(1:(w-1))
                                            CDdbi <- matrix(dbi[1:k], nrow=k, ncol=k, byrow=TRUE)
                                        }

                                        # Matrix for future distribution
                                        if (futureDistrib) {
                                            dai <- summary(factor(LOOKUP[q,(w+1):nc], levels=c(1:k), exclude=NULL))/length((w+1):nc)
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

                                        # Type transformation of the columns of
                                        # CDi The first values of CDi must be
                                        # of type factor (categorical values)
                                        CDi[,(1:(1+np_temp+nf))] <- lapply(CDi[,(1:(1+np_temp+nf))],factor, levels=c(1:k,NA), exclude=NULL)
                                        # The last values of CDi must be of type
                                        # numeric (distributions)
                                        if (pastDistrib | futureDistrib) {
                                            CDi[,(1+np_temp+nf+1):ncol(CDi)] <- lapply(CDi[,(1+np_temp+nf+1):ncol(CDi)],as.numeric)
                                        }
                                        # Eventually concatenating CDi with COi
                                        # (the matrix containing the covariates)
                                        if (all(is.na(CO))==FALSE) {
                                            # Checking if CO is
                                            # NOT completely
                                            # empty
                                            # Creation of the matrix COi used
                                            # in 3.3
                                            COi <- do.call(rbind, replicate(k, as.data.frame(CO[i,]), simplify=FALSE))
                                            # Concatenating CDi and COi into CDi
                                            CDi <- cbind(CDi, COi)
                                            # Transformation of the names of the
                                            # columns of CDi
                                            # (called V1, V2, ..., "VtotV")
                                            colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                                        }
                                        # Else, in case CO is empty (i.e. we
                                        # don't consider any covariate)
                                        # simply continue with the current CDi

                                        # Eventually concatenating CDi with
                                        # COtselected_i (the matrix containing
                                        # the current time-dependent covariates)
                                        # Checking if COt is NOT completely
                                        # empty
                                        if (ncot > 0) {
                                            COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
                                            for (d in 1:(ncot/nc)) {
                                                COtselected_i <- cbind(COtselected_i, COt[q,(w) + (d-1)*nc])
                                            }
                                            COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
                                            # Concatenating CDi and
                                            # COtselected_i into CDi
                                            CDi <- cbind(CDi, COtselected_i)
                                            # Transformation of the names of the
                                            # columns of CDi
                                            # (called V1, V2, ..., "VtotV")
                                            colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                                        }


                                        # Check for missing-values among
                                        # predictors
                                        if (max(is.na(CDi[1,2:totV_temp]))==0){


                                            #*******************************************************************
                                            # Initially "imputeValue.R"
                                            #*******************************************************************

                                            if (regr == "mlogit") {
                                                ## Case of MULTINOMIAL
                                                ## REGRESSION MODEL


                                                pred <- predict(reglog_6_left,newdata=CDi)
                                                # Example of value returned by
                                                # pred:
                                                # (Sytematically, the upper line
                                                # represents the
                                                # possible categories of
                                                # the variable
                                                # (here, k=2, so the possible
                                                # categories are either
                                                # 1" or "2"))
                                                #            1            2
                                                # 1.000000e+00 2.111739e-22
                                                #
                                                # Cumulative pred
                                                pred <- cumsum(pred)
                                                # Corresponding example value
                                                # returned by in the
                                                # "cumulative pred":
                                                # 1 2
                                                # 1 1
                                                #
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Checking that the components
                                                # of vector "pred" are
                                                # still included in the
                                                # interval [0,1]
                                                pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                                pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                                # Imputation
                                                alea <- runif(1)
                                                # Example value returned in
                                                # "alea":
                                                # [1] 0.2610005
                                                #
                                                sel <- which(pred>=alea)
                                                # Corresponding example value
                                                # returned in sel:
                                                # 1 2
                                                # 1 2
                                                #















                                            } else if (regr == "lm") {
                                                ## Case of LINEAR REGRESSION
                                                ## MODEL

                                                # Since we are performing a
                                                # linear regression, each
                                                # element of CDi are
                                                # numbers and have to be
                                                # considered as class "numeric"
                                                CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                                pred <- predict(reglog_6_left, CDi)
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Rounding pred to the nearest
                                                # integer
                                                pred <- round(pred)
                                                # Restricting pred to its
                                                # lowest value: 1
                                                pred <- ifelse(pred<1,1,pred)
                                                # Restricting pred to its
                                                # highest value: k
                                                pred <- ifelse(pred>k,k,pred)
                                                sel <- pred













                                            } else { # meaning (regr == "lrm")
                                                ## Case of ORDINAL REGRESSION
                                                ## MODEL

                                                # Since we are performing an
                                                # ordinal regression, each
                                                # element of CDi are numbers and
                                                # have to be considered as
                                                # class "numeric"
                                                CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                                pred <- predict(reglog_6_left, CDi, type="fitted.ind")
                                                # Testing if we are in case
                                                # where k=2 (if this is the
                                                # case, we need to create
                                                # the second complementary
                                                # probility by hand since lrm
                                                # returns only the first
                                                # probability)
                                                if (k == 2) {
                                                    pred <- c(pred,(1-pred))
                                                }
                                                # Cumulative pred
                                                pred <- cumsum(pred)
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Checking that the components
                                                # of vector "pred" are
                                                # still included in the
                                                # interval [0,1]
                                                pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                                pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                                # Imputation
                                                # Since we have introduce a
                                                # noise on the variance, it
                                                # might occur that "alea" is
                                                # greater than the greatest
                                                # value of "pred". We have then
                                                # to restrict "alea" to the last
                                                # value of "pred"
                                                alea <- runif(1)
                                                if (alea > pred[length(pred)]) {
                                                    alea <- pred[length(pred)]
                                                }
                                                sel <- which(pred>=alea)



                                            }

                                            #*******************************************************************
                                            #*******************************************************************


                                            ODi[q,w] <- sel[1]
                                        }

                                    }


                                } else {
                                    # meaning np_temp>0 and nf>0 and that, thus,
                                    # PAST as well as FUTURE VIs do exist
                                    for (u in 1:nr_REFORD) {
                                        q <- REFORDSLGLeft[u,1]
                                        # taking out the first
                                        # coordinate (row number
                                        # in ORDER) from
                                        # REFORDSLGLeft
                                        w <- REFORDSLGLeft[u,2]
                                        # taking out the second
                                        # coordinate (column
                                        # number in ORDER) from
                                        # REFORDSLGLeft

                                        CDi <- matrix(NA, nrow=k, ncol=1)

                                        # Matrix for past values
                                        vect <- LOOKUP[q,(w-shift-np_temp):(w-shift-1)]
                                        CDpi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)

                                        # Matrix for future values
                                        vect <- LOOKUP[q,(w-shift+MaxGapSLGLeft-order+1):(w-shift+MaxGapSLGLeft-order+nf)]
                                        CDfi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)

                                        # Matrix for past distribution
                                        if (pastDistrib) {
                                            dbi <- summary(factor(LOOKUP[q,1:(w-1)], levels=c(1:k), exclude=NULL))/length(1:(w-1))
                                            CDdbi <- matrix(dbi[1:k], nrow=k, ncol=k, byrow=TRUE)
                                        }

                                        # Matrix for future distribution
                                        if (futureDistrib) {
                                            dai <- summary(factor(LOOKUP[q,(w+1):nc], levels=c(1:k), exclude=NULL))/length((w+1):nc)
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

                                        # Type transformation of the columns of
                                        # CDi
                                        # The first values of CDi must be of
                                        # type factor (categorical values)
                                        CDi[,(1:(1+np_temp+nf))] <- lapply(CDi[,(1:(1+np_temp+nf))],factor, levels=c(1:k,NA), exclude=NULL)
                                        # The last values of CDi must be of type
                                        # numeric (distributions)
                                        if (pastDistrib | futureDistrib) {
                                            CDi[,(1+np_temp+nf+1):ncol(CDi)] <- lapply(CDi[,(1+np_temp+nf+1):ncol(CDi)],as.numeric)
                                        }
                                        # Eventually concatenating CDi with COi
                                        # (the matrix containing the covariates)
                                        if (all(is.na(CO))==FALSE) {
                                            # Checking if CO
                                            # is NOT completely
                                            # empty
                                            # Creation of the matrix COi used
                                            # in 3.3
                                            COi <- do.call(rbind, replicate(k, as.data.frame(CO[i,]), simplify=FALSE))
                                            # Concatenating CDi and COi into CDi
                                            CDi <- cbind(CDi,COi)
                                            # Transformation of the names of the
                                            # columns of CDi
                                            # (called V1, V2, ..., "VtotV")
                                            colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                                        }
                                        # Else, in case CO is empty (i.e. we
                                        # don't consider any covariate)
                                        # simply continue with the current CDi

                                        # Eventually concatenating CDi with
                                        # COtselected_i (the matrix containing
                                        # the current time-dependent covariates)
                                        # Checking if COt is NOT completely
                                        # empty
                                        if (ncot > 0) {
                                            COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
                                            for (d in 1:(ncot/nc)) {
                                                COtselected_i <- cbind(COtselected_i, COt[q,(w) + (d-1)*nc])
                                            }
                                            COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
                                            # Concatenating CDi and
                                            # COtselected_i into CDi
                                            CDi <- cbind(CDi, COtselected_i)
                                            # Transformation of the names of the
                                            # columns of CDi
                                            # (called V1, V2, ..., "VtotV")
                                            colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                                        }


                                        # Check for missing-values among
                                        # predictors (i.e. we won't impute any
                                        # value on the current MD if there is
                                        # any NA among the VIs)
                                        if (max(is.na(CDi[1,2:totV_temp]))==0){
                                            # checking that there is no NA among
                                            # the current VI (otherwise no data
                                            # will be imputed for the current
                                            # NA)


                                            #*******************************************************************
                                            # Initially "imputeValue.R"
                                            #*******************************************************************

                                            if (regr == "mlogit") {
                                                ## Case of MULTINOMIAL
                                                ## REGRESSION MODEL


                                                pred <- predict(reglog_6_left,newdata=CDi)
                                                # Example of value returned by
                                                # pred:
                                                # (Sytematically, the upper line
                                                # represents the
                                                # possible categories of
                                                # the variable
                                                # (here, k=2, so the possible
                                                # categories are either
                                                # 1" or "2"))
                                                #            1            2
                                                # 1.000000e+00 2.111739e-22
                                                #
                                                # Cumulative pred
                                                pred <- cumsum(pred)
                                                # Corresponding example value
                                                # returned by in the
                                                # "cumulative pred":
                                                # 1 2
                                                # 1 1
                                                #
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Checking that the components
                                                # of vector "pred" are
                                                # still included in the
                                                # interval [0,1]
                                                pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                                pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                                # Imputation
                                                alea <- runif(1)
                                                # Example value returned in
                                                # "alea":
                                                # [1] 0.2610005
                                                #
                                                sel <- which(pred>=alea)
                                                # Corresponding example
                                                # value returned in sel:
                                                # 1 2
                                                # 1 2
                                                #















                                            } else if (regr == "lm") {
                                                ## Case of LINEAR REGRESSION
                                                ## MODEL

                                                # Since we are performing a
                                                # linear regression, each
                                                # element of CDi
                                                # are numbers and have to be
                                                # considered as class "numeric"
                                                CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                                pred <- predict(reglog_6_left, CDi)
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Rounding pred to
                                                # the nearest integer
                                                pred <- round(pred)
                                                # Restricting pred to
                                                # its lowest value: 1
                                                pred <- ifelse(pred<1,1,pred)
                                                # Restricting pred to
                                                # its highest value: k
                                                pred <- ifelse(pred>k,k,pred)
                                                sel <- pred













                                            } else { # meaning (regr == "lrm")
                                                ## Case of ORDINAL REGRESSION
                                                ## MODEL

                                                # Since we are performing an
                                                # ordinal regression, each
                                                # element of CDi
                                                # are numbers and have to be
                                                # considered as class "numeric"
                                                CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                                pred <- predict(reglog_6_left, CDi, type="fitted.ind")
                                                # Testing if we are in case
                                                # where k=2 (if this is the
                                                # case, we need
                                                # to create the second
                                                # complementary probility by
                                                # hand since lrm returns only
                                                # the first probability)
                                                if (k == 2) {
                                                    pred <- c(pred,(1-pred))
                                                }
                                                # Cumulative pred
                                                pred <- cumsum(pred)
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Checking that the components
                                                # of vector "pred" are
                                                # still included in the
                                                # interval [0,1]
                                                pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                                pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                                # Imputation
                                                # Since we have introduce a
                                                # noise on the variance, it
                                                # might occur that "alea"
                                                # is greater than the greatest
                                                # value of "pred". We have
                                                # then to restrict "alea"
                                                # to the last value of "pred"
                                                alea <- runif(1)
                                                if (alea > pred[length(pred)]) {
                                                    alea <- pred[length(pred)]
                                                }
                                                sel <- which(pred>=alea)



                                            }

                                            #*******************************************************************
                                            #*******************************************************************


                                            ODi[q,w] <- sel[1]
                                        }
                                    }

                                }

                            }

                        }

                    }

                }











                if (max(ORDERSLGRight)!=0) {
                    # Checking if we have to impute right-hand
                    # side SLG

                    warning("/!\\ Specially Located Gaps (SLG) have been detected on the right-hand side of OD.","\n","For certain missing data groups close to the border of OD, nf may have been automatically reduced.","\n","If you don't want this to happen, please choose a lower value for nf.")

                    # 6.2.RIGHT Computation of the order of imputation of each MD ---------------------------------------------------------------------------------

                    # Creation of the temporary SLG matrices for the right-hand
                    # side of OD
                    for (h in (nc-1):(nc-nf+1)) {
                        if (max(ORDERSLGRight[,h])>0) {
                            # Checking if a gap begins
                            # somewhere on the current
                            # column.
                            # If that is not the case, there is no need to
                            # perform the entire following procedure and we
                            # simply can go to the next column of ORDERSLGRight

                            # initialization of a new temporary
                            # ORDERSLGRight_temp matrix
                            ORDERSLGRight_temp <- matrix(0,nrow=nr,ncol=nc)

                            for (i in 1:nr) {
                                j <- h      # storing the current column we are
                                # checking (and reinitializing j)

                                if (ORDERSLGRight[i,j]>0 & ORDERSLGRight[i,j+1]==0) {

                                    while (ORDERSLGRight[i,j]>0) {
                                        ORDERSLGRight_temp[i,j] <- ORDERSLGRight[i,j]
                                        j <- j-1
                                    }
                                }
                            }

                            if (max(ORDERSLGRight_temp)==0) {
                                next
                            }


                            # Update of nf and totV
                            nf_temp <- nc-h
                            totV_temp <- 1+np+nf_temp+nco+(ncot/nc)
                            if (pastDistrib) {
                                totV_temp <- totV_temp + k
                            }
                            if (futureDistrib) {
                                totV_temp <- totV_temp + k
                            }





                            # (i.e. updating matrix ORDERSLGRight_temp with
                            # coherent value of "order" (going from 1 to
                            # MaxGapSLGRight))

                            # Adjusting the matrix ORDERSLGRight_temp and
                            # storing the coordinates of the NA to impute

                            # In a similar manner to part 2.4., we go here one
                            # single time through the reduced version
                            # ORDERSLGRight_temp of ORDERSLG and we create
                            # MaxGapSLGLRight "REFORDSLGRight_" matrices
                            # collecting the coordinates of
                            # each corresponding values in
                            # ORDERSLGRight_temp which are
                            # greater than 0


                            # REFORDSLGRight matrices
                            # Initialization of the REFORDSLGRight matrices
                            MaxGapSLGRight <- max(ORDERSLGRight_temp[ORDERSLGRight_temp!=0])-(min(ORDERSLGRight_temp[ORDERSLGRight_temp!=0]) - 1)

                            for(order in 1:MaxGapSLGRight) {
                                assign(paste("REFORDSLGRight_",order,sep=''),matrix(nrow=0,ncol=2))
                            }

                            ORDERSLGRight_temp_init <- ORDERSLGRight_temp

                            for (i in 1:nr) {
                                for (j in 1:nc) {
                                    if (ORDERSLGRight_temp[i,j] > 0) {

                                        # Updating ORDERSLGRight_temp so that it
                                        # becomes a matrix with positive values
                                        # going from 1 to MaxGapSLGRight
                                        ORDERSLGRight_temp[i,j] <- ORDERSLGRight_temp[i,j] - (min(ORDERSLGRight_temp_init[ORDERSLGRight_temp_init!=0]) - 1)

                                        # Collecting the coordinates for each
                                        # value of "order"
                                        coord <- t(matrix(c(i,j)))

                                        tempObject = get(paste0("REFORDSLGRight_",ORDERSLGRight_temp[coord]))
                                        update <- rbind(tempObject,coord)

                                        assign( (paste("REFORDSLGRight_",ORDERSLGRight_temp[coord],sep='')),update )
                                    }
                                }
                            }









                            # 6.3.RIGHT Imputation of the missing data listed by ORDERSLGLeft_temp and ORDERSLGRight_temp using a specific model ------------------

                            # 6.3.1.RIGHT Building of the different data matrices CD ------------------------------------------------------------------------------
                            # for the computation of the model for every SLG
                            # on the right-hand side of OD

                            for (order in 1:MaxGapSLGRight) {
                                # /!\ "order" corresponds
                                # to the values of the
                                # components of
                                # ORDERSLGRight_temp
                                # (i.e. the number of the iteration,
                                # the order in which
                                # the values are going to be imputed)


                                # Building of a data matrix for the computation
                                # of the model

                                # number of usable data for each row of OD
                                ud <- nc-(MaxGapSLGRight-order+np+nf_temp)

                                # size of the current mobile caracteristic frame
                                # (that changes according to "order") which is
                                # equal to nc-ud+1
                                frameSize <- MaxGapSLGRight-order+np+nf_temp+1

                                # Structure and building of the data matrix CD
                                # The first column of CD is the dependent
                                # variable (VD, response variable)
                                # The following columns are the independent
                                # variables (VIs, explanatory variables) coming
                                # from the past (np>0) or the future (nf_temp>0)
                                # ordered by time and the distribution of the
                                # possible values (i.e. all possible categorical
                                # variables numbered from 1 to k and of course
                                # the value NA) respectively Before and After
                                # the NA to impute
                                #
                                #           VD   Past VIs   Future VIS    Past distribution    Future distribution
                                #
                                #        /                                                                           \
                                #       |                                                                            |
                                # CD =  |  CD      CDp        CPf               CDdb                  CDda          |
                                #       \                                                                          /
                                #        \                                                                       /
                                #
                                # We are then going to create a specific CD
                                # according to the value of np and nf_temp

                                # initialization of the current very left part
                                # of the predictive model matrix ("matrice de
                                # modele de prediction") with NA everywhere
                                # (/!\ for each "order", we are going to build
                                # such a CD)
                                CD <- matrix(NA,nrow=nr*ud,ncol=1)

                                # Dealing with the change of shape of the
                                # prediction frame (according to whether the
                                # imputed data is located at the beginning
                                # (left) of a gap or at the end (right))
                                # The purpose of this if statement is to detect
                                # if we have to insert a shift (to jump at the
                                # end of the gap) or not
                                if ( (np > 0 & nf_temp > 0) & ( (MaxGapSLGRight%%2==0 & order%%2==0) | (MaxGapSLGRight%%2!=0 & order%%2!=0) )){
                                    shift <- MaxGapSLGRight - order
                                    # jumping at the end
                                    # of the gap
                                } else {
                                    shift <- 0
                                    # no shift is needed (note that no shift
                                    # is needed for the case of model 2 (only
                                    # past) and model 3 (only future))
                                }

                                iter <- 1      # initialisation of the number of
                                # iterations of the following for loops


                                # For right SLG, naturally only the cases
                                # "only FUTURE" and "PAST and FUTURE" are
                                # possible to meet (because nf_temp has to be
                                # greater than 0, otherwise, it would mean that
                                # we are not in the case of a SLG and that we
                                # can impute it as a usual internal gap)

                                # Only FUTURE
                                if (np==0 & nf_temp>0) {
                                    # only FUTURE VIs do exist
                                    # initialisation of matrix CDf
                                    CDf <- matrix(NA, nrow=nr*ud, ncol=nf_temp)

                                    if (ncot > 0) {
                                        # initialisation of matrix COtselected
                                        COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
                                    }

                                    # initialisation of matrix CDdb (for past
                                    # distribution analysis)
                                    # (Distribution Before)
                                    if (pastDistrib) {
                                        CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
                                        db <- matrix(NA, nrow=nr, ncol=k)
                                        # submatrix of CDdb: CDdb is composed
                                        # of ud matrix db on top of each other
                                    }

                                    # initialisation of matrix CDda (for future
                                    # distribution analysis)
                                    # (Distribution After)
                                    if (futureDistrib) {
                                        # CDda has same dimensions as CDdb
                                        CDda <- matrix(NA, nrow=nr*ud, ncol=k)

                                        # da has same dimensions as db
                                        da <- matrix(NA, nrow=nr, ncol=k)
                                    }


                                    for (j in frameSize:nc){
                                        t1 <- (nr*(iter-1)+1)
                                        t2 <- nr*iter
                                        # VD
                                        CD[t1:t2,1] <- OD[,j-frameSize+np+1]
                                        # /!\ current pointer on column is thus:
                                        # "j-frameSize+np+1"

                                        # Future VIs
                                        CDf[t1:t2,] <- OD[,(j-nf_temp+1):j]

                                        # Eventually considering time-dependent
                                        # covariates
                                        if (ncot > 0) {
                                            COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                                            for (d in 1:(ncot/nc)) {
                                                COttemp <- cbind(COttemp, COt[,(j-frameSize+np+1) + (d-1)*nc])
                                            }
                                            COtselected[t1:t2,] <- COttemp
                                        }

                                        # Past distribution (i.e. Before)
                                        if (pastDistrib) {
                                            ODt <- t(OD)
                                            ODt <- as.data.frame(ODt)
                                            tempOD <- lapply(ODt[(1:(j-frameSize+np)),], factor, levels=c(1:k,NA), exclude=NULL)
                                            # because: j-frameSize+np+1 - 1 =
                                            # j-frameSize+np
                                            db_list <- lapply(tempOD,summary)
                                            db_matrix <- do.call(rbind,db_list)
                                            CDdb[t1:t2,] <- db_matrix[,1:k]/length(1:(j-frameSize+np))
                                        }

                                        # Future distribution (i.e. After)
                                        if (futureDistrib) {
                                            ODt <- t(OD)
                                            ODt <- as.data.frame(ODt)
                                            tempOD <- lapply(ODt[((j-frameSize+np+2):nc),], factor, levels=c(1:k,NA), exclude=NULL)
                                            # because: j-frameSize+np+1 + 1 =
                                            # j-frameSize+np+2
                                            da_list <- lapply(tempOD,summary)
                                            da_matrix <- do.call(rbind,da_list)
                                            CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np+2):nc)
                                        }

                                        iter <- iter+1
                                    }

                                    # Concatenating CD
                                    CD <- cbind(CD, CDf)

                                    if (pastDistrib) {
                                        CD <- cbind(CD, CDdb)
                                    }

                                    if (futureDistrib) {
                                        CD <- cbind(CD, CDda)
                                    }

                                    # Conversion of CD into a data frame
                                    CD <- as.data.frame(CD)

                                    # Eventually concatenating CD with COs
                                    # (the matrix containing the covariates)
                                    if (all(is.na(CO))==FALSE) {
                                        # Checking if CO is NOT
                                        # completely empty
                                        # Creation of the stacked
                                        # covariates matrix for 3.1
                                        COs <- do.call("rbind", rep(list(CO), ud))
                                        # Concatenating CD and COs into CD
                                        CD <- cbind(CD, COs)
                                    }
                                    # Else, in case CO is empty
                                    # (i.e. we don't consider
                                    # any covariate)
                                    # simply continue with the current CD

                                    # Eventually concatenating CD with
                                    # COtselected (the matrix containing the
                                    # current time-dependent covariates)
                                    # Checking if COt is NOT completely empty
                                    if (ncot > 0) {
                                        # Concatenating CD and COtselected into
                                        # CD
                                        CD <- cbind(CD, as.data.frame(COtselected))
                                    }


                                    # PAST and FUTURE
                                } else {
                                    # meaning np>0 and nf_temp>0 and that, thus,
                                    # PAST as well as FUTURE VIs do exist
                                    # initialisation of matrices CDp and CDf
                                    CDp <- matrix(NA, nrow=nr*ud, ncol=np)
                                    CDf <- matrix(NA, nrow=nr*ud, ncol=nf_temp)

                                    if (ncot > 0) {
                                        # initialisation of matrix COtselected
                                        COtselected <- do.call(rbind, replicate(ud, COtsample, simplify=FALSE))
                                    }

                                    # initialisation of matrix CDdb (for past
                                    # distribution analysis)
                                    # (Distribution Before)
                                    if (pastDistrib) {
                                        CDdb <- matrix(NA, nrow=nr*ud, ncol=k)
                                        db <- matrix(NA, nrow=nr, ncol=k)
                                        # submatrix of CDdb: CDdb is composed
                                        # of ud matrix db on top of each other
                                    }

                                    # initialisation of matrix CDda (for future
                                    # distribution analysis)
                                    # (Distribution After)
                                    if (futureDistrib) {
                                        # CDda has same dimensions as CDdb
                                        CDda <- matrix(NA, nrow=nr*ud, ncol=k)

                                        # da has same dimensions as db
                                        da <- matrix(NA, nrow=nr, ncol=k)
                                    }


                                    for (j in frameSize:nc){
                                        t1 <- (nr*(iter-1)+1)
                                        t2 <- nr*iter
                                        # VD
                                        CD[t1:t2,1] <- OD[,j-frameSize+np+1+shift]
                                        # /!\ current pointer on column is thus:
                                        # "j-frameSize+np+1+shift"

                                        # Past VIs
                                        CDp[t1:t2,] <- OD[,(j-frameSize+1):(j-frameSize+np)]
                                        # Future VIs
                                        CDf[t1:t2,] <- OD[,(j-nf_temp+1):j]

                                        # Eventually considering time-dependent
                                        # covariates
                                        if (ncot > 0) {
                                            COttemp <- as.data.frame(matrix(nrow=nr,ncol=0))
                                            for (d in 1:(ncot/nc)) {
                                                COttemp <- cbind(COttemp, COt[,(j-frameSize+np+1+shift) + (d-1)*nc])
                                            }
                                            COtselected[t1:t2,] <- COttemp
                                        }

                                        # Past distribution (i.e. Before)
                                        if (pastDistrib) {
                                            ODt <- t(OD)
                                            ODt <- as.data.frame(ODt)
                                            tempOD <- lapply(ODt[(1:(j-frameSize+np+shift)),], factor, levels=c(1:k,NA), exclude=NULL)
                                            # because:
                                            # j-frameSize+np+1+shift - 1 =
                                            # j-frameSize+np+shift
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
                                            # j-frameSize+np+1+shift + 1 =
                                            # j-frameSize+np+shift+2
                                            da_list <- lapply(tempOD,summary)
                                            da_matrix <- do.call(rbind,da_list)
                                            CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np+shift+2):nc)
                                        }

                                        iter <- iter+1
                                    }


                                    # Concatenating CD
                                    CD <- cbind(CD, CDp, CDf)

                                    if (pastDistrib) {
                                        CD <- cbind(CD, CDdb)
                                    }

                                    if (futureDistrib) {
                                        CD <- cbind(CD, CDda)
                                    }

                                    # Conversion of CD into a data frame
                                    CD <- as.data.frame(CD)

                                    # Eventually concatenating CD with
                                    # COs (the matrix containing the covariates)
                                    if (all(is.na(CO))==FALSE) {
                                        # Checking if CO is NOT
                                        # completely empty
                                        # Creation of the stacked
                                        # covariates matrix for 3.1
                                        COs <- do.call("rbind", rep(list(CO), ud))
                                        # Concatenating CD and COs into CD
                                        CD <- cbind(CD, COs)
                                    }
                                    # Else, in case CO is empty
                                    # (i.e. we don't consider
                                    # any covariate)
                                    # simply continue with the current CD

                                    # Eventually concatenating CD with
                                    # COtselected (the matrix containing the
                                    # current time-dependent covariates)
                                    # Checking if COt is NOT completely empty
                                    if (ncot > 0) {
                                        # Concatenating CD and COtselected into
                                        # CD
                                        CD <- cbind(CD, as.data.frame(COtselected))
                                    }


                                }







                                # 6.3.2.RIGHT Computation of the model (Dealing with the LOCATIONS of imputation) -------------------------------------------------


                                #*********************************************************************************
                                # Initially "computeModel.R"
                                #*********************************************************************************




                                # ==>> Manipulation of parameters

                                # Conversion of CD in a data frame
                                CD <- as.data.frame(CD)

                                # Transformation of the names of the columns
                                # of CD (called V1, V2, ..., "VtotV_temp")
                                colnames(CD) <- paste("V", 1:ncol(CD), sep = "")










                                if (regr == "mlogit") {
                                    ## Case of MULTINOMIAL REGRESSION MODEL

                                    # Linking the package mlogit

                                    # By default, every column of CD are of
                                    # class "numeric".
                                    # Thus, there is no need to convert the
                                    # columns containing the distribution
                                    # data to class "numeric".
                                    # Moreover the class of the covariates
                                    # columns at the very end are ALREADY
                                    # set correctly and we don't need to
                                    # update them.
                                    # On the other hand, the first columns of
                                    # CD (1 up to 1+np+nf) have to be of
                                    # class "factor" (because they are the
                                    # columns containing our categorical
                                    # data coming from OD).

                                    # Transformation of the first columns
                                    # (i.e. the categorical values) of CD
                                    # (column 1 up to column 1+np+nf) into
                                    # factor
                                    CD[,(1:(1+np+nf_temp))] <- lapply(CD[,(1:(1+np+nf_temp))],factor, levels=c(1:k,NA), exclude=NULL)


                                    # Dataframe for mlogit
                                    NCD <- mlogit.data(CD, varying=NULL, choice="V1", shape="wide")


                                    # Computation of the multinomial model
                                    if(totV_temp==1){
                                        # First case is evaluated aside
                                        reglog_6_right <- mlogit(V1~0, data=NCD, reflevel="1")
                                    }

                                    if(totV_temp>1){
                                        # creation of "V2" ... "VtotV_temp"
                                        # (to use them in the formula)
                                        factors_character <- paste("V", 2:totV_temp, sep = "")
                                        # Transformation of this object from
                                        # character to vector (in order to be
                                        # able to access its components)
                                        factors <- as.vector(factors_character)
                                        # Creation of a specific formula
                                        # according to the value of totV_temp
                                        fmla <- as.formula(paste("V1~0|", paste(factors, collapse="+")))
                                        reglog_6_right <- mlogit(fmla, data=NCD, reflevel="1")
                                    }









                                } else if (regr == "lm") {
                                    ## Case of LINEAR REGRESSION MODEL

                                    # Since we are performing a linear
                                    # regression, each element of CD are numbers
                                    # and have then to remain as class
                                    # "numeric" (we don't have to perform some
                                    # class adjustment as for the case
                                    # of the creation of a multinomial model).

                                    # Computation of the linear regression model
                                    if(totV_temp==1){
                                        reglog_6_right <- lm(V1~0, data=CD)
                                        # first case is evaluated aside
                                    }

                                    if(totV_temp>1){
                                        # creation of "V2" ... "VtotV_temp"
                                        # (to use them in the formula)
                                        factors_character <- paste("V", 2:totV_temp, sep = "")
                                        # Transformation of this object from
                                        # character to vector (in order to be
                                        # able to access its components)
                                        factors <- as.vector(factors_character)
                                        # Creation of a specific formula
                                        # according to the value of totV_temp
                                        fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                                        reglog_6_right <- lm(fmla, data=CD)
                                    }











                                } else { # meaning (regr == "lrm")
                                    ## Case of ORDINAL REGRESSION MODEL

                                    # Linking to the package rms to use the
                                    # function "lrm"

                                    # Since we are performing an ordinal
                                    # regression, each element of CD are numbers
                                    # and have then to remain as class "numeric"
                                    # (we don't have to perform some
                                    # class adjustment as for the case of the
                                    # creation of a multinomial model).

                                    # Computation of the ordinal model
                                    if(totV_temp==1){
                                        reglog_6_right <- lrm(V1~0, data=CD)
                                        # first case is evaluated aside
                                    }

                                    if(totV_temp>1){
                                        # creation of "V2" ... "VtotV_temp"
                                        # (to use them in the formula)
                                        factors_character <- paste("V", 2:totV_temp, sep = "")
                                        # Transformation of this object from
                                        # character to vector (in order to be
                                        # able to access its components)
                                        factors <- as.vector(factors_character)
                                        # Creation of a specific formula
                                        # according to the value of totV_temp
                                        fmla <- as.formula(paste("V1~0+", paste(factors, collapse="+")))
                                        reglog_6_right <- lrm(fmla, data=CD)
                                    }

                                }


                                #*********************************************************************************
                                #*********************************************************************************



                                # 6.3.3.RIGHT Imputation using the just created model (Dealing with the actual VALUES to impute) ----------------------------------

                                # Structure and building of the data matrix CDi
                                # The first column of CDi is the dependent
                                # variable (VD, response variable) that we have
                                # to implement during the current iteration
                                # (i.e. automatically a NA).
                                # The following columns are the corresponding
                                # independent variables (VIs, explanatory
                                # variables) coming from the past (np>0) or the
                                # future (nf_temp>0) (ordered by time and
                                # corresponding to the current predictive
                                # pattern) and the distribution of the possible
                                # values (i.e. all possible categorical
                                # variables numbered from 1 to k and of course
                                # the value NA) respectively Before and After
                                # the NA to impute (The fact that every lines of
                                # CDi are identical is related to the working of
                                # the function mlogit that has to have as much
                                # lines in CDi as there are categories of the
                                # variable (see the parameter "k")
                                # --> so, CDi is composed of k identical lines)
                                #
                                #            VD   Past VIs   Future VIS    Past distribution     Future distribution
                                #
                                #         /                                                                             \
                                #        |                                                                              |
                                # CDi =  |  CDi    CDpi       CPfi              CDdbi                   CDdai           |
                                #        \                                                                             /
                                #         \                                                                          /
                                #
                                # We are then going to create a
                                # specific CDi according
                                # to the value of np and nf_temp





                                # Analysing the value of parameter available
                                if (available==TRUE){
                                    # we take the previously imputed
                                    # data into account
                                    LOOKUP <- ODi
                                } else {
                                    # that is available == FALSE and thus we
                                    # don't take the previously imputed
                                    # data into account
                                    LOOKUP <- OD
                                }



                                # Assigning the current "REFORDSLGRight_order"
                                # matrix to the variable matrix REFORDSLGRight
                                # (according to the current value of "order")
                                tempObject = get(paste0("REFORDSLGRight_",order))
                                REFORDSLGRight <- as.matrix(tempObject)
                                if (ncol(REFORDSLGRight) == 1) {
                                    REFORDSLGRight <- t(REFORDSLGRight)
                                }
                                nr_REFORD <- nrow(REFORDSLGRight)



                                if (np==0 & nf_temp>0){
                                    # only FUTURE VIs do exist
                                    for (u in 1:nr_REFORD) {
                                        q <- REFORDSLGRight[u,1]
                                        # taking out the first
                                        # coordinate (row
                                        # number in ORDER) from
                                        # REFORDSLGLeft
                                        w <- REFORDSLGRight[u,2]
                                        # taking out the second
                                        # coordinate (column
                                        # number in ORDER) from
                                        # REFORDSLGLeft

                                        CDi <- matrix(NA, nrow=k, ncol=1)

                                        # Matrix for future values
                                        vect <- LOOKUP[q,(w+1):(w+nf_temp)]
                                        CDfi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)

                                        # Matrix for past distribution
                                        if (pastDistrib) {
                                            dbi <- summary(factor(LOOKUP[q,1:(w-1)], levels=c(1:k), exclude=NULL))/length(1:(w-1))
                                            CDdbi <- matrix(dbi[1:k], nrow=k, ncol=k, byrow=TRUE)
                                        }

                                        # Matrix for future distribution
                                        if (futureDistrib) {
                                            dai <- summary(factor(LOOKUP[q,(w+1):nc], levels=c(1:k), exclude=NULL))/length((w+1):nc)
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

                                        # Type transformation of the columns
                                        # of CDi
                                        # The first values of CDi
                                        # must be of type factor
                                        # (categorical values)
                                        CDi[,(1:(1+np+nf_temp))] <- lapply(CDi[,(1:(1+np+nf_temp))],factor, levels=c(1:k,NA), exclude=NULL)
                                        # The last values of CDi must
                                        # be of type numeric
                                        # (distributions)
                                        if (pastDistrib | futureDistrib) {
                                            CDi[,(1+np+nf_temp+1):ncol(CDi)] <- lapply(CDi[,(1+np+nf_temp+1):ncol(CDi)],as.numeric)
                                        }
                                        # Eventually concatenating CDi with COi
                                        # (the matrix containing the covariates)
                                        if (all(is.na(CO))==FALSE) {
                                            # Checking if CO is
                                            # NOT completely
                                            # empty
                                            # Creation of the matrix COi used
                                            # in 3.3
                                            COi <- do.call(rbind, replicate(k, as.data.frame(CO[i,]), simplify=FALSE))
                                            # Concatenating CDi and COi into CDi
                                            CDi <- cbind(CDi, COi)
                                            # Transformation of the
                                            # names of the columns of CDi
                                            # (called V1, V2, ..., "VtotV")
                                            colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                                        }
                                        # Else, in case CO is empty
                                        # (i.e. we don't consider any covariate)
                                        # simply continue with the current CDi

                                        # Eventually concatenating CDi with
                                        # COtselected_i (the matrix containing
                                        # the current time-dependent covariates)
                                        # Checking if COt is NOT completely
                                        # empty
                                        if (ncot > 0) {
                                            COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
                                            for (d in 1:(ncot/nc)) {
                                                COtselected_i <- cbind(COtselected_i, COt[q,(w) + (d-1)*nc])
                                            }
                                            COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
                                            # Concatenating CDi and
                                            # COtselected_i into CDi
                                            CDi <- cbind(CDi, COtselected_i)
                                            # Transformation of the names of the
                                            # columns of CDi
                                            # (called V1, V2, ..., "VtotV")
                                            colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                                        }


                                        # Check for missing-values among
                                        # predictors
                                        if (max(is.na(CDi[1,2:totV_temp]))==0){


                                            #*******************************************************************
                                            # Initially "imputeValue.R"
                                            #*******************************************************************

                                            if (regr == "mlogit") {
                                                ## Case of MULTINOMIAL
                                                ## REGRESSION MODEL


                                                pred <- predict(reglog_6_right,newdata=CDi)
                                                # Example of value returned by
                                                # pred:
                                                # (Sytematically, the upper
                                                # line represents the
                                                # possible categories of
                                                # the variable
                                                # (here, k=2, so the possible
                                                # categories are either
                                                # 1" or "2"))
                                                #            1            2
                                                # 1.000000e+00 2.111739e-22
                                                #
                                                # Cumulative pred
                                                pred <- cumsum(pred)
                                                # Corresponding example
                                                # value returned by in the
                                                # "cumulative pred":
                                                # 1 2
                                                # 1 1
                                                #
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Checking that the components
                                                # of vector "pred" are
                                                # still included in
                                                # the interval [0,1]
                                                pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                                pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                                # Imputation
                                                alea <- runif(1)
                                                # Example value
                                                # returned in "alea":
                                                # [1] 0.2610005
                                                #
                                                sel <- which(pred>=alea)
                                                # Corresponding example
                                                # value returned in sel:
                                                # 1 2
                                                # 1 2
                                                #















                                            } else if (regr == "lm") {
                                                ## Case of LINEAR REGRESSION
                                                ## MODEL

                                                # Since we are performing a
                                                # linear regression, each
                                                # element of CDi are numbers and
                                                # have to be considered as class
                                                # "numeric"
                                                CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                                pred <- predict(reglog_6_right, CDi)
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Rounding pred to the
                                                # nearest integer
                                                pred <- round(pred)
                                                # Restricting pred to its
                                                # lowest value: 1
                                                pred <- ifelse(pred<1,1,pred)
                                                # Restricting pred to its
                                                # highest value: k
                                                pred <- ifelse(pred>k,k,pred)
                                                sel <- pred













                                            } else { # meaning (regr == "lrm")
                                                ## Case of ORDINAL REGRESSION
                                                ## MODEL

                                                # Since we are performing an
                                                # ordinal regression, each
                                                # element of CDi
                                                # are numbers and have to be
                                                # considered as class "numeric"
                                                CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                                pred <- predict(reglog_6_right, CDi, type="fitted.ind")
                                                # Testing if we are in case
                                                # where k=2 (if this is the
                                                # case, we need to create the
                                                # second complementary
                                                # probility by hand since lrm
                                                # returns only the first
                                                # probability)
                                                if (k == 2) {
                                                    pred <- c(pred,(1-pred))
                                                }
                                                # Cumulative pred
                                                pred <- cumsum(pred)
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Checking that the components
                                                # of vector "pred" are
                                                # still included in the
                                                # interval [0,1]
                                                pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                                pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                                # Imputation
                                                # Since we have introduce a
                                                # noise on the variance,
                                                # it might occur that "alea"
                                                # is greater than the greatest
                                                # value of "pred". We have
                                                # then to restrict "alea"
                                                # to the last value of "pred"
                                                alea <- runif(1)
                                                if (alea > pred[length(pred)]) {
                                                    alea <- pred[length(pred)]
                                                }
                                                sel <- which(pred>=alea)



                                            }

                                            #*******************************************************************
                                            #*******************************************************************


                                            ODi[q,w] <- sel[1]
                                        }


                                    }


                                } else {
                                    # meaning np>0 and nf_temp>0 and that, thus,
                                    # PAST as well as FUTURE VIs do exist
                                    for (u in 1:nr_REFORD) {
                                        q <- REFORDSLGRight[u,1]
                                        # taking out the first
                                        # coordinate (row
                                        # number in ORDER) from
                                        # REFORDSLGLeft
                                        w <- REFORDSLGRight[u,2]
                                        # taking out the second
                                        # coordinate (column
                                        # number in ORDER) from
                                        # REFORDSLGLeft

                                        CDi <- matrix(NA, nrow=k, ncol=1)

                                        # Matrix for past values
                                        vect <- LOOKUP[q,(w-shift-np):(w-shift-1)]
                                        CDpi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)

                                        # Matrix for future values
                                        vect <- LOOKUP[q,(w-shift+MaxGapSLGRight-order+1):(w-shift+MaxGapSLGRight-order+nf_temp)]
                                        CDfi <- matrix(vect, nrow=k, ncol=length(vect), byrow=TRUE)

                                        # Matrix for past distribution
                                        if (pastDistrib) {
                                            dbi <- summary(factor(LOOKUP[q,1:(w-1)], levels=c(1:k), exclude=NULL))/length(1:(w-1))
                                            CDdbi <- matrix(dbi[1:k], nrow=k, ncol=k, byrow=TRUE)
                                        }

                                        # Matrix for future distribution
                                        if (futureDistrib) {
                                            dai <- summary(factor(LOOKUP[q,(w+1):nc], levels=c(1:k), exclude=NULL))/length((w+1):nc)
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

                                        # Type transformation of the columns of
                                        # CDi
                                        # The first values of CDi must
                                        # be of type factor
                                        # (categorical values)
                                        CDi[,(1:(1+np+nf_temp))] <- lapply(CDi[,(1:(1+np+nf_temp))],factor, levels=c(1:k,NA), exclude=NULL)
                                        # The last values of CDi must
                                        # be of type numeric
                                        # (distributions)
                                        if (pastDistrib | futureDistrib) {
                                            CDi[,(1+np+nf_temp+1):ncol(CDi)] <- lapply(CDi[,(1+np+nf_temp+1):ncol(CDi)],as.numeric)
                                        }
                                        # Eventually concatenating CDi with COi
                                        # (the matrix containing the covariates)
                                        if (all(is.na(CO))==FALSE) {
                                            # Checking if CO is
                                            # NOT completely
                                            # empty
                                            # Creation of the matrix COi used
                                            # in 3.3
                                            COi <- do.call(rbind, replicate(k, as.data.frame(CO[i,]), simplify=FALSE))
                                            # Concatenating CDi and COi into CDi
                                            CDi <- cbind(CDi,COi)
                                            # Transformation of the names
                                            # of the columns of CDi
                                            # (called V1, V2, ..., "VtotV")
                                            colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                                        }
                                        # Else, in case CO is empty (i.e. we
                                        # don't consider any covariate)
                                        # simply continue with the current CDi

                                        # Eventually concatenating CDi with
                                        # COtselected_i (the matrix containing
                                        # the current time-dependent covariates)
                                        # Checking if COt is NOT completely
                                        # empty
                                        if (ncot > 0) {
                                            COtselected_i <- as.data.frame(matrix(nrow=1,ncol=0))
                                            for (d in 1:(ncot/nc)) {
                                                COtselected_i <- cbind(COtselected_i, COt[q,(w) + (d-1)*nc])
                                            }
                                            COtselected_i <- do.call(rbind, replicate(k, as.data.frame(COtselected_i[1,]), simplify=FALSE))
                                            # Concatenating CDi and
                                            # COtselected_i into CDi
                                            CDi <- cbind(CDi, COtselected_i)
                                            # Transformation of the names of the
                                            # columns of CDi
                                            # (called V1, V2, ..., "VtotV")
                                            colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
                                        }


                                        # Check for missing-values among
                                        # predictors (i.e. we won't impute any
                                        # value on the current MD if there is
                                        # any NA among the VIs)
                                        if (max(is.na(CDi[1,2:totV_temp]))==0){
                                            # checking that there is no NA among
                                            # the current VI (otherwise no data
                                            # will be imputed for the current
                                            # NA)


                                            #*******************************************************************
                                            # Initially "imputeValue.R"
                                            #*******************************************************************

                                            if (regr == "mlogit") {
                                                ## Case of MULTINOMIAL
                                                ## REGRESSION MODEL


                                                pred <- predict(reglog_6_right,newdata=CDi)
                                                # Example of value returned by
                                                # pred:
                                                # (Sytematically, the upper line
                                                # represents the
                                                # possible categories
                                                # of the variable
                                                # (here, k=2, so the
                                                # possible categories
                                                # are either
                                                # 1" or "2"))
                                                #            1            2
                                                # 1.000000e+00 2.111739e-22
                                                #
                                                # Cumulative pred
                                                pred <- cumsum(pred)
                                                # Corresponding example value
                                                # returned by in the
                                                # "cumulative pred":
                                                # 1 2
                                                # 1 1
                                                #
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Checking that the components
                                                # of vector "pred" are
                                                # still included in the
                                                # interval [0,1]
                                                pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                                pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                                # Imputation
                                                alea <- runif(1)
                                                # Example value returned in
                                                # "alea":
                                                # [1] 0.2610005
                                                #
                                                sel <- which(pred>=alea)
                                                # Corresponding example value
                                                # returned in sel:
                                                # 1 2
                                                # 1 2
                                                #















                                            } else if (regr == "lm") {
                                                ## Case of LINEAR REGRESSION
                                                ## MODEL

                                                # Since we are performing a
                                                # linear regression, each
                                                # element of CDi are numbers and
                                                # have to be considered as class
                                                # "numeric"
                                                CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                                pred <- predict(reglog_6_right, CDi)
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Rounding pred to the
                                                # nearest integer
                                                pred <- round(pred)
                                                # Restricting pred to its
                                                # lowest value: 1
                                                pred <- ifelse(pred<1,1,pred)
                                                # Restricting pred to its
                                                # highest value: k
                                                pred <- ifelse(pred>k,k,pred)
                                                sel <- pred













                                            } else { # meaning (regr == "lrm")
                                                ## Case of ORDINAL REGRESSION
                                                ## MODEL

                                                # Since we are performing an
                                                # ordinal regression, each
                                                # element of CDi
                                                # are numbers and have to be
                                                # considered as class "numeric"
                                                CDi <- as.data.frame(lapply(CDi[,1:ncol(CDi)],as.numeric))

                                                pred <- predict(reglog_6_right, CDi, type="fitted.ind")
                                                # Testing if we are in case
                                                # where k=2 (if this is the
                                                # case, we need to create the
                                                # second complementary
                                                # probility by hand since lrm
                                                # returns only the first
                                                # probability)
                                                if (k == 2) {
                                                    pred <- c(pred,(1-pred))
                                                }
                                                # Cumulative pred
                                                pred <- cumsum(pred)
                                                # Introducing a variance "noise"
                                                pred <- rnorm(length(pred),pred,noise)
                                                # Checking that the components
                                                # of vector "pred" are
                                                # still included in
                                                # the interval [0,1]
                                                pred <- unlist(lapply(pred, function(x) {if (x < 0) 0 else x}))
                                                pred <- unlist(lapply(pred, function(x) {if (x > 1) 1 else x}))
                                                # Imputation
                                                # Since we have introduce a
                                                # noise on the variance, it
                                                # might occur that "alea" is
                                                # greater than the
                                                # greatest value of "pred".
                                                # We have then to restrict
                                                # "alea" to the last value of
                                                # "pred"
                                                alea <- runif(1)
                                                if (alea > pred[length(pred)]) {
                                                    alea <- pred[length(pred)]
                                                }
                                                sel <- which(pred>=alea)



                                            }

                                            #*******************************************************************
                                            #*******************************************************************


                                            ODi[q,w] <- sel[1]
                                        }
                                    }

                                }

                            }

                        }

                    }

                }









                # Trying to catch the potential singularity error (part 2/2)
                # (comment this part of code to debug more easily and see the
                # full message of any occuring error)
                #********************************************************************************************************************************
            },




            error=function(error_message) {
                # Catching the error that we want...
                if ( substr(error_message[1],1,14) == "Lapack routine" | substr(error_message[1],1,34) == "system is computationally singular" ) {
                    message("/!\\ Our model is too specified! It consists of too many independant variables and","\n",
                            "is thus 100% predictable. By inversing the data matrix, R finds a determinant very","\n",
                            "close to 0 and we meet a singularity.","\n",
                            "In order to avoid this error, you need to lower the value of np and/or nf.","\n",
                            "(By the way, you can also try to relaunch the program with the same settings and it","\n",
                            "might work on a second or third run... (this is due to randomly generated numbers","\n",
                            "to complete the prediction))","\n\n",
                            "Below is the error message from R:")
                    message(error_message)

                    # ... or simply displays the other error types
                } else {
                    message(error_message)
                }

            }



        )
        #****************************************************************************************************************************************





        # Updating the matrix RESULT used to store the multiple imputations
        RESULT <- rbind(RESULT,cbind(replicate(nr,o),ODi))











    }





















    # X. Final conversions ----------------------------------------------------------------------------------------------------------------------------------------






    ## Adjustment of the rendered form of RESULT
    if (mi.return == 1) {    # case mi.return == 1 (not including initial data
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
    # Else (meaning that we are in the case mi.return == 2 (including initial
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
        if (mi == 1 & mi.return == 1) {
            RESULT <- as.data.frame( sapply(RESULT, mapvalues,
                                            from = as.character(as.vector(1:length(ODlevels))),
                                            to = ODlevels) )
        } else {
            # Taking account ot the special notation of RESULT that has an extra
            # column on the left of RESULT (as soon as mi > 1 or in any case if
            # mi.return == 2)
            RESULT[,2:ncol(RESULT)] <- as.data.frame( sapply(RESULT[,2:ncol(RESULT)],
                                                             mapvalues,
                                                             from = as.character(as.vector(1:length(ODlevels))),
                                                             to = ODlevels) )
        }
    }
    #*************************************












    # Returning the matrix composed of every imputations
    return(RESULT)





}
