#' Visualizing missing patterns among a dataset
#'
#' \code{seqMissPatt.R} offers a quick visualization of the repeating missing
#' patterns among the dataset especially thanks to the package \code{TraMineR}.
#'
#' @param OD \code{matrix} object containing sequences of a variable with missing data (coded as \code{NA}).
#' @param k \code{numeric} object corresponding to the number of categories of the multinomial variable numbered from \code{1} to \code{k}.
#' @param numbOfMostFreqSeq \code{numeric} object defining the number of most frequent sequences \code{TraMineR} plots are going to display (default \code{10}).
#' @param pbarw \code{logical} object making the width of the bars of the \code{TraMineR} plots proportional to the sequence frequency in the dataset if \code{pbarw=TRUE} (default \code{FALSE}).
#' @param clustNumb \code{numeric} object corresponding to the number of clusters used to characterize the original dataset \code{OD} (default \code{3}) (min value: \code{2}, max value: \code{9}).
#'
#' @author Andre Berchtold <andre.berchtold@@unil.ch>
#'
#' @return Various plots indexing the missing patterns among the binary matrix
#' of the inputed dataset:
#' \itemize{
#'   \item Sequences ordered according to ascending size of initial/terminal gaps
#'   \item Most frequent sequences
#'   \item Clustered sequences
#'   \item Etc.
#' }
#'
#' @examples
#' data(OD)
#'
#' seqMissPatt(OD=OD, k=2, numbOfMostFreqSeq=5, pbarw=TRUE, clustNumb=3)
#'
#'
#' @export


seqMissPatt <- function(OD, k, numbOfMostFreqSeq=10, pbarw=FALSE, clustNumb=3) {

# test

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



    # 1.4 Making sure that numbOfMostFreqSeq is greater than 0 and lower or equal to the number of row of the dataset -----------------------------------
    numbOfMostFreqSeq <- abs(numbOfMostFreqSeq)
    if (numbOfMostFreqSeq<1 | numbOfMostFreqSeq>nrow(OD)) {
        stop("Please choose an integer value situated between 1 and the total number of rows of your original dataset for the parameter 'numbOfMostFreqSeq' (the number of most frequent sequences).")
    }


    # 1.5 Verifying that clustNumb is not smaller than 2 nor greater than 9 -----------------------------------------------------------------------------
    if (clustNumb%%1!=0 | clustNumb<2 | clustNumb>9) {
        stop("Please choose an integer value situated between 2 and 9 for the parameter clustNumb (the number of clusters used in the characterization of OD).")
    }













    # 2. Detailed map of the NA in OD -------------------------------------------------------------------------------------------------------------------
    # 2.1 Overall NA indexing ---------------------------------------------------------------------------------------------------------------------------

    # Replacing OD as a single vector
    # We have then like one single column of nrxnc rows
    # (every initial column of OD are then stacked on each other)
    # i.e. every nr+1 (i.e. 1, nr+1, (2*nr)+1, (3*nr)+1, etc.)
    # element we are situated like
    # at the top beginning of a new column of OD
    ODvector <- as.vector(OD)

    NaReplacedBy1s_col <- ifelse(is.na(ODvector),1,0)







    # 2.2 Total number of NAs and occurrence percentage -------------------------------------------------------------------------------------------------
    totSumNA <- sum(NaReplacedBy1s_col)

    percNA <- paste(round(totSumNA/(nr*nc)*100,3),"%")









    # 2.3 Number of NAs by column and occurrence percentage ---------------------------------------------------------------------------------------------
    # sumNA vectors (one for every column of OD)
    shift <- 1
    for (i in 1:nc) {
        assign(paste("sumNA_col_",i,sep=''),sum(NaReplacedBy1s_col[shift:(i*nr)]))
        tempObject = get(paste0("sumNA_col_",i,sep=''))
        assign(paste("percNA_col_",i,sep=''),round(tempObject/nr,3))
        shift <- shift+nr
    }









    # 2.4 Number of NAs by row and occurence percentage -------------------------------------------------------------------------------------------------
    # sumNA vectors (one for every row of OD)
    ODtvector <- as.vector(t(OD)) # transposing OD and taking its
    # corresponding vector composed of its rows
    NaReplacedBy1s_row <- ifelse(is.na(ODtvector),1,0)
    # Every nc+1 (i.e. 1, nc+1, (2*nc)+1, (3*nc)+1, etc.) element we are
    # situated like at the beginning of a new row of OD
    shift <- 1
    for (i in 1:nr) {
        assign(paste("sumNA_row_",i,sep=''),sum(NaReplacedBy1s_row[shift:(i*nc)]))
        tempObject = get(paste0("sumNA_row_",i,sep=''))
        assign(paste("percNA_row_",i,sep=''),round(tempObject/nc,3))
        shift <- shift+nc
    }








    # 2.5 Number of sequences (rows) wiht and without NAs
    numbOfSeqWithNA <- na.count(apply(OD,1,max))
    numbOfSeqWithoutNA <- nr - numbOfSeqWithNA






    # Display result with Amelia
    #*************************************************
    
    # Creating the label for the y axis
    y.labels <- c()
    x.labels <- c()
    # y.labels
    for (i in 1:nrow(OD)) {
        tempObject <- get(paste0("sumNA_row_",i))
        tempObject_2 <- get(paste0("percNA_row_",i))
        y.labels_current = paste('row.',i,'.NA.',tempObject,'.perc.',tempObject_2*100,'.',sep='')
        y.labels <- rbind(y.labels,y.labels_current)
    }
    y.labels <- as.vector(y.labels)
    y.at <- c(1:nr)
    # x.labels
    for (i in 1:nc) {
        tempObject <- get(paste0("sumNA_col_",i))
        tempObject_2 <- get(paste0("percNA_col_",i))
        x.labels_current = paste('col',i,' NA.',tempObject,' perc.',tempObject_2*100,sep='')
        x.labels <- rbind(x.labels,x.labels_current)
    }
    colnames(OD) <- x.labels
    for (i in 1:nc) {
        colnames(OD)[i] <- paste(colnames(OD)[i],"%",sep='')
    }
    OD <- data.frame(OD)
    missmap(OD, vars=1:ncol(OD), legend=TRUE, col=c("black","grey"), main="Overview of the missing data", x.cex=0.6, y.cex=0.7, y.labels, y.at, rank.order=FALSE, margins=c(10,10))
    #*************************************************




















    # 3. Ascending visualization of initial gaps --------------------------------------------------------------------------------------------------------

    # Creation of matrix ORDER (= corresponding binary matrix of OD
    # ('0' = observed data, '1' = missing data))
    # Other way to create such a corresponding binary matrix:
    #ODb <- apply( OD, 2, function(x) {ifelse((is.na(x)), 1, 0)} )
    # initialization of matrix ORDER with 0 everywhere
    ORDER <- matrix(0,nr,nc)
    # creation of matrix SEL, constituted of TRUE where there is MD in OD and
    # of FALSE everywhere else
    SEL <- is.na(OD)==TRUE
    # setting some 1 in ORDER at the location where in SEL we have some TRUE
    ORDER[SEL] <- 1


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


    # Combining InitGapSize with a vector containing the position of the
    # rows in OD
    RowAndInitGapSize <- as.data.frame(cbind(c(1:nr),InitGapSize))

    InitGapSizeOrdered <- RowAndInitGapSize[order(RowAndInitGapSize$InitGapSize),]

    # Reoordering the row of OD according to the ascending number of NA among
    # initial gaps
    ODInitOrdered <- matrix(NA,nrow=nr,ncol=nc)
    for (i in 1:nr) {
        ODInitOrdered[i,] <- as.matrix(OD[InitGapSizeOrdered[i,1],])
    }


    # Display result with Amelia
    #*************************************************
    
    # Creating the label for the y axis
    y.labels <- c()
    # y.labels
    for (i in 1:nr) {
        y.labels_current = paste('row',InitGapSizeOrdered[i,1],collapse=' ')
        y.labels <- rbind(y.labels,y.labels_current)
    }
    y.labels <- as.vector(y.labels)
    y.at <- c(1:nr)
    # Saving ODInitOrdered is possible by uncommenting the lines
    # "jpeg('OD.jpg')" and "dev.off()"
    ODInitOrdered <- data.frame(ODInitOrdered)
    missmap(ODInitOrdered, vars=1:ncol(ODInitOrdered), legend=TRUE, col=c("black","grey"), main="Ascending visualization of initial gaps", x.cex=0.6, y.cex=0.7, y.labels, y.at, rank.order=FALSE, margins=c(10,10))
    #*************************************************
































    # 4. Ascending visualization of terminal gaps -------------------------------------------------------------------------------------------------------

    # Creation of vector TermGapSize (i.e. a vector containing the size of
    # the terminal gaps of each line)
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


    # Combining TermGapSize with a vector containing the position of the
    # rows in OD
    RowAndTermGapSize <- as.data.frame(cbind(c(1:nr),TermGapSize))

    TermGapSizeOrdered <- RowAndTermGapSize[order(RowAndTermGapSize$TermGapSize),]

    # Reoordering the row of OD according to the ascending number of NA among
    # initial gaps
    ODTermOrdered <- matrix(NA,nrow=nr,ncol=nc)
    for (i in 1:nr) {
        ODTermOrdered[i,] <- as.matrix(OD[TermGapSizeOrdered[i,1],])
    }


    # Display result with Amelia
    #*************************************************
    
    # Creating the label for the y axis
    y.labels <- c()
    # y.labels
    for (i in 1:nr) {
        y.labels_current = paste('row',TermGapSizeOrdered[i,1],collapse=' ')
        y.labels <- rbind(y.labels,y.labels_current)
    }
    y.labels <- as.vector(y.labels)
    y.at <- c(1:nr)
    # Saving ODTermOrdered is possible by uncommenting the lines
    # "jpeg('OD.jpg')" and "dev.off()"
    ODTermOrdered <- data.frame(ODTermOrdered)
    missmap(ODTermOrdered, vars=1:ncol(ODTermOrdered), legend=TRUE, col=c("black","grey"), main="Ascending visualization of terminal gaps", x.cex=0.6, y.cex=0.7, y.labels, y.at, rank.order=FALSE, margins=c(10,10))
    #*************************************************






















    # 5. Exploring patterns with TraMineR ---------------------------------------------------------------------------------------------------------------
    #
    #
    #
    #
    # Info taken from: http://traminer.unige.ch/preview.shtml
    #
    # Examples of possible plots:
    #       - Index of first 10 sequences
    #       - Index plot of all sequences sorted by ...
    #       - /!| Frequency lot of 10 most frequent sequences
    #       - /!| State distribution plot by time point
    #       - Legend
    #       - Histogram of sequences turbulence
    #

    # Defining 'idxs' based on numbOfMostFreqSeq
    idxs <- 1:numbOfMostFreqSeq

    ORDER.labels <- c("observed","missing")
    ORDER.scodes <- c("P","M")
    ORDER.seq <- seqdef(ORDER, states=ORDER.scodes, labels=ORDER.labels)
    #
    #**************************
    # /!| Draw the sequence frequency plot of the 10 most frequent sequences
    # with bar width proportional to the frequencies
    seqfplot(ORDER.seq, group=NULL, idxs=get("idxs"), pbarw, with.legend=T)
    #
    #**************************
    # /!| Plot the state distribution by time points
    seqdplot(ORDER.seq, with.legend=T)
    #
    # Compute the optimal matching distances using substitution
    # costs based on transition
    # rates observed in the data and a 1 indel cost.
    # The resulting distance matrix is
    # stored in the dist.om1 object
    submat <- seqsubm(ORDER.seq, method="TRATE")
    dist.om1 <- seqdist(ORDER.seq, method="OM", indel=1, sm=submat)
    #
    # Make a typology of the trajectories:
    # load the cluster library, build a Ward
    # hierarchical clustering of the sequences
    # from the optimal matching distances
    # and retrieve for each individual sequence
    # the cluster membership of the 3 class
    # solution
    clusterward1 <- agnes(dist.om1, diss=T, method="ward")
    plot(clusterward1)
    # cl1.3 <- cutree(clusterward1, k=3)
    # cl1.3fac <- factor(cl1.3, labels=c("Type 1", "Type 2", "Type 3"))
    cl1.3 <- cutree(clusterward1, get("clustNumb"))
    labelsVect <- c()
    for (i in 1:clustNumb) {
        labelsVect[i] <- paste0("Type ",i,sep='')
    }
    cl1.3fac <- factor(cl1.3, labels=labelsVect)
    #
    # /!| Plot the 10 most frequent sequences of each cluster
    seqfplot(ORDER.seq, group=cl1.3fac, idxs=get("idxs"), pbarw=T)

    #**************************
    # /!| Plot the state distribution within each cluster
    seqdplot(ORDER.seq, group=cl1.3fac)
    #
    #**************************






}







