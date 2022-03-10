# 6. Imputation of SLG NAs
#
#
################################################################################
# Left-hand side SLG imputation

LSLGNAsImpute <- function(OD, ODi, CO, COt, COtsample, ORDERSLG, pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud, ncot, nco, k, noise, available,num.trees,min.node.size,max.depth,timing){     # Checking if we have to impute
  # left-hand side SLG
  
  warning("/!\\ Specially Located Gaps (SLG) have been detected on the left-hand side of OD.","\n","For certain missing data groups close to the border of OD, np may have been automatically reduced.","\n","If you don't want this to happen, please choose a lower value for np.")
  
  # 6.2.LEFT Computation of the order of imputation of each MD ----

  # Initialization of a list to take all the variable returned by the functions
  ParamList = list()
  # Creation of the temporary SLG matrices for the left-hand
  # side of OD
  
  for (h in 2:np) {
    if (max(ORDERSLG[,h])>0) {
      # Checking if a gap begins
      # somewhere on the current column
      # If that is not the case, there is
      # no need to perform
      # the entire following procedure
      # and we simply can go
      # to the next column of ORDERSLGLeft
      
      ParamList[c("ORDERSLG_temp","totV_temp","np_temp")] <- SLGMatrix_temp(nr, nc, nf, h, ORDERSLG, nco, ncot, pastDistrib, futureDistrib, k)
      
      if (max(ParamList$ORDERSLG_temp)==0) {
        next
      }
      
      # In a similar manner to part 2.4., we go here one
      # single time through the reduced version
      # ORDERSLGLeft_temp of ORDERSLG and we create
      # MaxGapSLGLeft "REFORDSLGRLeft_" matrices
      # collecting the coordinates of each corresponding
      # values in ORDERSLGLeft_temp which are greater
      # than 0
      
      
      # REFORDSLGLeft matrices
      # Initialization of the REFORDSLGLeft matrices
      ParamList[c("MaxGap","REFORD_L","ORDERSLG_temp")]  <- REFORDInit(ParamList$ORDERSLG_temp, nr, nc)
      
      # MaxGapSLGLeft <- REFORDOD[[1]]   
      # REFORDSLG_L <- REFORDOD[[2]]
      
      # 6.3.LEFT Imputation of the missing data listed by ORDERSLGLeft_temp and ORDERSLGRight_temp using a specific model ----
      
      # 6.3.1.LEFT Building of the different data matrices CD ----
      # for the computation of the model for every SLG
      # on the left-hand side of OD
      for (order in 1:ParamList$MaxGap) {
        ParamList[c("CD","shift")]  <- SLGCDMatBuild(CO, COt, OD, order, ParamList$MaxGap, ParamList$np_temp, ncot, nr, nc, nf, COtsample, pastDistrib, futureDistrib, k,timing)
        # 6.3.2.LEFT Computation of the model (Dealing with the LOCATIONS of imputation) ----
        log_CD <- list()
        log_CD[c("reglog","CD")] <- ComputeModel(ParamList$CD, regr, ParamList$totV_temp, ParamList$np_temp,nf, k,num.trees,min.node.size,max.depth,timing)
        # 6.3.3.LEFT Imputation using the just created model (Dealing with the actual VALUES to impute) ----
        ODi <- SLGCreatedModelImpute(CO, OD, log_CD$CD,  ODi, COt, ncot, nf, nc, regr, k, ParamList$totV_temp, log_CD$reglog, noise, available, ParamList$REFORD_L, ParamList$np_temp, pastDistrib, futureDistrib, order, ParamList$MaxGap, ParamList$shift,timing)
      }
     
    }
    
  }
  return(ODi)
}
  


##############################################################################
#Right-hand side SLG imputation

RSLGNAsImpute <- function(OD, ODi, CO, COt, COtsample, ORDERSLGRight, pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud, ncot, nco, k, noise, available,num.trees,min.node.size,max.depth,timing){
  # Checking if we have to impute right-hand
  # side SLG
  
  warning("/!\\ Specially Located Gaps (SLG) have been detected on the right-hand side of OD.","\n","For certain missing data groups close to the border of OD, nf may have been automatically reduced.","\n","If you don't want this to happen, please choose a lower value for nf.")
  
  # 6.2.RIGHT Computation of the order of imputation of each MD ----------------
  
  # Initialization of a list to take all the variable returned by the functions
  ParamList = list()
  
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
      
      ParamList[c("ORDERSLGRight_temp","totV_temp","nf_temp")] <- SLGMatrixRight_temp(nr, nc, np, h, ORDERSLGRight, nco, ncot, pastDistrib, futureDistrib, k)
      
      if (max(ParamList$ORDERSLGRight_temp)==0) {
        next
      }
      
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
      
      
      ParamList[c("MaxGap","REFORD_L","ORDERSLGRight_temp")] <- REFORDInit(ParamList$ORDERSLGRight_temp, nr, nc)
   
      
      
      # 6.3.RIGHT Imputation of the missing data listed by ORDERSLGLeft_temp and ORDERSLGRight_temp using a specific model ----
      
      for (order in 1:ParamList$MaxGap) {
        # 6.3.1.RIGHT Building of the different data matrices CD ---------------
        # for the computation of the model for every SLG
        # on the right-hand side of OD
        ParamList[c("CD","shift")] <- SLGCDMatBuild(CO, COt, OD, order, ParamList$MaxGap, np, ncot, nr, nc, ParamList$nf_temp, COtsample, pastDistrib, futureDistrib, k,timing)

        # 6.3.2.RIGHT Computation of the model (Dealing with the LOCATIONS of imputation) ----
        log_CD <- list()
        log_CD[c("reglog","CD")] <- ComputeModel(ParamList$CD, regr, ParamList$totV_temp, np,ParamList$nf_temp, k,num.trees,min.node.size,max.depth,timing)
        

        # 6.3.3.RIGHT Imputation using the just created model (Dealing with the actual VALUES to impute) ----
        ODi <- SLGCreatedModelImpute(CO, OD, log_CD$CD, ODi, COt, ncot, ParamList$nf_temp, nc, regr, k, ParamList$totV_temp, log_CD$reglog, noise, available, ParamList$REFORD_L, np, pastDistrib, futureDistrib, order, ParamList$MaxGap, ParamList$shift,timing)
 
        }
    
      }
    
    }
    return(ODi)
}





SLGMatrix_temp <- function(nr, nc, nf, h, ORDERSLG, nco, ncot, pastDistrib, futureDistrib, k ){
  # Creation of the temporary SLG matrices for the left-hand
  # side of OD
  
  # initialization of a new temporary
  # ORDERSLG_temp matrix
  ORDERSLG_temp <- matrix(0,nrow=nr,ncol=nc)
  
  for (i in 1:nr) { # going from top to bottom
    j <- h        # storing the current column
    # we are checking (and
    # reinitializing j)
    
    if (ORDERSLG[i,j]>0 & ORDERSLG[i,j-1]==0) {
      
      while (ORDERSLG[i,j]>0) {
        # going from left to right
        ORDERSLG_temp[i,j] <- ORDERSLG[i,j]
        j <- j+1
      }
    }
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
  
  
  # (i.e. updating matrix ORDERSLG_temp with
  # coherent value of "order"
  # (going from 1 to MaxGap))
  
  # Adjusting the matrix ORDERSLG_temp and
  # storing the coordinates of the NA to impute
  return(list(ORDERSLG_temp, totV_temp, np_temp))
}
    
    


SLGMatrixRight_temp <- function(nr, nc, np, h, ORDERSLGRight, nco, ncot, pastDistrib, futureDistrib, k ){
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
  return(list(ORDERSLGRight_temp, totV_temp, nf_temp))
}


################################################################################
# Compute the CD matrix for SLG

SLGCDMatBuild <- function(CO,COt, OD, order, MaxGapSLGLeft, np, ncot, nr, nc, nf, COtsample, pastDistrib, futureDistrib, k,timing){
  # Building of a data matrix for the computation
  # of the model

  ud <- nc-(MaxGapSLGLeft-order+np+nf)
  # number of usable data for each row of OD
  # size of the current mobile caracteristic frame
  # (that changes according to "order") which is
  # equal to nc-ud+1
  frameSize <- MaxGapSLGLeft-order+np+nf+1
  # Structure and building of the data matrix CD
  # The first column of CD is the dependent
  # variable (VD, response variable)
  # The following columns are the independent
  # variables (VIs, explanatory variables) coming
  # from the past (np>0) or the future (nf>0)
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
  # according to the value of np and nf
  
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
  if ( (np > 0 & nf > 0) & ( (MaxGapSLGLeft%%2==0 & order%%2==0) | (MaxGapSLGLeft%%2!=0 & order%%2!=0) )){
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
  # meet (because np has to be greater than
  # 0, otherwise, it would mean that we are not
  # in the case of a SLG and that we
  # can impute it as a usual internal gap)
  
  # Only PAST
  if (np>0 & nf==0) {
    CD <- PastVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, COtsample, COt, pastDistrib, futureDistrib, k,timing)
    # PAST and FUTURE
  } else {
    CD <- PastFutureVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, COtsample, COt, pastDistrib, futureDistrib, k, nf, shift,timing)
  }
 
  return(list(CD, shift))
}



################################################################################
# Impute SLG using created model

SLGCreatedModelImpute <- function(CO, OD, CD, ODi, COt, ncot, nf, nc, regr, k, totV_temp, reglog_6, noise, available, REFORDSLG_L, np, pastDistrib, futureDistrib, order, MaxGap, shift,timing){
    
    # Structure and building of the data matrix CDi
    # The first column of CDi is the dependent
    # variable (VD, response variable) that we have
    # to implement during the current iteration
    # (i.e. automatically an NA).
    # The following columns are the corresponding
    # independent variables
    # (VIs, explanatory variables)
    # coming from the past (np>0) or the future
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
    # according to the value of np and nf
    
    
    
    
    
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
    
    REFORDSLGLeft <- as.matrix(REFORDSLG_L[[order]])
    if (ncol(REFORDSLGLeft) == 1) {
      REFORDSLGLeft <- t(REFORDSLGLeft)
    }
    nr_REFORD <- nrow(REFORDSLGLeft)
    
    
    
    if (np>0 & nf==0){
      # only PAST VIs do existe)
      ODi <- ODiImputePAST(CO, ODi, CD, COt, REFORDSLGLeft, nr_REFORD, pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV_temp, reglog_6, LOOKUP, regr, noise,timing)
        
  
    } else {
      # meaning np>0 and nf>0 and that, thus,
      # PAST as well as FUTURE VIs do exist
      ODi <- ODiImputePF(CO, ODi, CD, COt, REFORDSLGLeft, nr_REFORD, pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV_temp, reglog_6, LOOKUP, regr, noise, shift, MaxGap, order,timing)
    }
    return(ODi)
}








# 
# 
# ################################################################################
# # Long Gap imputation
# #
# LongGapImpute <- function(OD, ODi, CO, COtsample, ORDERSLGBoth, pastDistrib, futureDistrib, regr, np, nr, nc, ud, ncot, nco, k, noise, available){
# 
#   nfix<-1
#   # 6.2.Both Computation of the order of imputation of each MD ----------------------------------------------------------------------------------
# 
#   # Creation of the temporary SLG matrices for the Both-hand
#   # side of OD
#   for (h in 2:np) {
#     if (max(ORDERSLGBoth[,h])>0) {
#       # Checking if a gap begins
#       # somewhere on the current column
#       # If that is not the case, there is
#       # no need to perform
#       # the entire following procedure
#       # and we simply can go
#       # to the next column of ORDERSLGBoth
# 
#       ParamList[c("ORDERSLGBoth_temp","totV_temp","np_temp")] <- SLGMatrix_temp(nr, nc, nfix, h, ORDERSLGBoth, nco, ncot, pastDistrib, futureDistrib, k)
# 
# 
# 
#       # In a similar manner to part 2.4., we go here one
#       # single time through the reduced version
#       # ORDERSLGBoth_temp of ORDERSLG and we create
#       # MaxGapSLGBoth "REFORDSLGRBoth_" matrices
#       # collecting the coordinates of each corresponding
#       # values in ORDERSLGBoth_temp which are greater
#       # than 0
# 
# 
#       # REFORDSLGBoth matrices
# 
#       ParamList[c("MaxGap","REFORD_L", "ORDERSLGBoth_temp")] <- REFORDInit(ParamList$ORDERSLGBoth_temp, nr, nc)
# 
# 
# 
#       # 6.3.Both Imputation of the missing data listed by ORDERSLGBoth_temp and ORDERSLGRight_temp using a specific model -------------------
# 
# 
#       # 6.3.1.Both Building of the different data matrices CD -------------------------------------------------------------------------------
#       # for the computation of the model for every SLG
#       # on the Both-hand side of OD
# 
#       for (order in 1:ParamList$MaxGap) { # /!\ "order"
#         # corresponds to the values
#         # of the components of
#         # ORDERSLGBoth_temp
#         # (i.e. the number of the iteration, the order
#         # in which the values are going to be imputed)
#         ParamList[c("CD","shift")] <- SLGCDMatBuild(CO, OD, ODt, order, ParamList$MaxGap, ParamList$np, ncot, nr, nc, nfix, COtsample, pastDistrib, futureDistrib, k)
# 
#         # 6.3.2. Both Computation of the model (Dealing with the LOCATIONS of imputation)
#         reglog_6_both <- ComputeModel(ParamList$CD, regr, ParamList$totV_temp, ParamList$np_temp+nfix, k)
# 
#         # 6.3.3. lEFT and RIGHT  Imputation using the just created model (Dealing with the actual VALUES to impute)
#         ODi <- SLGCreatedModelImpute(CO, OD, ODi, COt, ncot, nfix, nc, regr, k, ParamList$totV_temp, reglog_6_both, noise, available, ParamList$REFORD_L, ParamList$np_temp, pastDistrib, futureDistrib, order, ParamList$MaxGap, ParamList$shift)
# 
# 
# 
# 
# 
#             }
#           }
# 
#         }
# 
#       }
