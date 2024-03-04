#' Generation of missing data under the form of gaps, which
#' is the typical form of missing data with longitudinal data.
#' A missing completely at random (MAR) mechanism is used.
#'
#' @param data either a data frame containing sequences of a multinomial variable with missing data (coded as \code{NA}) or
#' a state sequence object built with the TraMineR package
#' @param states.high list of states that will have a larger probability to trigger a
#' subsequent gap of missing data
#' @param pstart.high probability to start a missing data for the specified states
#' @param pstart.low probability to start a missing data for the other states
#' @param propdata proportion of the observations on which missing data will be simulated
#' @param maxgap maximum length of a gap of missing data
#'
#' @return Returns either a data frame or a state sequence object, depending
#' the type of data that was provided to the function
#'
#' @export
seqaddNA <- function(data, states.high = NULL, propdata = 1, pstart.high = 0.1, pstart.low = 0.005, maxgap = 3) {
  sizehalf <- round(propdata * nrow(data))
  rowsmiss <- sample(1:nrow(data), size = sizehalf, replace = FALSE)
  matrix.missing <- matrix(1, nrow(data), ncol(data))
  for (i in 1:length(rowsmiss)) {
    nmis <- ncol(data)
    length.gap <- 0
    while (nmis > floor(0.75 * ncol(data))) {
      for (j in 1:ncol(data)) {
        if (length.gap == maxgap) {
          matrix.missing[rowsmiss[i], j] <- 1
          length.gap <- 0
        } else {
          if (j == 1) {
            matrix.missing[rowsmiss[i], j] <- sample(x = c(0, 1), size = 1, prob = c(pstart.low, 1 - pstart.low))
          } else {
            if (matrix.missing[rowsmiss[i], j - 1] == 1) {
              if (data[rowsmiss[i], j - 1] %in% states.high) {
                matrix.missing[rowsmiss[i], j] <- sample(x = c(0, 1), size = 1, prob = c(pstart.high, 1 - pstart.high))
              } else {
                matrix.missing[rowsmiss[i], j] <- sample(x = c(0, 1), size = 1, prob = c(pstart.low, 1 - pstart.low))
              }
            } else {
              matrix.missing[rowsmiss[i], j] <- sample(x = c(0, 1), size = 1, prob = c(66, 34))
            }
          }
          if (matrix.missing[rowsmiss[i], j] == 0) {
            length.gap <- length.gap + 1
          } else {
            length.gap <- 0
          }
        }
      }
      nmis <- sum(matrix.missing[rowsmiss[i], ] == 0)
    }
  }
  if (inherits(data, "stslist")) {
    data[matrix.missing == 0] <- attr(data, "nr")
  } else {
    data[matrix.missing == 0] <- NA
  }
  return(data)
}
