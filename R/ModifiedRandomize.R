#' Randomize the data matrix
#'
#' To generate random profiles by providing random permutations of the data matrix
#' @param expDat the data matrix
#' @param num number of random profiles
#' @return a randomized matrix
ModifiedRandomize<-function(expDat, num=50) {

  if (num <= ncol(expDat)) {
    randDat = t(apply(expDat, 1, sample))
    randDat = apply(randDat, 2, sample)

    randDat = randDat[,sample(1:ncol(randDat), num)]

    colnames(randDat) = paste0(rep("rand_", num), 1:num)
    rownames(randDat) = rownames(expDat)
    randDat

    returnRandDat = randDat
  }
  else {
    returnRandDat = NULL
    while (num > ncol(expDat)) {
      randDat = t(apply(expDat, 1, sample))
      randDat = apply(randDat, 2, sample)

      randDat = randDat[, sample(1:ncol(randDat), ncol(randDat))]


      if (is.null(returnRandDat) == TRUE) {
        returnRandDat = randDat
      }
      else {
        returnRandDat = cbind(returnRandDat, randDat)
      }
      num = num - ncol(expDat)
    }

    randDat = t(apply(expDat, 1, sample))
    randDat = apply(randDat, 2, sample)

    randDat = randDat[, sample(1:ncol(randDat), num)]
    returnRandDat = cbind(returnRandDat, randDat)

    colnames(returnRandDat) = paste0(rep("rand_", ncol(returnRandDat)), 1:ncol(returnRandDat))
    rownames(returnRandDat) = rownames(expDat)
  }

  return(returnRandDat)
}
