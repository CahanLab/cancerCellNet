#' Randomize the data matrix
#'
#' To generate random profiles by providing random permutations of the data matrix
#' @param expDat the data matrix
#' @param num number of random profiles
#' @return a randomized matrix
#' @export
randomize<-function(expDat, num=50) {

  randDat = t(apply(expDat, 1, sample))
  randDat = apply(randDat, 2, sample)

  randDat = randDat[,sample(1:ncol(randDat), num)]

  colnames(randDat) = paste0(rep("rand_", num), 1:num)
  rownames(randDat) = rownames(expDat)

  return(randDat)
}
