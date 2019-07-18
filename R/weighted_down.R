#' @title
#' Down Sample Gene Expressions
#'
#' @description
#' Weighted subtraction from mapped reads and simulate expression profile of  _total_ mapped reads
#' for normalization.
#'
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#' @param dThresh the threshold at which anything lower than that is 0
#'
#' @return matrix of downsampled read mapped to genes/transcripts
#'
#' @export
weighted_down<-function(expRaw, total=1e5, dThresh=0) {
  expCountDnW<-apply(expRaw, 2, downSampleW, total=total, dThresh=dThresh)
  expCountDnW
}

#' Weighted subtraction from mapped reades
#'
#' Simulate expression profile of  _total_ mapped reads as a way to normalize the counts
#'
#' @param vector a vector of total mapped reads per gene/transcript
#' @param total the sum of the read counts after transformation
#' @param dThresh the threshold at which anything lower is 0 after transformation. Usually 0
#'
#' @return vector of downsampled read mapped to genes/transcripts
downSampleW<-function(vector, total=1e5, dThresh=0){

  totalSignal<-sum(vector) # get the sum of the vector
  wAve<-vector/totalSignal #

  resid<-totalSignal-total #num to subtract from sample
  residW<-wAve*resid # amount to substract from each gene

  ans<-vector-residW
  ans[which(ans<dThresh)]<-0
  ans
}
