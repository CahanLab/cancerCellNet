#' find best pairs
#'
#' find best pairs
#'
#' @param expDat expDat
#' @param cellLabels named vector of cell groups
#'
#' @return vector of pairs
#'
#' @export
gnrBP<-function(expDat, cellLabels,topX=50){

  ans<-vector()

  myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
  for(i in seq(length(myPatternG))){
    cat(i,"\n")
    xres<-sc_testPattern(myPatternG[[i]], expDat=expDat)
    tmpAns<-findBestPairs(xres, topX)
    ans<-append(ans, tmpAns)
  }
  unique(ans)
}
