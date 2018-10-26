#' makes vector of gene pairs, iterates over this and computes pairDat, sc_testPattern, then, at the end, findBestPairs
#'
#' @param expDat expDat
#' @param cell_labels named vector, value is grp, name is cell name
#' @param topX 50
#'
#' @return vector of gene-pair names
#'
#' @export
ptGetTop<-function(expDat, cell_labels, topX=50, sliceSize = 5e3){
  ans<-vector()

  # make a data frame of pairs of genes that will be sliced later
  cat("Making pairTable\n")
  ngenes<-nrow(expDat)
  genes<-rownames(expDat)
  genes1<-vector()
  genes2<-vector()
  for(i in 1:ngenes){
    for(j in 1:ngenes){
      if(j>i){
        genes1<-append(genes1, genes[i])
        genes2<-append(genes2, genes[j])
      }
    }
  }

  pairTab = data.frame(genes1=genes1, genes2=genes2)
  pairNames<-paste(pairTab[,1], "_",pairTab[,2], sep='')
  pairTab <- cbind(pairTab, pairName=pairNames)

  ###
  # setup tmp ans list of sc_testPattern
  cat("setup ans and make pattern\n")
  grps<-unique(cell_labels)
  myPatternG<-sc_sampR_to_pattern(as.character(cell_labels))
  statList<-list()
  for(grp in grps){
    statList[[grp]]<-data.frame()
  }

  # make the pairedDat, and run sc_testPattern
  cat("make pairDat on slice and test\n")
  nPairs = nrow(pairTab)
  cat("nPairs = ",nPairs,"\n")
  str = 1
  stp = min(c(sliceSize, nPairs))
  while(str <= nPairs){
    if(stp>nPairs){
      stp <- nPairs
    }
    cat(str,"-", stp,"\n")
    tmpTab<-pairTab[str:stp,]
    tmpPdat<-ptSmall(expDat, tmpTab)

    for(gi in seq(length(myPatternG))){

      grp<-grps[[gi]]
      statList[[grp]]<-rbind( statList[[grp]], sc_testPattern(myPatternG[[gi]], expDat=tmpPdat) )
    }


    str = stp+1
    stp = str + sliceSize - 1
  }
  cat("compile results\n")
  for(grp in grps){
    tmpAns<-findBestPairs(statList[[grp]], topX)
    ans<-append(ans, tmpAns)
  }
  unique(ans)
}

#'
#'
#'
ptSmall<-function
(expDat,
 pTab){
  npairs = nrow(pTab)
  ans<-matrix(0, nrow=npairs, ncol=ncol(expDat))
  genes1<-as.vector(pTab$genes1)
  genes2<-as.vector(pTab$genes2)

  for(i in seq(nrow(pTab))){
    #cat(genes1[i], ": ", genes2[i],"\n")
    ans[i,]<-as.numeric(expDat[genes1[i],]>expDat[genes2[i],])
  }
  colnames(ans)<-colnames(expDat)
  rownames(ans)<-as.vector(pTab$pairName)
  ans
}

#' finds the best pairs to use
#'
#' @param xdiff xdiff
#' @param n number of pairs
#' @param maxPer indicates the number of pairs that a gene is allowed to be in
#'
#' @return vector of good pairs
#'
#' @export
findBestPairs<-function # find best and diverse set of pairs
(xdiff,
 n=50,
 maxPer=3){

  xdiff<-xdiff[order(xdiff$cval, decreasing=TRUE),]
  genes<-unique(unlist(strsplit(rownames(xdiff), "_")))
  countList <- rep(0, length(genes))
  names(countList) <- genes

  i<-0
  ans<-vector()
  xdiff_index<-1
  pair_names<-rownames(xdiff)
  while(i < n ){
    tmpAns<-pair_names[xdiff_index]
    tgp <- unlist(strsplit(tmpAns, "_"))
    if( (countList[ tgp[1] ] < maxPer) & (countList[ tgp[2] ] < maxPer )){
      ans<-append(ans, tmpAns)
      countList[ tgp[1] ] <- countList[ tgp[1] ]+ 1
      countList[ tgp[2] ] <- countList[ tgp[2] ]+ 1
      i<-i+1
    }
    xdiff_index <- xdiff_index + 1
  }
  ans
}
