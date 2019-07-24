#' @title
#' Find the best gene pairs for training
#' @description
#' Find the gene pairs that most distinguish a cancer group from the rest
#'
#' @param expDat raw expression data
#' @param cell_labels named vector, value is grp, name is cell name
#' @param topX number of
#' @param sliceSize
#'
#' @return vector of gene-pair names
#'
#' @export
ptGetTop_old<-function(expDat, cell_labels, topX=50, sliceSize = 5e3){
  ans<-vector()

  # make a data frame of pairs of genes that will be sliced later
  cat("Making pairTable\n")
  ngenes<-nrow(expDat)
  genes<-rownames(expDat)

  genes1<-vector()
  genes2<-vector()

  #there is probably a faster way to do this since this is O(n^2)
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
  myPatternG<-sc_sampR_to_pattern(as.character(cell_labels)) #map each cancer categories into 0 and 1
  statList<-list()

  for(grp in grps){
    statList[[grp]]<-data.frame()
  }

  # make the pairedDat, and run sc_testPattern
  cat("make pairDat on slice and test\n")
  nPairs = nrow(pairTab)
  cat("nPairs = ",nPairs,"\n")
  str = 1
  stp = min(c(sliceSize, nPairs)) #get the smaller one of number of pairs vs slice size

  while(str <= nPairs){

    #if the slice size is greater than nPairs due to the error from user's end
    if(stp>nPairs){
      stp <- nPairs
    }

    cat(str,"-", stp,"\n")
    tmpTab<-pairTab[str:stp,]

    tmpPdat<-ptSmall(expDat, tmpTab)

    # linearly fit the expDat with Patterns to find good gene pairs
    # add those into the correspnding groups for the subsetted tempPdat
    for(gi in seq(length(myPatternG))){
      grp<-grps[[gi]]
      statList[[grp]]<-rbind( statList[[grp]], sc_testPattern(myPatternG[[gi]], expDat=tmpPdat) )
    }


    str = stp+1 #current starting position
    stp = str + sliceSize - 1 #the ending position
  }

  cat("compile results\n")

  for(grp in grps){
    tmpAns<-findBestPairs(statList[[grp]], topX)
    ans<-append(ans, tmpAns)
  }
  unique(ans)
}

#' @title
#' Find the best gene pairs for training (new code)
#' @description
#' Find the gene pairs that most distinguish a cancer group from the rest
#'
#' @param expDat expDat
#' @param cell_labels named vector, value is grp, name is cell name
#' @param cgenes_list the list of labelled cgenes 
#' @param topX number of genepairs for training 
#' @param sliceSize the size of the slice. Default at 5e3
#' @param quickPairs TRUE if wanting to select the gene pairs in a quick fashion 
#'
#' @import parallel
#' @return vector of gene-pair names
#'
#' @export
ptGetTop <-function(expDat, cell_labels, cgenes_list=NA, topX=50, sliceSize = 5e3, quickPairs = FALSE ){
  if(!quickPairs){
    ans<-vector()
    genes<-rownames(expDat)

    ncores <- parallel::detectCores() # detect the number of cores in the system
    mcCores <- 1
    if(ncores>1){
      mcCores <- ncores - 1
    }
    cat(ncores, "  --> ", mcCores,"\n")

    # make a data frame of pairs of genes that will be sliced later
    pairTab<-makePairTab(genes)

    if(topX > nrow(pairTab)) {
      stop(paste0("The data set has ", nrow(pairTab), " total combination of gene pairs. Please select a smaller topX."))
    }

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
    stp = min(c(sliceSize, nPairs)) # detect what is smaller the slice size or npairs

    while(str <= nPairs){
      if(stp>nPairs){
        stp <- nPairs
      }
      cat(str,"-", stp,"\n")
      tmpTab<-pairTab[str:stp,]
      tmpPdat<-ptSmall(expDat, tmpTab)

      ### new

      if (Sys.info()[['sysname']] == "Windows") {
        tmpAns<-lapply(myPatternG, sc_testPattern, expDat=tmpPdat)
      }
      else {
        tmpAns<-parallel::mclapply(myPatternG, sc_testPattern, expDat=tmpPdat, mc.cores=mcCores) # this code cannot run on windows
      }
      ### names(tmpAns) <- grps

      ###for(gi in seq(length(myPatternG))){
      ###	grp<-grps[[gi]]
      ###	statList[[grp]]<-rbind( statList[[grp]], sc_testPattern(myPatternG[[gi]], expDat=tmpPdat) )
      ### }

      for(gi in seq(length(myPatternG))){
        grp<-grps[[gi]]
        #cat(i, " grp: ",grp,"\n")
        statList[[grp]]<-rbind( statList[[grp]],  tmpAns[[grp]])
      }


      str = stp+1
      stp = str + sliceSize - 1
    }

    cat("compile results\n")
    for(grp in grps){
      tmpAns<-findBestPairs(statList[[grp]], topX)
      ans<-append(ans, tmpAns)
    }
    return(unique(ans))

  }else{
    myPatternG<-sc_sampR_to_pattern(as.character(cell_labels))
    ans<-vector()

    for(cct in names(cgenes_list)){
      genes<-cgenes_list[[cct]]
      pairTab<-makePairTab(genes)

      nPairs = nrow(pairTab)
      cat("nPairs = ", nPairs," for ", cct, "\n")

      tmpPdat<-ptSmall(expDat, pairTab)

      tmpAns<-findBestPairs( sc_testPattern(myPatternG[[cct]], expDat=tmpPdat), topX)
      ans<-append(ans, tmpAns)
    }

    return(unique(ans))
  }
}

#' @title
#' Make the pair tabs
#' @description
#' Generate all the combination of gene pairs
#'
#' @param genes a vector of all the genes in the expression matrix
makePairTab<-function(genes){
  pTab<-t(combn(genes, 2))
  colnames(pTab)<-c("genes1", "genes2")
  pTab<-cbind(pTab, pairName=paste(pTab[,1], "_",pTab[,2], sep=''))
  pTab
}



#' @title
#' Pair Transform on small scale
#' @description
#' Performs gene pair comparison on a smaller subset to conserve RAM
#'
#' @param expDat the gene expression dataframe
#' @param pTab the gene pair table generated as one of the intermediate step from \code{\link{ptGetTop}}
#'
#' @return a dataframe with gene pairs as rows and samples as columns
ptSmall<-function(expDat, pTab){
  npairs = nrow(pTab)
  ans<-matrix(0, nrow=npairs, ncol=ncol(expDat))
  genes1<-as.vector(pTab[, "genes1"])
  genes2<-as.vector(pTab[, "genes2"])

  for(i in seq(nrow(pTab))){
    #cat(genes1[i], ": ", genes2[i],"\n")
    ans[i,]<-as.numeric(expDat[genes1[i],]>expDat[genes2[i],])
  }
  colnames(ans)<-colnames(expDat)
  rownames(ans)<-as.vector(pTab[, "pairName"])
  ans
}

#' @title
#' Find best pairs
#' @description
#' Perform finding the best and diverse set of gene pairs for training
#'
#' @param xdiff statList of a certain group generated as an intermediate step from \code{\link{ptGetTop}}
#' @param n the number of top pairs
#' @param maxPer indicates the maximum number of pairs that a gene is allowed to be in
#'
#' @return vector of suitable gene pairs
findBestPairs<-function(xdiff, n=50,maxPer=3){

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

  #return
  ans
}



