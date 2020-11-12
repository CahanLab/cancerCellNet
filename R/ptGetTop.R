#' @title
#' Find the best gene pairs for training
#' @description
#' Find the gene pairs that most distinguish a cancer group from the rest
#'
#' @param expDat expDat
#' @param cell_labels named vector, value is grp, name is cell name
#' @param cgenes_list the list of labelled cgenes
#' @param topX number of genepairs for training
#' @param sliceSize the size of the slice for pair transform. Default at 5e3
#' @param quickPairs TRUE if wanting to select the gene pairs in a quick fashion
#' @param coreProportion the proportion of logical cores for finding top scoring gene pairs (not applicable for quick pairs)

#' @import parallel
#' @return vector of top gene-pair names
#'
#' @export
ptGetTop<-function(expDat, cell_labels, cgenes_list=NA, topX=50, sliceSize = 5e3, quickPairs = FALSE, coreProportion = 1/4){

  ncores = parallel::detectCores(logical = TRUE) # detect the number of cores in the system
  mcCores = 1
  if(ncores * coreProportion > 1){
    mcCores = round(ncores * coreProportion)
  }

  if(!quickPairs){
    #ans<-vector()
    ans = list()
    genes = rownames(expDat)

    cat(ncores, "logical cores in total", "  --> ", mcCores, "cores running to find top scoring gene pairs...","\n")

    # make a data frame of pairs of genes that will be sliced later
    pairTab = makePairTab(genes)

    if(topX > nrow(pairTab)) {
      stop(paste0("The data set has ", nrow(pairTab), " total combination of gene pairs. Please select a smaller topX."))
    }

    # setup tmp ans list of sc_testPattern
    cat("setup ans and make pattern\n")
    grps = unique(cell_labels)
    myPatternG = sc_sampR_to_pattern(as.character(cell_labels))
    statList = list()
    for(grp in grps){
      statList[[grp]] = data.frame()
    }

    # make the pairedDat, and run sc_testPattern
    cat("make pairDat on slice and test\n")
    nPairs = nrow(pairTab)
    cat("nPairs = ",nPairs,"\n")
    str = 1
    stp = min(c(sliceSize, nPairs)) # detect what is smaller the slice size or npairs


    if(mcCores == 1) {
      while(str <= nPairs){
        if(stp>nPairs){
          stp = nPairs
        }
        tmpTab = pairTab[str:stp,]
        tmpPdat = ptSmall(expDat, tmpTab)

        tmpAns = lapply(X = myPatternG, FUN = sc_testPattern, expDat=tmpPdat)
        cat(str,"-", stp,"\n")

        for(gi in seq(length(myPatternG))){
          grp = grps[[gi]]
          statList[[grp]] = rbind( statList[[grp]],  tmpAns[[grp]])
        }


        str = stp+1
        stp = str + sliceSize - 1

      }
    }
     else {
       cl = snow::makeCluster(mcCores, type="SOCK")
       while(str <= nPairs){
         if(stp>nPairs){
           stp = nPairs
         }
         tmpTab = pairTab[str:stp,]
         tmpPdat = ptSmall(expDat, tmpTab)

         tmpAns = snow::parLapply(cl = cl, x = myPatternG, fun = sc_testPattern, expDat=tmpPdat)
         cat(str,"-", stp,"\n")

         for(gi in seq(length(myPatternG))){
           grp = grps[[gi]]
           statList[[grp]] = rbind( statList[[grp]],  tmpAns[[grp]])
         }


         str = stp+1
         stp = str + sliceSize - 1

       }
       snow::stopCluster(cl)
     }


    cat("compile results\n")
    for(grp in grps){
      tmpAns = findBestPairs(statList[[grp]], topX)
      ans[[grp]] = tmpAns
    }
    return(ans)
  }
  else {

    myPatternG = sc_sampR_to_pattern(as.character(cell_labels))
    ans = list()
    for(cct in names(cgenes_list)){
      genes = cgenes_list[[cct]]
      pairTab = makePairTab(genes)

      nPairs = nrow(pairTab)
      cat("nPairs = ", nPairs," for ", cct, "\n")

      tmpPdat = ptSmall(expDat, pairTab)

      tmpAns = findBestPairs(sc_testPattern(myPatternG[[cct]], expDat=tmpPdat), topX)

      ans[[cct]] = tmpAns
    }

    return(ans)
  }
}

#' @title
#' function for parallel computation of quick pairs (non-windows...)
#' @description
#' this function was written to compute the gene pairs in parallel for quick pairs
#' @param inputPackage a list of 2 things: template vectors and gene pair matrix
#' @param topX number of top pairs selection
#' @return top scoring gene pairs
parallel_quickPairs <- function(inputPackage, topX) {

  ans = findBestPairs(sc_testPattern(inputPackage[[1]], expDat=inputPackage[[2]]), topX)
  return(ans)
}


#' @title
#' Make the pair tabs
#' @description
#' Generate all the combination of gene pairs
#' @param genes a vector of all the genes in the expression matrix
#' @return a gene pair table
makePairTab <- function(genes) {
  pTab = t(combn(genes, 2))
  colnames(pTab) = c("genes1", "genes2")
  pTab = cbind(pTab, pairName=paste(pTab[,1], "_",pTab[,2], sep=''))

  return(pTab)
}



#' @title
#' Pair Transform on small scale
#' @description
#' Performs gene pair comparison on a smaller subset to conserve RAM
#' @param expDat the gene expression dataframe
#' @param pTab the gene pair table generated as one of the intermediate step from \code{\link{ptGetTop}}
#' @return a dataframe with gene pairs as rows and samples as columns
ptSmall <- function(expDat, pTab) {
  npairs = nrow(pTab)
  ans = matrix(0, nrow=npairs, ncol=ncol(expDat))
  genes1 = as.vector(pTab[, "genes1"])
  genes2 = as.vector(pTab[, "genes2"])

  for(i in seq(nrow(pTab))) {
    ans[i,] = as.numeric(expDat[genes1[i],]>expDat[genes2[i],])
  }
  colnames(ans) = colnames(expDat)
  rownames(ans) = as.vector(pTab[, "pairName"])

  return(ans)
}

#' @title
#' Find best pairs
#' @description
#' Perform finding the best and diverse set of gene pairs for training
#' @param xdiff statList of a certain group generated as an intermediate step from \code{\link{ptGetTop}}
#' @param n the number of top pairs
#' @param maxPer indicates the maximum number of pairs that a gene is allowed to be in
#' @return vector of suitable gene pairs
findBestPairs <- function(xdiff, n=50,maxPer=3) {

  # error catching in case the number of pairs wanted is more than pairs generated
  if(nrow(xdiff) < n) {
    cat("there are only", nrow(xdiff), "genepairs generated.", "\n")

    ans = as.vector(xdiff$cval)
    names(ans) = rownames(xdiff)
    return(ans)

  }
  else {
    xdiff = xdiff[order(abs(xdiff$cval), decreasing=TRUE),]

    genes = unique(unlist(strsplit(rownames(xdiff), "_")))
    countList = rep(0, length(genes))
    names(countList) = genes

    i = 0
    ans_names = vector()
    ans_signs = vector()

    xdiff_index = 1
    pair_names = rownames(xdiff)

    backup_vector = c()
    backup_vector_sign = c()

    while(i < n ){
      tmpAns = pair_names[xdiff_index]
      tmpSigns = sign(as.numeric(xdiff[tmpAns, "cval"]))

      tgp = unlist(strsplit(tmpAns, "_"))

      if((countList[tgp[1]] < maxPer) & (countList[tgp[2]] < maxPer )){

        # record down the gene pair name and sign
        ans_names = append(ans_names, tmpAns)
        ans_signs = append(ans_signs, tmpSigns)

        countList[ tgp[1] ] = countList[ tgp[1] ]+ 1
        countList[ tgp[2] ] = countList[ tgp[2] ]+ 1

        i = i+1
      }

      else {
        backup_vector = c(backup_vector, tmpAns) # place into backup vector
        backup_vector_sign = c(backup_vector_sign, tmpSigns)
      }


      xdiff_index = xdiff_index + 1

      # in the case where the original list is exhausted, dig into the backup vector
      if(xdiff_index > length(pair_names)) {
        additional_pairs = backup_vector[1:(n - i)]
        additional_pairs_sign = backup_vector_sign[1:(n - i)]

        ans_names = c(ans_names, additional_pairs)
        ans_signs = c(ans_signs, additional_pairs_sign)
        i = length(ans_signs)
      }

    }

    # assign the signs and names to the return answer
    ans = ans_signs
    names(ans) = ans_names

    ans = na.omit(ans) # just in case there were NA
    return(ans)
  }

}
