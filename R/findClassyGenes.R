#' @title
#' Find Classifier-Worthy Gene Candidates
#' @description
#' Find classifier-worthy genes for each cancer category to train the classifier
#'
#' @param expDat a matrix of normalized expression data from \code{\link{trans_prop}}
#' @param sampTab a dataframe of the sample table
#' @param dLevel a string indicating the column name in sample table that contains the cancer category
#' @param topX an integer indicating the number of top positive classification genes for each category to select for training. Will also select topX number of negative classification genes.
#' @param dThresh a number representing the detection threshold
#' @param alpha1 a number representing proportion of cells in which a gene must be considered detected (as defined in geneStats)
#' @param alpha2 a number representing lower proportion of cells for genes that must have higher expression level
#' @param mu a number represeting threshold for average expression level of genes passing the lower proportion criteria
#'
#' @return a list containing two lists: a list of classifier worthy genes named 'cgenes' and a list of cancer category named 'grps'
#' @export
findClassyGenes<-function(expDat, sampTab, dLevel, topX=25, dThresh=0, alpha1=0.05, alpha2=.001, mu=2, sliceSize = 1000) {
  if((dLevel %in% colnames(sampTab)) == FALSE) {
    stop("Please enter the correct column name for sampTab that indicates the categories")
  }

  if((topX * 3 > nrow(expDat)) == TRUE) {
    stop(paste0("Please enter a topX value smaller than ", as.integer(nrow(expDat) / 3)))
  }

  # remove duplicated genes
  expDat = expDat[!duplicated(rownames(expDat)), ]

  gsTrain<-sc_statTab(expDat, dThresh=dThresh)
  ggenes<-sc_filterGenes(gsTrain, alpha1=alpha1, alpha2=alpha2, mu=mu)
  grps<-as.vector(sampTab[,dLevel])
  names(grps)<-rownames(sampTab)

  # more robust way to avoid subscript out of bound
  if(length(ggenes) != nrow(expDat)) {
    expDat = expDat[ggenes,]
  }

  ncores<-parallel::detectCores(logical = FALSE) # detect the number of cores in the system
  mcCores<-1
  if(ncores/2>1){
    mcCores<- round(ncores /2)
  }

  cat(ncores, "cores in total", "  --> ", mcCores, "cores running to find classification genes...","\n")
  xdiff<-gnrAll(expDat, grps, sliceSize)

  if (Sys.info()[['sysname']] == "Windows") {
    cl<-snow::makeCluster(mcCores, type="SOCK")
    cgenes<-snow::parLapply(cl = cl, x = xdiff, fun = getClassGenes, topX = topX)
    stopCluster(cl)
  }
  else {
    cgenes<-parallel::mclapply(xdiff, getClassGenes, topX=topX, mc.cores=mcCores)
  }
  labelled_cgenes <- cgenes
  cgenes<-unique(unlist(cgenes))
  list(cgenes=cgenes, grps=grps, labelled_cgenes=labelled_cgenes)
}

#' @title
#' Find Genes that Pass Criteria
#' @description
#' Based on idea that reliably detected genes will either be detected in many cells, or highly expressed in a small number of cells or both
#' @param geneStats a matrix containing stats of genes generated from running \code{\link{sc_statTab}}
#' @param alpha1 a number representing proportion of cells in which a gene must be considered detected (as defined in geneStats)
#' @param alpha2 a number representing lower proportion of cells for genes that must have higher expression level
#' @param mu a number represeting threshold for average expression level of genes passing the lower proportion criteria
#' @return a vector of gene symbols
#' @export
sc_filterGenes<-function(geneStats, alpha1=0.1, alpha2=0.01, mu=2){

  #return
  rownames(geneStats[(geneStats$alpha > alpha1) | (geneStats$alpha > alpha2 & geneStats$mu > mu), ])
}

#' @title
#' Find Genes in Clusters compare to other categories
#' @description
#' Find genes higher in a cluster compared to all other cells
#' @param expDat a matrix of normalized gene expression taken from \code{\link{trans_prop}}
#' @param cellLabels a vector of all the named vector of cancer categories
#' @return list of dataFrames containing pval, cval and holm for each gene in each cancer category
gnrAll<-function(expDat, cellLabels, sliceSize){
  ncores<-parallel::detectCores(logical = FALSE) # detect the number of cores in the system
  mcCores<-1
  if(ncores/2>1){
    mcCores <- round(ncores / 2)
  }

  nGenes = nrow(expDat)
  cat("nGenes = ",nGenes,"\n")
  str = 1
  stp = min(c(sliceSize, nGenes)) # detect what is smaller the slice size or nGenes

  myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
  statList<-list()
  grps<-unique(cellLabels)

  for(grp in grps){
    statList[[grp]]<-data.frame()
  }

  if (Sys.info()[['sysname']] == "Windows") {
    cl<-snow::makeCluster(mcCores, type="SOCK")

    while(str <= nGenes){
      if(stp>nGenes){
        stp <- nGenes
      }
      cat(str,"-", stp,"\n")

      tempExpDat = expDat[str:stp, ]
      tmpAns<-snow::parLapply(cl = cl, x = myPatternG, fun = sc_testPattern, expDat=tempExpDat)

      for(gi in seq(length(myPatternG))){
        grp<-grps[[gi]]
        statList[[grp]]<-rbind( statList[[grp]],  tmpAns[[grp]])
      }

      str<-stp+1
      stp<-str + sliceSize - 1

    }
    stopCluster(cl)

  }
  else {
    while(str <= nGenes){
      if(stp>nGenes){
        stp <- nGenes
      }
      cat(str,"-", stp,"\n")

      tempExpDat = expDat[str:stp, ]
      tmpAns<-parallel::mclapply(myPatternG, sc_testPattern, expDat=tempExpDat, mc.cores=mcCores) # this code cannot run on windows
      for(gi in seq(length(myPatternG))){
        grp<-grps[[gi]]
        statList[[grp]]<-rbind( statList[[grp]],  tmpAns[[grp]])
      }

      str<-stp+1
      stp<-str + sliceSize - 1
    }
  }

  cat("Done testing\n")

  # return
  return(statList)
}

#' @title
#' Find Classy Genes
#' @description
#' Extract genes suitable for training classifier
#' @param diffRes a dataframe with pval, cval, holm, and rownames as the gene names
#' @param topX a number dicataing the number of genes to select for training classifier
#' @param bottom logic if true use the top x genes with - cvals
#'
#' @return a vector of genes that are good for training classifier for that category
getClassGenes<-function(diffRes, topX=25, bottom=TRUE, leastDiff = TRUE) {
  #exclude NAs
  xi<-which(!is.na(diffRes$cval))
  diffRes<-diffRes[xi,] # exclude the NAs. Select the rows that does not have NAs

  diffRes<-diffRes[order(diffRes$cval, decreasing=TRUE),] #order based on classification value largest to lowest
  ans<-rownames(diffRes[1:topX,]) # get the top 20 genes

  if(bottom){
    ans<-append(ans, rownames( diffRes[nrow(diffRes) - ((topX-1):0),]))
  }

  # get the least differentially expressed genes as house holders
  if(leastDiff) {
    sameRes<-diffRes[order(abs(diffRes$cval), decreasing = FALSE), ]
    ans<-append(ans, rownames(sameRes[1:topX, ]))
  }

  #return
  ans
}
