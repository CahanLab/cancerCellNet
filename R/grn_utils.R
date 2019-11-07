# utility functions

#' @title logRank gene expression
#' @description rank the gene expression values and log10 the ranks
#' @param expTrain the expression matrix
#' @return the log10 rank of genes
#' @export
logRank <- function(expTrain) {
  expTrain = apply(expTrain, FUN = rank, MARGIN = 2)

  return(log10(expTrain))
}

#' @title process classifier importance
#' @description process the importance scores of gene pairs into enriched genes
#' @param classifier the Random Forest classifier
#' @param xpairs the xpair_list from training
#' @param prune blooln indicate whether to select genes exclusive to cancer category
#'
#' @return a list of genes with their importance divided by groups
#' @exports
processImportance <- function(classifier, xpairs, prune = TRUE) {
  genePairImportance = classifier$importance

  gene_importanceList = list()

  ignoreList = vector()
  # loop through individual grousp
  for(cancerGroup in names(xpairs)) {
    temp_pairList = xpairs[[cancerGroup]]

    group_geneImportance = vector()
    # loop through individual gene pair
    for(tempPair in names(temp_pairList)){
      if(temp_pairList[[tempPair]] == 1) {
        # assign the first gene with the importance of the the genepair
        group_geneImportance[[strsplit(tempPair, split = "_")[[1]][1]]] = genePairImportance[tempPair, 1]
      }
      else {
        # if the value is -1, assign the second gene with the importance of the genepair
        group_geneImportance[[strsplit(tempPair, split = "_")[[1]][2]]] = genePairImportance[tempPair, 1]
      }
    }

    gene_importanceList[[cancerGroup]] = group_geneImportance
    ignoreList = c(ignoreList, names(group_geneImportance))
  }

  ignoreList = ignoreList[duplicated(ignoreList)]
  if(prune == TRUE) {
    for(cancer_group in names(gene_importanceList)) {
      tempGene = gene_importanceList[[cancer_group]]
      gene_importanceList[[cancer_group]] = tempGene[!(tempGene %in% ignoreList)]
    }
  }

  return(gene_importanceList)
}
