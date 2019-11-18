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
#' @export
processImportance <- function(classifier, xpairs, prune = TRUE) {
  genePairImportance = classifier$importance

  gene_importanceList = list()

  ignoreList = vector()
  # loop through individual grousp
  for(cancerGroup in names(xpairs)) {
    temp_pairList = names(xpairs[[cancerGroup]])

    geneNames = unique(unlist(strsplit(x = temp_pairList, split = "_"))) # get unique genes from gene pairs
    group_geneImportance = rep(0, length(geneNames))
    names(group_geneImportance) = geneNames

    # loop through individual gene pair
    for(tempPair in temp_pairList){
      gene1 = strsplit(tempPair, split = "_")[[1]][1]
      gene2 = strsplit(tempPair, split = "_")[[1]][2]

      tempImportance = genePairImportance[tempPair, 1]
      if(group_geneImportance[[gene1]] < tempImportance){
        group_geneImportance[[gene1]] = tempImportance
      }

      if(group_geneImportance[[gene2]] < tempImportance){
        group_geneImportance[[gene2]] = tempImportance
      }

    }

    gene_importanceList[[cancerGroup]] = group_geneImportance
    ignoreList = c(ignoreList, names(group_geneImportance))
  }

  # find the duplicated genes that are across different
  ignoreList = ignoreList[duplicated(ignoreList)]
  if(prune == TRUE) {
    for(cancer_group in names(gene_importanceList)) {
      tempGene = gene_importanceList[[cancer_group]]
      gene_importanceList[[cancer_group]] = tempGene[!(tempGene %in% ignoreList)]
    }
  }

  return(gene_importanceList)
}
