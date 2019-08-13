#' @title
#' Make Gene Comparison Table
#' @description
#' To compile an expression table for comparison
#'
#' @param queryExpTab a matrix of normalized expression query data from \code{\link{trans_prop}}
#' @param avgGeneTab a matrix of averaged expression of training data from \code{\link{avgGeneCat}}.
#' @param querySample a vector or string indicating the query samples
#' @param trainingCat a vector or string indicating the categories of the training data
#' @param geneSamples a vector or string indicating the genes of interest
#'
#' @return a matrix that combines query and training data with genes of interest for comparison
#' @export

makeGeneCompareTab<-function(queryExpTab, avgGeneTab, querySample = NULL, trainingCat=NULL, geneSamples) {

  if(is.null(querySample) == TRUE) {
    filteredQuery = queryExpTab[geneSamples, ]
    querySample = colnames(filteredQuery)
  }
  else if(all(querySample %in% colnames(queryExpTab)) == FALSE) {
    filteredQuery = queryExpTab[geneSamples, ]
    querySample = colnames(filteredQuery)
    print("Please enter sample names that are in the query table")

  }
  else {
    filteredQuery = queryExpTab[geneSamples, querySample]
  }

  if(is.null(trainingCat) == TRUE) {
    filteredGeneTab = avgGeneTab[geneSamples, ]
    trainingCat = colnames(filteredGeneTab)
  }
  else if(all(trainingCat %in% colnames(avgGeneTab)) == FALSE) { #maybe modified it later
    filteredGeneTab = avgGeneTab[geneSamples, ]
    trainingCat = colnames(filteredGeneTab)

    print("Please enter proper category names that are in the filtered gene table")

  }
  else {
    filteredGeneTab = avgGeneTab[geneSamples, trainingCat]
  }

  if(all(rownames(filteredQuery) == rownames(filteredGeneTab)) == TRUE) {
    print("All Good")
  }

  returnLabel = c(querySample, trainingCat)
  returnMatrix = cbind(filteredQuery, filteredGeneTab)
  colnames(returnMatrix) = returnLabel

  #return
  returnMatrix
}
