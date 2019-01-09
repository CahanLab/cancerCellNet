#' @title
#' Make Gene Comparison Table
#' @description
#' To compile an expression table for comparison
#'
#' @param queryExpTab a matrix of normalized expression query data from \code{\link{trans_prop}}
#' @param avgGeneTab a matrix of averaged expression of training data from \code{\link{avgGeneCat}}
#' @param querySample a vector or string indicating the query samples
#' @param trainingCat a vector or string indicating the categories of the training data
#' @param geneSamples a vector or string indicating the genes of interest
#'
#' @return a matrix that combines query and training data with genes of interest for comparison
#' @export

makeGeneCompareTab<-function(queryExpTab, avgGeneTab, querySample = NULL, trainingCat=NULL, geneSamples) {

  if(is.null(querySample) == TRUE) {
    filteredQuery = queryExpTab[geneSamples, ]

  }
  else if(all(querySample %in% colnames(queryExpTab)) == FALSE) {
    filteredQuery = queryExpTab[geneSamples, ]
  }
  else {
    filteredQuery = queryExpTab[geneSamples, querySample]
  }

  #
  if(is.null(trainingCat) == TRUE) {
    filteredGeneTab = avgGeneTab[geneSamples, ]
  }
  else if(all(trainingCat %in% colnames(avgGeneTab)) == FALSE) { #maybe modified it later
    filteredGeneTab = avgGeneTab[geneSamples, ]
  }
  else {
    filteredGeneTab = avgGeneTab[geneSamples, trainingCat]
  }

  if(all(rownames(filteredQuery) == rownames(filteredGeneTab)) == TRUE) {
    print("All Good")
  }
  returnMatrix = cbind(filteredQuery, filteredGeneTab)
  #return
  returnMatrix
}
