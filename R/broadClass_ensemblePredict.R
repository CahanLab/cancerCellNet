#' @title
#' Predict query using an ensemble of broad class classifier
#' @description
#' The function predicts the query data using an ensemble broad class classifier
#'
#' @param broadReturnList a list of broadReturn from training
#' @param queryData the expression data of query data
#' @param nrand the number of random profiles generate for evaluation process
#' @return an average classification score matrix
#'
#' @export
broadClass_ensemblePredict <- function(broadReturnList, queryData, nrand = 2) {

  list_classMatrix = list()
  for(broadClassName in names(broadReturnList)) {
    broadReturn = broadReturnList[[broadClassName]]
    broadCnProc = broadReturn$cnProc

    tempClassMatrix = broadClass_predict(cnProc = broadCnProc, expDat = queryData, nrand = nrand)
    list_classMatrix[[broadClassName]] = tempClassMatrix
  }

  catOrder = NULL
  avgClassMatrix = NULL
  for(classMatrix in list_classMatrix) {

    if(is.null(catOrder) == TRUE) {
      catOrder = rownames(classMatrix)
    }

    if(is.null(avgClassMatrix) == TRUE) {
      avgClassMatrix = classMatrix[catOrder, ]
    }
    else {
      avgClassMatrix = avgClassMatrix[catOrder, ] + classMatrix[catOrder, ]
    }
  }

  avgClassMatrix = avgClassMatrix / length(list_classMatrix)

  #return
  avgClassMatrix
}

