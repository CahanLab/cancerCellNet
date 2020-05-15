#' @title
#' Split training set and validation set
#'
#' @description
#' Split the sample table into training and validation by proportion
#'
#' @param sampTab the sample table
#' @param proportion the fraction of all available data in a category that you want to use for training
#' @param dLevel the name of the sample table column with the classification categories
#' @return a list containing training sample table and validation sample table
#'
#' @export
splitCommon_proportion <- function(sampTab, proportion, dLevel) {
  trainSampleTable = NULL
  for(myCat in unique(sampTab[, dLevel])) {

    tempSt = sampTab[sampTab[, dLevel] == myCat, ]
    tempSt = tempSt[sample(rownames(tempSt), as.integer((proportion)*nrow(tempSt))), ]

    if(is.null(trainSampleTable) == TRUE) {
      trainSampleTable = tempSt
    }
    else {
      trainSampleTable = rbind(trainSampleTable, tempSt)
    }
  }

  valSampleTable = sampTab[!(rownames(sampTab) %in% rownames(trainSampleTable)), ]

  return(list("trainingSet" = trainSampleTable, "validationSet" = valSampleTable))
}
