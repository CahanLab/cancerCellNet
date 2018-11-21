#' @title
#' Match Sample Table
#'
#' @description
#' To match the ordering of the sample names in expression matrix with the sample table.
#'
#' @param expDat the expression matrix
#' @param sampTab the sample table
#' @param sampID the name of the column in sampTab that contains the sample names
#' @param sampDesc the name of the column in sampTab that contains the sample groups/description
#'
#' @return a list of ordered ordered expression matrix, sample table. If there are unmatched expression matrix or sample tabble
#' the list will also contain unmatched expression matrix and unmatched sample table.
#'
#' @export
utils_matchSampTab <- function(expDat, sampTab, sampID, sampDesc) {

  # if all the samples in table corresponds to
  if(all(colnames(expDat) %in% sampTab[, sampID]) == TRUE) {

    print("All the samples in expression matrix is represented in sample table.")

    sortedSampTab = sampTab[order(sampTab[, sampDesc]), ]

    sortedSampleID = sortedSampTab[, sampID]

    sortedSampleID = as.vector(sortedSampleID)

    sortedExpDat = expDat[, sortedSampleID]
    returnList = list(matchedExpDat = sortedExpDat, matchedSampTab = sortedSampTab)
  }

  else {

    matchSamples = intersect(colnames(expDat), sampTab[, sampID])

    old_expDat = expDat
    expDat = expDat[, matchSamples] #only take out the ones that are matching

    old_sampTab = sampTab
    sampTab = sampTab[(sampTab[, sampID] %in% matchSamples), ]

    sortedSampTab = sampTab[order(sampTab[, sampDesc]), ]

    sortedSampleID = sortedSampTab[, sampID]
    sortedSampleID = as.vector(sortedSampleID)

    sortedExpDat = expDat[, sortedSampleID]

    # Select out the expDat samples that are not matched
    unmatchedExpDatID = setdiff(colnames(old_expDat), matchSamples)

    unmatchedExpDat = NULL

    if (length(unmatchedExpDatID) > 0) {
      print("There are unmatched samples in expression matrix. ")
      print(unmatchedExpDatID)

      unmatchedExpDat = old_expDat[, unmatchedExpDatID]
    }


    # select the samples in sampleTab that are not matched
    unmatchedSampleID = setdiff(old_sampTab[, sampID], matchSamples)
    unmatchedSampleTab = NULL
    if (length(unmatchedSampleID) > 0) {
      print("There are unmatched samples in sample table. ")
      print(unmatchedSampleID)

      unmatchedSampleTab = old_sampTab[old_sampTab[, sampID] %in% unmatchedSampleID, ]
    }

    returnList = list(matchedExpDat = sortedExpDat, matchedSampTab = sortedSampTab, unmatchedExpDat = unmatchedExpDat, unmatchedSampleTab = unmatchedSampleTab)
  }

  #return
  returnList
}
