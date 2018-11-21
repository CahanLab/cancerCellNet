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
utils_matchSampTab <- function(expDat, sampTab, sampID, sampDesc) {

  # if all the samples in table corresponds to
  if(all(colnames(expDat) %in% sampTab[, sampID]) == TRUE) {
    print("All the samples in expression matrix is represented in sample table.")
    sortedSampTab = sampTab[order(sampTab[, sampDesc]), ]
    sortedSampID = sortedSampTab[, sampID]
    sortedExpDat = expDat[, sortedSampleID]
    returnList = list(matchedExpDat = sortedExpDat, matchedSampTab = sortedSampTab)
  }

  else {
    matchSamples = intersect(colnames(expDat), sampTab[, sampID])
    expDat = expDat[, matchSamples] #only take out the ones that are matching
    sampTab = sampTab[(sampTab[, sampID] %in% matchSamples), ]

    sortedSampTab = sampTab[order(sampTab[, sampDesc]), ]
    sortedSampID = sortedSampTab[, sampID]
    sortedExpDat = expDat[, sortedSampleID]

    # Select out the expDat samples that are not matched
    unmatchedExpDatID = setdiff(matchSamples, colnames(expDat))
    unmatchedExpDat = epxDat[, unmatchedSample]

    if (length(unmatchedExpDatID) > 0) {
      print("There are unmatched samples in expression matrix. ")
    }

    # select the samples in sampleTab that are not matched
    unmatchedSampleID = setdiff(matchSamples, sampTab[, sampID])
    unmatchedSampleTab = sampTab[sampTab[, sampID] %in% unmatchedSampleID, ]

    if (length(unmatchedSampleTab) > 0) {
      print("There are unmatched samples in sample table. ")
    }

    returnList = list(matchedExpDat = sortedExpDat, matchedSampTab = sortedSampTab, unmatchedExpDat = unmatchedExpDat, unmatchedSampleTab = unmatchedSampleTab)
  }

  #return
  returnList
}

i = utils_matchSampTab(sampTab = test_df, sampDesc = "sample2")


newIndex = c("sample2", "sample1", "sample4", "sample3", "sample6", "sample5", "sample8", "sample7", "sample9", "sample10")

newTest_df = test_df[, newIndex]



set.seed(626) # set random seed to generate
test_df = data.frame(replicate(10, sample(1:4, 50, rep=TRUE)))
col_name = c()
for (i in 1:ncol(test_df)){
  col_name = c(col_name, paste0("sample", i))
}
colnames(test_df) = col_name
row_name = c()
for (i in 1:nrow(test_df)) {
  row_name = c(row_name, paste0("gene", i))
}
rownames(test_df) = row_name
