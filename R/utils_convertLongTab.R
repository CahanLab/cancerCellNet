#' @title
#' Convert Long Expression Table
#'
#' @description
#' Convert long gene expression/count table into wide dataframe to proceed to the next step
#' of training classifier or validating samples.
#'
#' @param longDf a dataframe of gene expression/count in long table format
#' @param SampCol the name of the column containing sample names
#' @param geneCol the name of the column containing gene names
#' @param geneExpcol the name of the column containing gene counts or expression
#' @param geneNaOmit TRUE if you want to omit the genes with NA
#'
#' @return a list containing the gene expression dataframe in wide table format and a dataframe of omitted genes if geneNaOmit is TRUE
#' Otherwise, just returns gene Expression matrix in wide table
#' @export
#'
#' @importFrom stringr str_replace
#'
utils_convertLongTab <- function(longDf, SampCol, geneCol, geneExpCol, geneNaOmit = TRUE) {
  longDf_extract = data.frame(longDf[, c(SampCol, geneCol, geneExpCol) ])
  exp_tab_rs = reshape(longDf_extract, idvar = geneCol, timevar = SampCol, direction = "wide")
  colnames(exp_tab_rs) = stringr::str_replace(colnames(exp_tab_rs), paste0(geneExpCol, "."), "")
  rownames(exp_tab_rs) = exp_tab_rs[, geneCol]
  exp_tab_rs[, geneCol] = NULL

  exp_tab_omitted <- data.frame()

  if (geneNaOmit == TRUE) {
    exp_tab_clean = na.omit(exp_tab_rs);

    exp_tab_omitted = exp_tab_rs[!(rownames(exp_tab_rs) %in% rownames(exp_tab_clean)), ]

    output = list(wideTable = exp_tab_clean, omittedGenesTable = exp_tab_omitted);

  }
  else {
    output = exp_tab_rs
  }

  output
}
