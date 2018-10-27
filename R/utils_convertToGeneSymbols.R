#' @title
#' Gene Symbol Conversion
#'
#' @description
#' The function can convert Ensembl Gene ID, Ensembl Transcript ID or mouse gene symbol into human gene symbol
#'
#' @param expdf the dataframe of experssion data with sample names as column names and gene names as row names
#' @param typeENSG TRUE if the gene identifier is in Ensembl Transcript ID
#' @param typeENSG TRUE if the gene identifier is in Ensembl Gene ID
#' @param typeMusGene TRUE if the gene identifier is in mouse gene symbols
#'
#' @return a list containing the converted experssion table as the first element and the expression table with genes that cannot be converted with our conversion dataset
utils_convertToGeneSymbols <- function(expdf, typeENST = FALSE, typeENSG= FALSE, typeMusGene = FALSE) {

  #catch errors
  if (typeENSG + typeENST + typeMusGene > 1) {
    message("Please only indicate one type of conversion");
  }
  #catch errors
  if (typeENSG + typeENST + typeMusGene == 0) {
    message("Please indicate one type of conversion");
  }


  if (typeENST == TRUE) {
    convertMatrix = conversionList[[1]]
    expdf_extract = expdf[rownames(expdf) %in% convertMatrix$ensembl_transcript_id, ]
    expdf_leftover = expdf[!(rownames(expdf_extract) %in% rownames(expdf)),  ]

    tempName = data.frame(ensembl_transcript_id = rownames(expdf_extract))
    tempName = merge(tempName, convertMatrix, "ensembl_transcript_id")

    rownames(expdf_extract) = tempName$gene_symbol

    returnList = list(convertedDF = expdf_extract, unableToConvert = expdf_leftover)

  } else if (typeENSG == TRUE) {
    convertMatrix = conversionList[[2]]

    expdf_extract = expdf[rownames(expdf) %in% convertMatrix$ensembl_gene_id, ]
    expdf_leftover = expdf[!(rownames(expdf_extract) %in% rownames(expdf)),  ]

    tempName = data.frame(ensembl_gene_id = rownames(expdf_extract))
    tempName = merge(tempName, convertMatrix, "ensembl_gene_id")

    rownames(expdf_extract) = tempName$gene_symbol

    returnList = list(convertedDF = expdf_extract, unableToConvert = expdf_leftover)
  } else if (typeMusGene == TRUE) {
    convertMatrix = conversionList[[3]]

    expdf_extract = expdf[rownames(expdf) %in% convertMatrix$mouse_gene, ]
    expdf_leftover = expdf[!(rownames(expdf_extract) %in% rownames(expdf)),  ]

    tempName = data.frame(mouse_gene = rownames(expdf_extract))
    tempName = merge(tempName, convertMatrix, "mouse_gene")

    rownames(expdf_extract) = tempName$gene_symbol

    returnList = list(convertedDF = expdf_extract, unableToConvert = expdf_leftover)
  }
  #return
  returnList
}
