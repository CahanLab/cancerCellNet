#' @title
#' Gene Symbol Conversion
#'
#' @description
#' The function can convert Ensembl Gene ID, Ensembl Transcript ID or mouse gene symbol into human gene symbol
#'
#' @param expTab the matrix of experssion data with sample names as column names and gene names as row names
#' @param typeENST TRUE if the gene identifier is in Ensembl Transcript ID
#' @param typeENSG TRUE if the gene identifier is in Ensembl Gene ID
#' @param typeMusGene TRUE if the gene identifier is in mouse gene symbols
#'
#' @return a converted expression matrix
#'
#' @export
utils_convertToGeneSymbols <- function(expTab, typeENST = FALSE, typeENSG= FALSE, typeMusGene = FALSE) {

  #catch errors
  if (typeENSG + typeENST + typeMusGene > 1) {
    stop("Please only indicate one type of conversion");
  }
  #catch errors
  if (typeENSG + typeENST + typeMusGene == 0) {
    stop("Please assign one type of conversion");
  }


  if (typeENST == TRUE) {
    convertMatrix = conversionList[[1]]

    nOGgenes = nrow(expTab) # the number of original genes
    rownames(convertMatrix) = as.vector(convertMatrix$ensembl_transcript_id)
    cgenes = intersect(as.vector(convertMatrix$ensembl_transcript_id), rownames(expTab))
    expTab = as.matrix(expTab[cgenes,])
    rownames(expTab) = as.vector(convertMatrix[cgenes,]$gene_symbol)

    nGenesFailed = nOGgenes - nrow(expTab)
    print(paste0("Could not convert ", nGenesFailed, " genes."))

    returnMatrix = expTab

  } else if (typeENSG == TRUE) {
    convertMatrix = conversionList[[2]]

    nOGgenes = nrow(expTab) # the number of original genes
    rownames(convertMatrix) = as.vector(convertMatrix$ensembl_gene_id)
    cgenes = intersect(as.vector(convertMatrix$ensembl_gene_id), rownames(expTab))
    expTab = as.matrix(expTab[cgenes,])
    rownames(expTab) = as.vector(convertMatrix[cgenes,]$gene_symbol)

    nGenesFailed = nOGgenes - nrow(expTab)
    print(paste0("Could not convert ", nGenesFailed, " genes."))

    returnMatrix = expTab

  } else if (typeMusGene == TRUE) {
    convertMatrix = conversionList[[3]]

    nOGgenes = nrow(expTab) # the number of original genes
    rownames(convertMatrix) = as.vector(convertMatrix$mouse_gene)
    cgenes = intersect(as.vector(convertMatrix$mouse_gene), rownames(expTab))

    expTab = as.matrix(expTab[cgenes,])
    rownames(expTab) = as.vector(convertMatrix[cgenes,]$gene_symbol)

    nGenesFailed = nOGgenes - nrow(expTab)
    print(paste0("Could not convert ", nGenesFailed, " genes."))

    returnMatrix = expTab
  }
  #return
  returnMatrix
}
