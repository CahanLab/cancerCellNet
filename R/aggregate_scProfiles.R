#' @title Aggregate scRNA Profiles
#' @description Aggregate scRNA expression profiles that have the same cluster label or cell type label
#'
#' @param scRNA_exp scRNA expression profiles
#' @param scRNA_st scRNA expression profiles sample table
#' @param cell_id_col column name for cell id
#' @param group_id_col column name for cluster label or cell type label
#'
#' @return a list consists of aggregated expression profiles and aggregated sample table
#' @export
aggregate_scProfiles <- function(scRNA_exp, scRNA_st, cell_id_col, group_id_col) {
  rownames(scRNA_st) = as.vector(scRNA_st[, cell_id_col])
  sumMatrix = matrix(data = NA, nrow = nrow(scRNA_exp), ncol = length(unique(scRNA_st[, group_id_col])))
  colnames(sumMatrix) = unique(scRNA_st[, group_id_col])
  for(cellType in unique(as.vector(scRNA_st[, group_id_col]))) {
    temp_st = scRNA_st[scRNA_st[, group_id_col] == cellType, ]
    temp_exp = scRNA_exp[, rownames(temp_st)]
    sumMatrix[, cellType] = apply(temp_exp, MARGIN = 1, FUN = sum)
  }

  rownames(sumMatrix) = names(apply(temp_exp, MARGIN = 1, FUN = sum))
  sumMatrix_st = data.frame(cell_type = colnames(sumMatrix),
                            cell_id = colnames(sumMatrix))

  return(list(agg_exp = sumMatrix, agg_st = sumMatrix_st))
}
