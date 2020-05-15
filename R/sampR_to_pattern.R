#' @title
#' Create Template
#'
#' @description
#' Create a template of the cancer sample in terms of cancer category.
#'
#' @param sampR a vector of the sample catgories
#'
#' @return a list of vectors indicating mapping of each cancer type in the sampR vector
#' @export
sc_sampR_to_pattern<-function(sampR){
  d_ids = unique(as.vector(sampR))
  nnnc = length(sampR)

  ans = list()
  for(d_id in d_ids){
    x = rep(0,nnnc)
    x[which(sampR==d_id)] = 1
    ans[[d_id]] = x
  }

  return(ans)
}
