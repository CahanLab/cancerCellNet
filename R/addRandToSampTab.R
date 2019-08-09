#' @title
#' Add Random Profiles to Sample Table
#'
#' @description
#' Add the random profiles generated from \code{\link{rf_classPredict}} to the sample table
#' for classification of the random profiles.
#' @param classRes the classification matrix generated from \code{\link{rf_classPredict}}
#' @param sampTab the meta sample table
#' @param desc the column name of the column containing the sample categories in the meta sample table
#' @param id the column name of the column containing the sample names in the meta sample table
#' @return meta sample table containing the random profiles
#' @examples
#' stValRand<-addRandToSampTab(classRes_val, stVal, "description2", "sample_name")
#' @export
addRandToSampTab<-function(classRes, sampTab, desc, id) {
  cNames<-colnames(classRes)
  snames<-rownames(sampTab)

  rnames<-setdiff(cNames, snames)

  cat("number of random samples: ",length(rnames), "\n")

  stNew<-data.frame(rid=rnames, rdesc=rep("rand", length(rnames)))
  stTop<-sampTab[,c(id, desc)]
  colnames(stNew)<-c(id, desc)

  ans<-rbind(stTop, stNew)
  rownames(ans)<-colnames(classRes)

  #return
  ans
}
