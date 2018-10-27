#' Load an R object
#'
#' Load an R object
#' @param fname The file name of the R object
#'
#' @return A R object
#' @example
#' utils_loadObject("ccn_classifier_Jun_29_2018.rda")
#'
#' @export
utils_loadObject<-function (fname){
  x<-load(fname);
  get(x);
}
