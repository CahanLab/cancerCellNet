#' Load an R object
#'
#' Load an R object.
#' @param fname The file name of the R object
#'
#' @return A R object
#'
#' @export
utils_loadObject<-function(fname) {
  x = load(fname)
  return(get(x))
}
