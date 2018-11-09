#' @title
#' @title
#' Split Sample Table Randomly
#'
#' @description
#' Split the sample table into training set and validation set in random fashion.
#'
#' @param sampTab sample table
#' @param ncells number of samples for training in each category
#' @param dLevel the column name with the classification categories
#' @return a list containing training sample table and validation sample table
#'
#' @export
splitCommon<-function(sampTab, ncells, dLevel="description1"){

  cts<-unique(as.vector(sampTab[,dLevel])) #receive the names of the categories
  trainingids<-vector()

  for(ct in cts){
    cat(ct,": ")
    stX<-sampTab[sampTab[,dLevel]==ct,]
    ccount<-nrow(stX)-3
    ccount<-min(ccount, ncells)
    cat(nrow(stX),"\n")
    trainingids<-append(trainingids, sample(rownames(stX), ccount)) # randomly samples ccount of training samples
  }

  val_ids<-setdiff(rownames(sampTab), trainingids) # the samples that are not used to training are used for validation
  list(train=sampTab[trainingids,], val=sampTab[val_ids,])
}
