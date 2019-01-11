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
splitCommon<-function(sampTab, ncells=NA, dLevel){
  cts<-unique(as.vector(sampTab[,dLevel])) #receive the names of the categories
  trainingids<-vector()

  # automatically pick a ncells if that value is missing
  if(is.na(ncells) == TRUE){
    smallestNCells = 1e9 #if someone has that many samples, good for them.

    for(ct in cts) {
      if(nrow(sampTab[sampTab[,dLevel]==ct,]) < smallestNCells) {
        smallestNCells <- nrow(sampTab[sampTab[,dLevel]==ct,])
      }
    }

    ncells <- as.integer(smallestNCells / 2) # select the ncell based on the smallest samples divided by 2
  }

  for(ct in cts){
    #cat(ct,": ") # this one can also be deleted
    stX<-sampTab[sampTab[,dLevel]==ct,]

    # Error catching mechanism
    if (nrow(stX) < ncells) {
      stop(paste0("Category ", ct, " has ", nrow(stX), " samples. Please pick a samller number for ncells. "))
    }
    if (ncells <= 0) {
      stop(paste0("Category ", ct, " has ", nrow(stX), " samples. Please pick a number above 0 for ncells. "))
    }
    # ccount<-nrow(stX)-3 # I don't think I need this line
    # ccount<-min(ccount, ncells) # I don't think I need this line as well

    cat(ct, "has ", nrow(stX), " samples.", "\n")

    trainingids<-append(trainingids, sample(rownames(stX), ncells)) # randomly samples ccount of training samples
  }

  val_ids<-setdiff(rownames(sampTab), trainingids) # the samples that are not used to training are used for validation
  list(train=sampTab[trainingids,], val=sampTab[val_ids,])
}
