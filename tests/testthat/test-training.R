context("Just Testing The functions necessary during trainng stage")
library(testthat)

test_that("whether the utils_convertToGeneSymbols", {
  conversionList = conversionList
  EnsToGene= conversionList[[1]][1:3, ]
  rownames(EnsToGene) = EnsToGene$ensembl_transcript_id

  #normal cases
  ans = utils_convertToGeneSymbols(EnsToGene, typeENST = TRUE)[[1]]
  expect_equal("OR4F5", as.character(rownames(ans)[1]))
  expect_equal("FO538757.3", as.character(rownames(ans)[2]))
  expect_equal("FO538757.2", as.character(rownames(ans)[3]))

  EnsgToGene = conversionList[[2]][1:3, ]
  rownames(EnsgToGene) = EnsgToGene$ensembl_gene_id

  #normal cases
  ans = utils_convertToGeneSymbols(EnsgToGene, typeENSG = TRUE)[[1]]
  expect_equal("RF00100", as.character(rownames(ans)[1]))
  expect_equal("RNU4-59P", as.character(rownames(ans)[2]))
  expect_equal("SNORD114-2", as.character(rownames(ans)[3]))

  musToGene = conversionList[[3]][1:3, ]
  rownames(musToGene) = musToGene$mouse_gene

  #normal cases
  ans = utils_convertToGeneSymbols(musToGene, typeMusGene = TRUE)[[1]]
  expect_equal("SERPINB10", as.character(rownames(ans)[1]))
  expect_equal("HYI", as.character(rownames(ans)[2]))
  expect_equal("GBP7", as.character(rownames(ans)[3]))
})

test_that("Let's test the splitCommons functions", {
  # set up for sample table for the lol
  cat1 = data.frame(sample = 1:10, category = "cat_1")
  cat2 = data.frame(sample = 1:20, category = "cat_2")
  cat3 = data.frame(sample = 1:10, category = "cat_3")

  big_cat = rbind(cat1, cat2, cat3)

  #normal case
  withoutNcells = splitCommon(big_cat, dLevel = "category")
  expect_equal(15, nrow(withoutNcells[[1]]))
  expect_equal(40 - 15, nrow(withoutNcells[[2]]))


  #normal case
  withoutNcells = splitCommon(big_cat, ncells = 3, dLevel = "category")
  expect_equal(9, nrow(withoutNcells[[1]]))
  expect_equal(40 - 9, nrow(withoutNcells[[2]]))


  #edge case
  expect_error(splitCommon(big_cat, ncells = 10, dLevel = "category"))

  #fuck users who do this
  expect_error(splitCommon(big_cat, ncells = 0, dLevel = "category"))

})


test_that("Test my implementation of computing Alpha with Patrick's ", {

  # Patrick's implementation had a tiny bug. When we set the threshold to 0
  # in this dataset, Patrick's implement will return a matrix with the WRONG alpha
  set.seed(88) # set random seed to generate
  test_df = data.frame(replicate(10, sample(1:4, 50, rep=TRUE)))

  col_name = c()
  for (i in 1:ncol(test_df)){
    col_name = c(col_name, paste0("sample", i))
  }
  colnames(test_df) = col_name

  row_name = c()
  for (i in 1:nrow(test_df)) {
    row_name = c(row_name, paste0("gene", i))
  }

  rownames(test_df) = row_name

  # now that the sample test_df is set up. Start Testing

  # Patrick's old code
  Patrick_compAlpha<-function(expMat, threshold=0,pseudo=FALSE){

    # identify the index of vectors greater than threshold
    indexFunction<-function(vector, threshold){
      names(which(vector>threshold));
    }

    indexes<-apply(expMat, 1, indexFunction, threshold);

    alphas<-unlist(lapply(indexes, length));

    ans<-alphas/ncol(expMat)

    if(pseudo){
      ans<-(alphas+1)/(ncol(expMat)+1)
    }

    ans
  }


  for (i in 1:3) {
    mytest = sc_compAlpha(test_df, i)
    PatrickTest = Patrick_compAlpha(test_df, i)

    for(a in 1:length(mytest)) {
      expect_equal(mytest[[a]], PatrickTest[[a]])
    }
  }


})

