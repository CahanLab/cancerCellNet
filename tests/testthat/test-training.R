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

test_that("Let's test the implementation of compMu", {
  set.seed(626) # set random seed to generate
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


  # test the sc_compmu aspect of statTab
  # i didn't test from 1:3 because my test is desigend such that NaN occurs. I don't want to completely
  # rewrite the function
  for (i in 1:2) {
    test_statTab = sc_statTab(test_df, i)

    myFunc <- function(vector, threshold) {
      mean(vector[which(vector > threshold)])
    }

    myFunc_apply = apply(test_df, 1, myFunc, threshold = i)
    for(a in 1:length(test_statTab$mu)) {

      expect_equal(test_statTab$mu[a], myFunc_apply[[a]])
    }
  }
})

test_that("Test SC_filterGenes", {
  set.seed(626) # set random seed to generate
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

  # patrick's implementation
  patrick_filterGenes<-function(geneStats, alpha1=0.1, alpha2=0.01, mu=2){
    passing1<-rownames(geneStats[geneStats$alpha>alpha1,])
    notPassing<-setdiff(rownames(geneStats), passing1)
    geneStats<-geneStats[notPassing,]
    c(passing1, rownames(geneStats[which(geneStats$alpha>alpha2 & geneStats$mu>mu),]))

    myTest = testStatTab[(testStatTab$alpha > 0.6) | (testStatTab$alpha > 0 & testStatTab$mu > 3), ]

  }

  testStatTab = sc_statTab(test_df, 1)
  test_filter = sc_filterGenes(testStatTab, alpha1 = 0.6, alpha2 = 0, mu = 3)

  myTest = patrick_filterGenes(test_df, alpha1 = 0.6, alpha2 = 0, mu = 3)

  for(i in 1:length(test_filter)) {
    test_filter = sort(test_filter)
    expect_equal(sort(rownames(myTest))[i], sort(test_filter)[i])
  }
})

# cannot think of meaningful tests
test_that("Test sampR_to_pattner", {

  # never mind, it's kind of trivial. Until I find more meaningful test, I will leave it blank
  cat1 = rep("cancer1", 4)
  cat2 = rep("cancer2", 7)
  cat3 = rep("caner3", 2)



})

# cannot think of meaningful tests as well
test_that("Test testPattern", {
  # almost trivial as well. I will come back when I think of a better test


})

test_that("Error case for findClassGenes", {

  set.seed(626) # set random seed to generate
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

  cat1 = rep("cancer1", 5)
  cat2 = rep("cancer2", 7)
  cat3 = rep("cancer3", 3)

  bigCat = c(cat1, cat2, cat3)
  test_sampTab = data.frame(categories = bigCat)

  expect_error(findClassyGenes(expDat = test_df, sampTab = test_sampTab))

})

test_that("ptGetTop comparison between the old one and the new one", {
  # one suggestion
  set.seed(626) # set random seed to generate
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

  grp = list("c1", "c1", "c2", "c2", "c3", "c3", "c4", "c4", "c5", "c5")
  names(grp) = colnames(test_df)


  result = ptGetTop_old(expDat = test_df, topX = 20, cell_labels = grp)
  result2 = ptGetTop(expDat = test_df, topX = 20, cell_labels = grp)

  expect_equal(sort(result), sort(result2))

  for (i in 1:length(result)) {
    expect_equal(sort(result)[i], sort(result2)[i])
  }

})

test_that("ptGetTop Error case. xTop TOO Large", {
  # one suggestion
  set.seed(626) # set random seed to generate
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

  grp = list("c1", "c1", "c2", "c2", "c3", "c3", "c4", "c4", "c5", "c5")
  names(grp) = colnames(test_df)

  expect_error(ptGetTop(expDat = test_df, topX = 9e100, cell_labels = grp))

})


