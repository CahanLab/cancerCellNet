context("Just Testing The functions necessary during validation stage")
library(testthat)

test_that("sampleTest", {
  expect_equal("hello", "hello")
})

test_that("addRandToSampleTab", {
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

  addRandomMatrix = data.frame(replicate(10, sample(1:4, 50, rep=TRUE)))
  col_name = c()
  for (i in 1:ncol(addRandomMatrix)){
    col_name = c(col_name, paste0("Rand", i))
  }
  colnames(addRandomMatrix) = col_name
  row_name = c()
  for (i in 1:nrow(addRandomMatrix)) {
    row_name = c(row_name, paste0("gene", i))
  }
  rownames(addRandomMatrix) = row_name

  bigMatrix = cbind(test_df, addRandomMatrix)


  grp_mat = data.frame(sample = names(grp), cat = unlist(grp))

  stValRand<-addRandToSampTab(bigMatrix, grp_mat, "cat", "sample")


  ranSample = colnames(addRandomMatrix)
  ManualAddRand = data.frame(sample = ranSample, cat = "rand")
  rownames(ManualAddRand) = colnames(addRandomMatrix)
  manualTest = rbind(grp_mat, ManualAddRand)

  expect_true(all.equal(manualTest, stValRand))
})
