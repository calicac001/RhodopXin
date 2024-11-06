library(RhodopXin)

test_that("loadSequence() returns the correct output", {
  filepath <- system.file("extdata", "rhodopsins.fasta", package = "RhodopXin")
  results <- loadSequence(inputStr = filepath)
  expect_s4_class(results, "AAStringSet")
  expect_length(results, 3)
})
