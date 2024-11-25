library(RhodopXin)

test_that("loadSequence() returns the correct output", {
  filepath <- system.file("extdata", "rhodopsins.fasta", package = "RhodopXin")
  results <- loadSequence(filePath = filepath)
  expect_s4_class(results, "AAStringSet")
  expect_length(results, 3)
})

test_that("loadSequence() throws an error when no sequences have been loaded
          from the FASTA file", {
  # Create an empty FASTA file
  empty_file <- "empty.fasta"
  file.create(empty_file)

  expect_error(loadSequence(filePath = empty_file), "Empty file: no sequences detected")

  # Remove empty test file
  file.remove(empty_file)
})

test_that("loadSequence() throws an error when the wrong file format is passed
          as an argument", {
  # Create a csv file
  csv_file <- "test.csv"
  file.create(csv_file)

  expect_error(loadSequence(filePath = csv_file), "Incorrect file format: must be one of '.fa', '.fasta', '.txt'")

  # Remove empty test file
  file.remove(csv_file)
})

test_that("loadSequence() throws an error when the file path provided is not a string", {
  expect_error(loadSequence(filePath = 12345), "Invalid parameter type: filePath must be of type character")
  expect_error(loadSequence(filePath = TRUE), "Invalid parameter type: filePath must be of type character")
  expect_error(loadSequence(filePath = data.frame()), "Invalid parameter type: filePath must be of type character")
})

test_that("loadSequence() throws an error when more than 1 input is provided", {
  expect_error(loadSequence(filePath = c("file", "path")), "More than 1 input: provide only one string for the file path")
})
