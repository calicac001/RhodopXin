library(RhodopXin)

# Test createHelixAlignments()
test_that("createHelixAlignments() returns a list of AAStringSet objects", {
  template <- template_rhodopsins[1]
  sequences <- sample_rhodopsins
  rcsb_id <- "1QHJ"

  results <- createHelixAlignments(template = template, sequences = sequences,
                                    rcsb_id = rcsb_id)
  expect_type(results, "list")
  expect_true(all(sapply(results, inherits, "AAStringSet")))
})

test_that("createHelixAlignments() returns the correct number of AAStringSet Objects", {
  template <- template_rhodopsins[1]
  sequences <- sample_rhodopsins
  rcsb_id <- "1QHJ"

  results <- createHelixAlignments(template = template, sequences = sequences,
                                    rcsb_id = rcsb_id)
  expect_length(results, 7)

  template <- template_rhodopsins[2]
  rcsb_id <- "3UG9"

  results <- createHelixAlignments(template = template, sequences = sequences,
                                    rcsb_id = rcsb_id)
  expect_length(results, 9)
})

test_that("createHelixAlignments() returns the AAStringSet objects with the correct
          number of sequences", {
  template <- template_rhodopsins[1]
  sequences <- sample_rhodopsins
  rcsb_id <- "1QHJ"

  results <- createHelixAlignments(template = template, sequences = sequences,
                            rcsb_id = rcsb_id)
  expect_true(all(sapply(results, length) == 4))
})

# Test helixSequences()
test_that("helixSequences() returns a list of AAStringSet objects", {
  template <- template_rhodopsins[1]
  rcsb_id <- "1QHJ"

  results <- helixSequences(template = template, rcsb_id = rcsb_id)
  expect_type(results, "list")
  expect_true(all(sapply(results, inherits, "AAStringSet")))
})

test_that("helixSequences() returns the correct number of AAStringSet Objects", {
  template <- template_rhodopsins[1]
  rcsb_id <- "1QHJ"

  results <- helixSequences(template = template, rcsb_id = rcsb_id)
  expect_length(results, 7)

  template <- template_rhodopsins[2]
  rcsb_id <- "3UG9"

  results <- helixSequences(template = template, rcsb_id = rcsb_id)
  expect_length(results, 9)
})

test_that("helixSequences() returns the AAStringSet objects with the correct
          number of sequences", {
  template <- template_rhodopsins[1]
  rcsb_id <- "1QHJ"

  results <- helixSequences(template = template, rcsb_id = rcsb_id)
  expect_true(all(sapply(results, length) == 1))
})
