test_that("Package loads correctly", {
  expect_true(require(CellHet))
})

test_that("Basic functions exist", {
  expect_true(exists("compareDEGs"))
  expect_true(exists("visualizeHeterogeneity"))
  expect_true(exists("runPathwayAnalysis"))
})