make_pca_data <- function() {
  set.seed(1)
  counts <- matrix(
    stats::rnbinom(1000, mu = 100, size = 1),
    nrow = 100, ncol = 10,
    dimnames = list(paste0("gene", 1:100), paste0("s", 1:10))
  )
  metadata <- data.frame(
    condition = factor(rep(c("control", "treated"), each = 5)),
    row.names = paste0("s", 1:10)
  )
  list(counts = counts, metadata = metadata)
}

test_that("PCAplot errors when counts/metadata are misaligned", {
  d <- make_pca_data()
  bad_metadata <- d$metadata
  rownames(bad_metadata) <- paste0("x", 1:10)

  expect_error(
    PCAplot(d$counts, bad_metadata, method = "DESeq2"),
    "colnames\\(counts\\) must match"
  )
})

test_that("PCAplot runs end-to-end with DESeq2", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  d <- make_pca_data()

  p <- PCAplot(d$counts, d$metadata, method = "DESeq2", color_by = "condition")

  expect_s3_class(p, "ggplot")
})

test_that("PCAplot runs end-to-end with edgeR", {
  skip_if_not_installed("edgeR")
  d <- make_pca_data()

  p <- PCAplot(d$counts, d$metadata, method = "edgeR", color_by = "condition")

  expect_s3_class(p, "ggplot")
})
