make_heatmap_data <- function() {
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
  results <- data.frame(
    Gene = paste0("gene", 1:100),
    log2FoldChange = stats::rnorm(100),
    padj = stats::runif(100)
  )
  list(counts = counts, metadata = metadata, results = results)
}

test_that("DEGheatmap errors when counts/metadata are misaligned", {
  d <- make_heatmap_data()
  bad_metadata <- d$metadata
  rownames(bad_metadata) <- paste0("x", 1:10)

  expect_error(
    DEGheatmap(d$counts, bad_metadata, d$results, method = "DESeq2"),
    "colnames\\(counts\\) must match"
  )
})

test_that("DEGheatmap errors when no genes in results match counts", {
  d <- make_heatmap_data()
  d$results$Gene <- paste0("notagene", 1:100)

  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  expect_error(
    DEGheatmap(d$counts, d$metadata, d$results, method = "DESeq2"),
    "None of the top genes"
  )
})

test_that("DEGheatmap runs end-to-end with DESeq2", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  d <- make_heatmap_data()

  p <- DEGheatmap(d$counts, d$metadata, d$results, method = "DESeq2", n = 10)

  expect_s3_class(p, "ggplot")
})

test_that("DEGheatmap runs end-to-end with annotate_by", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  d <- make_heatmap_data()

  p <- DEGheatmap(d$counts, d$metadata, d$results, method = "DESeq2", n = 10,
                   annotate_by = "condition")

  expect_s3_class(p, "ggplot")
})
