make_counts_metadata <- function() {
  set.seed(1)
  counts <- matrix(
    stats::rnbinom(500, mu = 100, size = 1),
    nrow = 50, ncol = 10,
    dimnames = list(paste0("gene", 1:50), paste0("s", 1:10))
  )
  metadata <- data.frame(
    condition = factor(rep(c("control", "treated"), each = 5),
                        levels = c("control", "treated")),
    row.names = paste0("s", 1:10)
  )
  list(counts = counts, metadata = metadata)
}

test_that("makeComp errors when counts/metadata are misaligned", {
  d <- make_counts_metadata()
  bad_metadata <- d$metadata
  rownames(bad_metadata) <- paste0("x", 1:10)

  expect_error(
    makeComp(d$counts, bad_metadata, method = "DESeq2"),
    "colnames\\(counts\\) must match"
  )
})

test_that("makeComp errors when formula is missing and no condition column exists", {
  d <- make_counts_metadata()
  names(d$metadata) <- "group"

  expect_error(
    makeComp(d$counts, d$metadata, method = "DESeq2"),
    "no 'condition' column"
  )
})

test_that("makeComp gives an informative error when DESeq2 isn't installed", {
  skip_if(requireNamespace("DESeq2", quietly = TRUE))
  d <- make_counts_metadata()

  expect_error(
    makeComp(d$counts, d$metadata, method = "DESeq2"),
    "BiocManager::install\\(\"DESeq2\"\\)"
  )
})

test_that("makeComp gives an informative error when edgeR isn't installed", {
  skip_if(requireNamespace("edgeR", quietly = TRUE))
  d <- make_counts_metadata()

  expect_error(
    makeComp(d$counts, d$metadata, method = "edgeR"),
    "BiocManager::install\\(\"edgeR\"\\)"
  )
})

test_that("makeComp runs end-to-end with DESeq2 and returns the expected shape", {
  skip_if_not_installed("DESeq2")
  d <- make_counts_metadata()

  res <- makeComp(d$counts, d$metadata, method = "DESeq2")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(d$counts))
  expect_equal(colnames(res)[2], "log2FoldChange")
  expect_equal(colnames(res)[7], "Gene")
  expect_true(all(c("Gene", "padj", "log2FoldChange") %in% colnames(res)))
})

test_that("makeComp runs end-to-end with edgeR and returns the expected shape", {
  skip_if_not_installed("edgeR")
  d <- make_counts_metadata()

  res <- makeComp(d$counts, d$metadata, method = "edgeR")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(d$counts))
  expect_equal(colnames(res)[2], "log2FoldChange")
  expect_equal(colnames(res)[7], "Gene")
  expect_true(all(c("Gene", "padj", "log2FoldChange") %in% colnames(res)))
})

test_that("makeComp accepts an explicit DESeq2 contrast", {
  skip_if_not_installed("DESeq2")
  d <- make_counts_metadata()

  res <- makeComp(d$counts, d$metadata, method = "DESeq2",
                   contrast = c("condition", "treated", "control"))

  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res)[2], "log2FoldChange")
})
