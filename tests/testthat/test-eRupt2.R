make_results <- function(n = 20, zero_padj = FALSE) {
  results <- data.frame(
    Gene = paste0("Gene", seq_len(n)),
    log2FoldChange = seq(-3, 3, length.out = n),
    padj = rep(c(0.001, 0.2), length.out = n)
  )
  if (zero_padj) {
    results$padj[1] <- 0
  }
  results
}

test_that("eRupt2 returns a ggplot object", {
  expect_s3_class(eRupt2(make_results()), "ggplot")
})

test_that("eRupt2 handles padj values of exactly 0 without error", {
  expect_s3_class(eRupt2(make_results(zero_padj = TRUE)), "ggplot")
})

test_that("eRupt2 does not error when fewer than 10 rows are provided", {
  small_results <- make_results(n = 4)
  expect_s3_class(eRupt2(small_results), "ggplot")
})

test_that("eRupt2 respects a custom geneID column", {
  results <- make_results()
  names(results)[names(results) == "Gene"] <- "symbol"
  expect_s3_class(eRupt2(results, geneID = "symbol"), "ggplot")
})
