make_ma_results <- function() {
  set.seed(1)
  data.frame(
    baseMean = stats::runif(50, 10, 1000),
    log2FoldChange = c(seq(-3, -0.2, length.out = 25), seq(0.2, 3, length.out = 25)),
    padj = rep(c(0.001, 0.2), length.out = 50),
    Gene = paste0("gene", 1:50)
  )
}

test_that("MAplot returns a ggplot object", {
  expect_s3_class(MAplot(make_ma_results()), "ggplot")
})

test_that("MAplot handles all-NS results without error", {
  results <- make_ma_results()
  results$padj <- 1
  expect_s3_class(MAplot(results), "ggplot")
})
