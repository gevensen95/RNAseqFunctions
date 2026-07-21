test_that("summarizeDEGs tallies a single results data.frame correctly", {
  results <- data.frame(
    Gene = paste0("gene", 1:6),
    log2FoldChange = c(2, -2, 0.1, -0.1, 3, -3),
    padj = c(0.01, 0.01, 0.01, 0.01, 0.2, 0.2)
  )

  summary <- summarizeDEGs(results, alpha = 0.05, FCsig = 1)

  expect_equal(nrow(summary), 1)
  expect_equal(summary$up, 1)
  expect_equal(summary$down, 1)
  expect_equal(summary$not_significant, 4)
  expect_equal(summary$total, 6)
})

test_that("summarizeDEGs handles a named list of comparisons", {
  results_up <- data.frame(
    Gene = paste0("gene", 1:4),
    log2FoldChange = c(2, 2, 2, 2),
    padj = c(0.01, 0.01, 0.01, 0.01)
  )
  results_down <- data.frame(
    Gene = paste0("gene", 1:4),
    log2FoldChange = c(-2, -2, -2, -2),
    padj = c(0.01, 0.01, 0.01, 0.01)
  )

  summary <- summarizeDEGs(list(comp_a = results_up, comp_b = results_down),
                            alpha = 0.05, FCsig = 1)

  expect_equal(nrow(summary), 2)
  expect_equal(summary$comparison, c("comp_a", "comp_b"))
  expect_equal(summary$up, c(4, 0))
  expect_equal(summary$down, c(0, 4))
})

test_that("summarizeDEGs auto-names an unnamed list", {
  results <- data.frame(
    Gene = paste0("gene", 1:4),
    log2FoldChange = c(2, 2, -2, -2),
    padj = c(0.01, 0.01, 0.01, 0.01)
  )

  summary <- summarizeDEGs(list(results, results))

  expect_equal(summary$comparison, c("comparison_1", "comparison_2"))
})
