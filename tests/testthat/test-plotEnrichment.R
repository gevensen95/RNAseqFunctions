make_ora_results <- function() {
  data.frame(
    ID = paste0("GO:000", 1:10, ": term ", 1:10),
    geneSet = paste0("GO:000", 1:10),
    description = paste0("term ", 1:10),
    enrichmentRatio = c(3, -3, 2.5, -2, 1.8, -1.5, 1.2, -1.1, 4, -4),
    FDR = c(0.001, 0.002, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.0005, 0.0007),
    stringsAsFactors = FALSE
  )
}

make_gsea_results <- function() {
  data.frame(
    ID = paste0("GO:000", 1:10, ": term ", 1:10),
    geneSet = paste0("GO:000", 1:10),
    description = paste0("term ", 1:10),
    NES = c(2, -2, 1.5, -1.5, 1.2, -1.2, 1.1, -1.1, 3, -3),
    FDR = c(0.001, 0.002, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.0005, 0.0007),
    stringsAsFactors = FALSE
  )
}

test_that("plotEnrichment auto-detects ORA (enrichmentRatio) results", {
  expect_s3_class(plotEnrichment(make_ora_results()), "ggplot")
})

test_that("plotEnrichment auto-detects GSEA (NES) results", {
  expect_s3_class(plotEnrichment(make_gsea_results()), "ggplot")
})

test_that("plotEnrichment errors when no value column can be found", {
  bad_results <- make_ora_results()
  bad_results$enrichmentRatio <- NULL

  expect_error(
    plotEnrichment(bad_results),
    "Could not auto-detect a value column"
  )
})

test_that("plotEnrichment errors when alpha filters out everything", {
  expect_error(
    plotEnrichment(make_ora_results(), alpha = 1e-10),
    "No enrichment terms remain"
  )
})

test_that("plotEnrichment respects an explicit valueCol", {
  results <- make_ora_results()
  names(results)[names(results) == "enrichmentRatio"] <- "customRatio"

  expect_s3_class(plotEnrichment(results, valueCol = "customRatio"), "ggplot")
})
