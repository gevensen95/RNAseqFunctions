# These tests mock WebGestaltR::WebGestaltR() via testthat's
# local_mocked_bindings(.package = ) (requires testthat >= 3.1.4) so they
# run offline and exercise WebGestalt_Pipe2()'s own logic — argument
# construction, NA handling, and result combination — without hitting the
# real WebGestalt web service.

test_that("WebGestalt_Pipe2 warns on an unrecognized mode", {
  results <- data.frame(
    Gene = c("Gene1", "Gene2"),
    log2FoldChange = c(1, -1),
    padj = c(0.01, 0.01)
  )

  expect_warning(
    WebGestalt_Pipe2(mode = "bogus", results = results),
    "mode not recognized"
  )
})

test_that("WebGestalt_Pipe2 ORA returns results when only one direction has hits", {
  # Regression test: previously this combination (down = 0 rows, up > 0
  # rows) fell through to `return(ORA)` with `ORA` never having been
  # assigned, raising "object 'ORA' not found".
  results <- data.frame(
    Gene = paste0("gene", 1:6),
    log2FoldChange = c(2, 2, 2, -2, -2, -2),
    padj = c(0.01, 0.01, 0.01, 0.2, 0.2, 0.2)
  )

  fake_webgestaltr <- function(interestGene, ...) {
    if (length(interestGene) == 0) {
      return(data.frame())
    }
    data.frame(
      geneSet = "GO:0000001",
      description = "fake term",
      enrichmentRatio = 2,
      FDR = 0.01,
      stringsAsFactors = FALSE
    )
  }
  testthat::local_mocked_bindings(WebGestaltR = fake_webgestaltr, .package = "WebGestaltR")

  res <- WebGestalt_Pipe2(mode = "ORA", results = results, alpha = 0.05, FC = 1)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_equal(colnames(res)[1], "ID")
})

test_that("WebGestalt_Pipe2 ORA messages and returns NULL when neither direction has hits", {
  results <- data.frame(
    Gene = paste0("gene", 1:4),
    log2FoldChange = c(0.1, -0.1, 0.2, -0.2),
    padj = c(0.9, 0.9, 0.9, 0.9)
  )

  fake_webgestaltr <- function(interestGene, ...) data.frame()
  testthat::local_mocked_bindings(WebGestaltR = fake_webgestaltr, .package = "WebGestaltR")

  expect_message(
    res <- WebGestalt_Pipe2(mode = "ORA", results = results, alpha = 0.05, FC = 1),
    "No ORA terms are returned"
  )
  expect_null(res)
})

test_that("WebGestalt_Pipe2 ORA combines up and down results and applies minEnrichmentRatio", {
  results <- data.frame(
    Gene = paste0("gene", 1:4),
    log2FoldChange = c(2, 2, -2, -2),
    padj = c(0.01, 0.01, 0.01, 0.01)
  )

  fake_webgestaltr <- function(interestGene, ...) {
    if (length(interestGene) == 0) return(data.frame())
    data.frame(
      geneSet = c("GO:1", "GO:2"),
      description = c("term1", "term2"),
      enrichmentRatio = c(2, 0.5),
      FDR = c(0.01, 0.02),
      stringsAsFactors = FALSE
    )
  }
  testthat::local_mocked_bindings(WebGestaltR = fake_webgestaltr, .package = "WebGestaltR")

  res <- WebGestalt_Pipe2(mode = "ORA", results = results, alpha = 0.05, FC = 1,
                           minEnrichmentRatio = 1)

  # GO:2 (enrichmentRatio 0.5) should be filtered out by minEnrichmentRatio = 1
  # in both the up and down (flipped to -0.5) sets, leaving only GO:1 twice
  # -> deduplicated to one row per unique geneSet.
  expect_true(all(abs(res$enrichmentRatio) >= 1))
})

test_that("WebGestalt_Pipe2 drops rows with NA padj or log2FoldChange before building gene sets", {
  results <- data.frame(
    Gene = paste0("gene", 1:5),
    log2FoldChange = c(2, 2, NA, -2, -2),
    padj = c(0.01, NA, 0.01, 0.01, 0.01)
  )

  captured_genes <- list()
  fake_webgestaltr <- function(interestGene, ...) {
    captured_genes[[length(captured_genes) + 1]] <<- interestGene
    if (length(interestGene) == 0) return(data.frame())
    data.frame(geneSet = "GO:1", description = "d", enrichmentRatio = 2, FDR = 0.01)
  }
  testthat::local_mocked_bindings(WebGestaltR = fake_webgestaltr, .package = "WebGestaltR")

  WebGestalt_Pipe2(mode = "ORA", results = results, alpha = 0.05, FC = 1)

  all_genes <- unlist(captured_genes)
  expect_false(any(is.na(all_genes)))
  expect_false("gene2" %in% all_genes)
  expect_false("gene3" %in% all_genes)
})

test_that("WebGestalt_Pipe2 GSEA selects gene ID and log2FoldChange by column name, not position", {
  # Regression test: the original code used results[c(7, 2)], which would
  # error outright on a data.frame with fewer than 7 columns, or silently
  # grab the wrong columns if their order differed from the assumed layout.
  results <- data.frame(
    Gene = paste0("gene", 1:5),
    padj = rep(0.01, 5),
    log2FoldChange = c(1, -1, 2, -2, 0.5)
  )

  captured_input <- NULL
  fake_webgestaltr <- function(interestGene, ...) {
    captured_input <<- interestGene
    data.frame(geneSet = "GO:1", description = "d", NES = 1.5, FDR = 0.01)
  }
  testthat::local_mocked_bindings(WebGestaltR = fake_webgestaltr, .package = "WebGestaltR")

  res <- WebGestalt_Pipe2(mode = "GSEA", results = results)

  expect_equal(colnames(captured_input), c("Gene", "log2FoldChange"))
  expect_equal(captured_input$log2FoldChange, results$log2FoldChange)
  expect_s3_class(res, "data.frame")
})

test_that("WebGestalt_Pipe2 GSEA messages and returns NULL when no terms are returned", {
  results <- data.frame(
    Gene = paste0("gene", 1:5),
    padj = rep(0.01, 5),
    log2FoldChange = c(1, -1, 2, -2, 0.5)
  )

  fake_webgestaltr <- function(interestGene, ...) data.frame()
  testthat::local_mocked_bindings(WebGestaltR = fake_webgestaltr, .package = "WebGestaltR")

  expect_message(
    res <- WebGestalt_Pipe2(mode = "GSEA", results = results),
    "No GSEA terms are returned"
  )
  expect_null(res)
})
