test_that("convert_gene_names errors on unsupported species", {
  expect_error(
    convert_gene_names(genes = "ENSG00000141510", species = "zebrafish"),
    "Species not supported"
  )
})

test_that("convert_ensembl_to_symbol errors on unsupported species", {
  expect_error(
    convert_ensembl_to_symbol(ensembl_ids = "ENSMUSG00000064341", species = "zebrafish"),
    "Species not supported"
  )
})
