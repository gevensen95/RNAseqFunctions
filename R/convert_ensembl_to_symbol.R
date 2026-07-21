#' Convert Ensembl IDs to gene symbols
#'
#' Convenience wrapper around [convert_gene_names()] that looks up gene
#' symbols for a vector of Ensembl gene IDs.
#'
#' @param ensembl_ids Vector of Ensembl IDs
#' @param species Species ("mouse" or "human")
#' @param attributes Additional attributes to retrieve
#' @return A dataframe of Ensembl IDs, gene symbols, and any additional
#' attributes (if chosen)
#' @export
convert_ensembl_to_symbol <- function(ensembl_ids, species = 'mouse',
                                      attributes = NULL) {
  convert_gene_names(
    genes = ensembl_ids,
    species = species,
    gene_filters = "ensembl_gene_id",
    extra_attributes = attributes
  )
}
