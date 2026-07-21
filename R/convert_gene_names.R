#' Convert Ensembl IDs to gene symbols
#'
#' This function is automatically converts human or mouse Ensembl IDs to
#'
#' @param genes Vector of genes
#' @param species Species
#' @param gene_filters Filters to query
#' @param attributes Additional attributes to retrieve
#' @return A dataframe of Ensembl IDs, gene symbols, and any additional
#' attributes (if chosen)
#' @export
convert_gene_names <- function(genes, species = 'mouse',
                               gene_filters = 'ensembl_gene_id',
                               extra_attributes = NULL) {
  # Choose the dataset based on species
  dataset <- switch(species,
                    "mouse" = "mmusculus_gene_ensembl",
                    "human" = "hsapiens_gene_ensembl",
                    stop("Species not supported. Use 'mouse' or 'human'."))

  # Connect to the Ensembl database
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)

  # Choose the appropriate symbol attribute
  species_attr <- if (species == "mouse") "mgi_symbol" else "hgnc_symbol"

  # Retrieve gene symbols
  gene_info <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", species_attr, extra_attributes),
    filters = gene_filters,
    values = genes,
    mart = ensembl
  )

  # Return the results
  return(gene_info)
}
