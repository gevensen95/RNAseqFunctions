#' WebGestalt Enrichment
#'
#' Runs GO/pathway enrichment (ORA or GSEA) on a DEG results data.frame
#' using WebGestaltR.
#'
#' @param mode "ORA" or "GSEA". Note: GSEA mode runs GSEA on the log2FC values.
#' @param results A [makeComp()] output data.frame. Any data frame with a
#' gene identifier column, "padj", and "log2FoldChange" also works — genes
#' are looked up by name via `geneID`, not by column position. Rows with NA
#' padj or log2FoldChange are dropped before running.
#' @param alpha pAdj cut-off for DEGs (ORA only)
#' @param FC absolute log2FC cut-off for DEGs (ORA only)
#' @param DB WebGestalt enrichDatabase to query against, see WebGestaltR::listGeneSet() for a full list.
#' @param species WebGestalt organism, see WebGestaltR::listOrganism() for a full list.
#' @param projectName Name of raw data folder. If `NULL` (the default), WebGestaltR runs without writing output files.
#' @param geneID Column for genes to use for enrichment
#' @param outputTerms The threshold for the top method. The default is 100.
#' @param minEnrichmentRatio Absolute enrichmentRatio cut-off applied after
#' combining up/down ORA results (ORA only). Terms with weaker enrichment
#' than this are dropped. Default 1 (no enrichment signal).
#' @return A dataframe of enriched terms (ORA or GSEA) with an "ID" column
#' ("geneSet: description") in column 1. Returns `NULL` (invisibly) with a
#' message if no terms are returned.
#' @export
WebGestalt_Pipe2 <- function(mode = "ORA", results, alpha = 0.05, FC = 0.5, DB = NULL,
                             species = "mmusculus", projectName = NULL,
                             geneID = 'Gene', outputTerms = 100,
                             minEnrichmentRatio = 1){

  if (is.null(DB)) {
    DB <- ifelse(mode == "ORA", "geneontology_Biological_Process_noRedundant", "pathway_KEGG")
  }

  results <- results[!is.na(results$padj) & !is.na(results$log2FoldChange), ]

  if (mode == "ORA") {

    run_ora <- function(genes, direction_suffix) {
      if (length(genes) == 0) {
        return(data.frame())
      }
      WebGestaltR::WebGestaltR(
        enrichMethod = "ORA", organism = species, enrichDatabase = DB,
        enrichDatabaseType = "genesymbol",
        isOutput = !is.null(projectName),
        projectName = if (!is.null(projectName)) paste0(projectName, direction_suffix) else NULL,
        interestGene = genes, interestGeneType = "genesymbol",
        sigMethod = "top", referenceGene = results[[geneID]],
        referenceGeneType = "genesymbol", topThr = outputTerms,
        hostName = "https://www.webgestalt.org"
      )
    }

    down_genes <- results[[geneID]][results$padj <= alpha & results$log2FoldChange < -1 * FC]
    up_genes   <- results[[geneID]][results$padj <= alpha & results$log2FoldChange > 1 * FC]

    ORA_down <- run_ora(down_genes, "_ORA_downreg")
    ORA_up   <- run_ora(up_genes, "_ORA_upreg")

    if (nrow(ORA_down) == 0 && nrow(ORA_up) == 0) {
      message("No ORA terms are returned!")
      return(invisible(NULL))
    }

    if (nrow(ORA_down) > 0) {
      ORA_down$enrichmentRatio <- ORA_down$enrichmentRatio * -1
    }

    ORA <- rbind(ORA_down, ORA_up)
    ORA <- ORA[abs(ORA$enrichmentRatio) >= minEnrichmentRatio, ]
    ORA <- ORA[order(abs(ORA$enrichmentRatio), decreasing = TRUE), ]
    ORA <- ORA[order(abs(ORA$FDR)), ]
    ORA <- ORA[!duplicated(ORA$geneSet), ]
    ORA$ID <- paste0(ORA$geneSet, ": ", ORA$description)
    ORA <- ORA[c(ncol(ORA), 1:(ncol(ORA) - 1))]

    return(ORA)

  } else if (mode == "GSEA") {

    gsea_input <- results[c(geneID, "log2FoldChange")]

    GSEA_output <- WebGestaltR::WebGestaltR(
      enrichMethod = "GSEA", organism = species, enrichDatabase = DB,
      isOutput = !is.null(projectName),
      projectName = if (!is.null(projectName)) paste0("GSEA_", projectName) else NULL,
      interestGene = gsea_input, interestGeneType = "genesymbol",
      sigMethod = "top", topThr = outputTerms, perNum = 1000,
      hostName = "https://www.webgestalt.org"
    )

    if (nrow(GSEA_output) == 0) {
      message("No GSEA terms are returned!")
      return(invisible(NULL))
    }

    GSEA_output$ID <- paste0(GSEA_output$geneSet, ": ", GSEA_output$description)
    GSEA_output <- GSEA_output[c(ncol(GSEA_output), 1:(ncol(GSEA_output) - 1))]

    return(GSEA_output)

  } else {
    warning("mode not recognized")
  }

}
