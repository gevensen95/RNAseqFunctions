#' Heatmap of top differentially expressed genes
#'
#' Normalizes raw counts, picks the top `n` genes by padj from a DEG results
#' data.frame (e.g. from [makeComp()]), and plots a z-scored expression
#' heatmap across samples.
#'
#' @param counts Matrix or data.frame of raw counts, genes as rows, samples
#' as columns (same orientation as [makeComp()]).
#' @param metadata Data.frame of sample metadata. `rownames(metadata)` must
#' match `colnames(counts)` exactly, in the same order.
#' @param results DEG results data.frame (e.g. from [makeComp()]) used to
#' pick and order the top genes by padj.
#' @param method Which engine to use for normalization: `"DESeq2"`
#' (default, via `vst()`) or `"edgeR"` (via `cpm(log = TRUE)`).
#' @param n Number of top genes (by padj) to include.
#' @param annotate_by Optional metadata column name to group and label
#' samples by (columns are reordered and split into facets by this
#' variable).
#' @param geneID Column in `results` containing gene identifiers matching
#' `rownames(counts)`.
#' @return A ggplot2 heatmap object.
#' @export
DEGheatmap <- function(counts, metadata, results, method = c("DESeq2", "edgeR"),
                       n = 30, annotate_by = NULL, geneID = "Gene") {

  method <- match.arg(method)

  if (!identical(colnames(counts), rownames(metadata))) {
    stop("colnames(counts) must match rownames(metadata) exactly, ",
         "in the same order.")
  }

  norm_mat <- normalize_counts(counts, metadata, method)

  ordered_results <- results[order(results$padj), ]
  top_genes <- utils::head(ordered_results[[geneID]], n)
  top_genes <- top_genes[top_genes %in% rownames(norm_mat)]

  if (length(top_genes) == 0) {
    stop("None of the top genes in `results[[geneID]]` were found in ",
         "rownames(counts). Check that `geneID` matches the gene ",
         "identifiers used as row names in `counts`.")
  }

  mat <- norm_mat[top_genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))

  plot_df <- as.data.frame(as.table(mat_z), stringsAsFactors = FALSE)
  names(plot_df) <- c("gene", "sample", "z")
  plot_df$gene <- factor(plot_df$gene, levels = rev(top_genes))

  if (!is.null(annotate_by)) {
    group_map <- stats::setNames(metadata[[annotate_by]], rownames(metadata))
    plot_df$group <- group_map[as.character(plot_df$sample)]
    sample_order <- rownames(metadata)[order(metadata[[annotate_by]])]
    plot_df$sample <- factor(plot_df$sample, levels = sample_order)
  }

  plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = gene, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    Ol_Reliable() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "z-score")

  if (!is.null(annotate_by)) {
    plot <- plot + ggplot2::facet_grid(~group, scales = "free_x", space = "free_x")
  }

  plot
}
