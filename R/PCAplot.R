#' Plots a sample PCA
#'
#' Normalizes raw counts (via DESeq2's variance-stabilizing transform or
#' edgeR's log-CPM) and plots a PCA of samples using the top variable genes,
#' for quick QC of sample clustering and batch effects.
#'
#' @param counts Matrix or data.frame of raw counts, genes as rows, samples
#' as columns (same orientation as [makeComp()]).
#' @param metadata Data.frame of sample metadata. `rownames(metadata)` must
#' match `colnames(counts)` exactly, in the same order.
#' @param method Which engine to use for normalization: `"DESeq2"`
#' (default, via `vst()`) or `"edgeR"` (via `cpm(log = TRUE)`).
#' @param color_by Optional metadata column name to color points by.
#' @param shape_by Optional metadata column name to shape points by.
#' @param ntop Number of most-variable genes to use for the PCA.
#' @return A ggplot2 PCA plot object.
#' @importFrom rlang .data
#' @export
PCAplot <- function(counts, metadata, method = c("DESeq2", "edgeR"),
                    color_by = NULL, shape_by = NULL, ntop = 500) {

  method <- match.arg(method)

  if (!identical(colnames(counts), rownames(metadata))) {
    stop("colnames(counts) must match rownames(metadata) exactly, ",
         "in the same order.")
  }

  norm_mat <- normalize_counts(counts, metadata, method)

  ntop <- min(ntop, nrow(norm_mat))
  row_var <- apply(norm_mat, 1, stats::var)
  select <- order(row_var, decreasing = TRUE)[seq_len(ntop)]

  pca <- stats::prcomp(t(norm_mat[select, , drop = FALSE]), scale. = FALSE)
  percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

  plot_df <- as.data.frame(pca$x[, 1:2])
  plot_df <- cbind(plot_df, metadata[rownames(plot_df), , drop = FALSE])

  plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]])) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    Ol_Reliable() +
    ggplot2::labs(
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance")
    )

  # Added conditionally (rather than inline in the initial aes()) so an
  # unused color_by/shape_by never evaluates to a literal NULL aesthetic.
  if (!is.null(color_by)) {
    plot <- plot + ggplot2::aes(colour = .data[[color_by]])
  }
  if (!is.null(shape_by)) {
    plot <- plot + ggplot2::aes(shape = .data[[shape_by]])
  }

  plot
}
